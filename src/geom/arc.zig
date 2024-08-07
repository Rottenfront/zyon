const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;
const PathEl = mod.PathEl;
const Point = mod.Point;
const Vec2 = mod.Vec2;

/// A single arc segment.
pub const Arc = struct {
    /// The arc's centre point.
    center: Point,
    /// The arc's radii, where the vector's x-component is the radius in the
    /// positive x direction after applying `x_rotation`.
    radii: Vec2,
    /// The start angle in radians.
    start_angle: f64,
    /// The angle between the start and end of the arc, in radians.
    sweep_angle: f64,
    /// How much the arc is rotated, in radians.
    x_rotation: f64,

    /// Create a new `Arc`
    pub fn new(center: Point, radii: Vec2, start_angle: f64, sweep_angle: f64, x_rotation: f64) Arc {
        return Arc{
            .center = center,
            .radii = radii,
            .start_angle = start_angle,
            .sweep_angle = sweep_angle,
            .x_rotation = x_rotation,
        };
    }

    /// Create an iterator generating Bezier path elements.
    ///
    /// The generated elements can be appended to an existing bezier path.
    pub fn appendIter(self: *const Arc, tolerance: f64) ArcAppendIter {
        const sign = math.sign(self.sweep_angle);
        const scaled_err = @max(self.radii.x, self.radii.y) / tolerance;
        // Number of subdivisions per ellipse based on error tolerance.
        // Note: this may slightly underestimate the error for quadrants.
        const n_err = @max(math.pow(1.1163 * scaled_err, 1.0 / 6.0), 3.999_999);
        const n_f64 = math.ceil(n_err * @abs(self.sweep_angle) * (1.0 / (2.0 * math.pi)));
        const angle_step = self.sweep_angle / n_f64;
        const n: usize = math.round(n_f64);
        const arm_len = (4.0 / 3.0) * math.tan(@abs(0.25 * angle_step)) * sign;
        const angle0 = self.start_angle;
        const p0 = sample_ellipse(self.radii, self.x_rotation, angle0);

        return ArcAppendIter{
            .idx = 0,

            .center = self.center,
            .radii = self.radii,
            .x_rotation = self.x_rotation,
            .n = n,
            .arm_len = arm_len,
            .angle_step = angle_step,

            .p0 = p0,
            .angle0 = angle0,
        };
    }
};

pub const ArcAppendIter = struct {
    idx: usize,

    center: Point,
    radii: Vec2,
    x_rotation: f64,
    n: usize,
    arm_len: f64,
    angle_step: f64,

    p0: Vec2,
    angle0: f64,

    pub fn next(self: *ArcAppendIter) ?PathEl {
        if (self.idx >= self.n) {
            return null;
        }

        const angle1 = self.angle0 + self.angle_step;
        const p0 = self.p0;
        const p1 = p0 + sample_ellipse(self.radii, self.x_rotation, self.angle0 + math.pi / 2.0).mul(self.arm_len);
        const p3 = sample_ellipse(self.radii, self.x_rotation, angle1);
        const p2 = p3 - sample_ellipse(self.radii, self.x_rotation, self.angle1 + math.pi / 2.0).mul(self.arm_len);

        self.angle0 = angle1;
        self.p0 = p3;
        self.idx += 1;

        return PathEl{ .curve_to = struct { .p0 = self.center.sum(p1), .p1 = self.center.sum(p2), .end = self.center.sum(p3) } };
    }
};

/// Take the ellipse radii, how the radii are rotated, and the sweep angle, and return a point on
/// the ellipse.
fn sample_ellipse(radii: Vec2, x_rotation: f64, angle: f64) Vec2 {
    const angle_sin = math.sin(angle);
    const angle_cos = math.cos(angle);
    const u = radii.x * angle_cos;
    const v = radii.y * angle_sin;
    return rotate_pt(Vec2.new(u, v), x_rotation);
}

/// Rotate `pt` about the origin by `angle` radians.
fn rotate_pt(pt: Vec2, angle: f64) Vec2 {
    const angle_sin = math.sin(angle);
    const angle_cos = math.cos(angle);
    return Vec2.new(
        pt.x * angle_cos - pt.y * angle_sin,
        pt.x * angle_sin + pt.y * angle_cos,
    );
}
