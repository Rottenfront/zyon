const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;
const BezPath = mod.BezPath;
const util = mod.util;
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
    pub fn new(center: Point, radii: Vec2, start_angle: f64, sweep_angle: f64, x_rotation: f64) @This() {
        return @This(){
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
    pub fn appendIter(self: *const @This(), tolerance: f64) AppendIter {
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
        const p0 = util.sample_ellipse(self.radii, self.x_rotation, angle0);

        return AppendIter{
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

    /// Converts an Arc into a list of cubic bezier segments.
    pub fn to_cubic_beziers(
        self: @This(),
        tolerance: f64,
        allocator: ?std.mem.Allocator,
    ) !std.ArrayList(BezPath.Element) {
        var path = self.appendIter(tolerance);
        var result = std.ArrayList(BezPath.Element).init(alloc: {
            if (allocator != null) {
                break :alloc allocator;
            }
            break :alloc std.heap.page_allocator;
        });

        var current_elem = path.next();
        while (current_elem != null) {
            result.append(current_elem);
            current_elem = path.next();
        }

        return result;
    }

    pub const AppendIter = struct {
        idx: isize,

        center: Point,
        radii: Vec2,
        x_rotation: f64,
        n: usize,
        arm_len: f64,
        angle_step: f64,

        p0: Vec2,
        angle0: f64,

        pub fn next(self: *@This()) ?BezPath.Element {
            if (self.idx == -1) {
                self.idx += 1;
                return BezPath.Element{
                    .move_to = struct {
                        .p = self.p0,
                    },
                };
            }

            if (self.idx >= self.n) {
                return null;
            }

            const angle1 = self.angle0 + self.angle_step;
            const p0 = self.p0;
            const p1 = p0 + util.sample_ellipse(
                self.radii,
                self.x_rotation,
                self.angle0 + math.pi / 2.0,
            ).mul(self.arm_len);
            const p3 = util.sample_ellipse(self.radii, self.x_rotation, angle1);
            const p2 = p3 - util.sample_ellipse(
                self.radii,
                self.x_rotation,
                self.angle1 + math.pi / 2.0,
            ).mul(self.arm_len);

            self.angle0 = angle1;
            self.p0 = p3;
            self.idx += 1;

            return BezPath.Element{ .curve_to = struct {
                .p0 = self.center.sum(p1),
                .p1 = self.center.sum(p2),
                .end = self.center.sum(p3),
            } };
        }
    };

    pub fn pathElements(self: *const @This(), tolerance: f64) AppendIter {
        var iter = self.appendIter(tolerance);
        iter.idx = -1;
        return iter;
    }

    /// Note: shape isn't closed so area is not well defined.
    pub fn area(self: *const @This()) f64 {
        return math.pi * self.radii.x * self.radii.y;
    }

    /// The perimeter of the arc.
    ///
    /// For now we just approximate by using the bezier curve representation.
    pub fn perimeter(self: *const @This(), accuracy: f64) f64 {
        return self.pathElements(0.1).perimeter(accuracy);
    }
};
