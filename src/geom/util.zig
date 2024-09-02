const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;
const Point = mod.Point;
const Vec2 = mod.Vec2;

pub const ArrayAccessError = error{OutOfBounds};

pub const util = struct {
    /// An approximation to $\int (1 + 4x^2) ^ -0.25 dx$
    ///
    /// This is used for flattening curves.
    pub fn approx_parabola_integral(x: f64) f64 {
        const D: f64 = 0.67;
        return x / (1.0 - D + @sqrt(@sqrt(math.powi(D, 4) + 0.25 * x * x)));
    }

    /// An approximation to the inverse parabola integral.
    pub fn approx_parabola_inv_integral(x: f64) f64 {
        const B: f64 = 0.39;
        return x * (1.0 - B + @sqrt(B * B + 0.25 * x * x));
    }

    /// Take the ellipse radii, how the radii are rotated, and the sweep angle, and return a point on
    /// the ellipse.
    pub fn sample_ellipse(radii: Vec2, x_rotation: f64, angle: f64) Vec2 {
        const angle_sin = math.sin(angle);
        const angle_cos = math.cos(angle);
        const u = radii.x * angle_cos;
        const v = radii.y * angle_sin;
        return rotate_pt(Vec2.new(u, v), x_rotation);
    }

    /// Rotate `pt` about the origin by `angle` radians.
    pub fn rotate_pt(pt: Vec2, angle: f64) Vec2 {
        const angle_sin = math.sin(angle);
        const angle_cos = math.cos(angle);
        return Vec2.new(
            pt.x * angle_cos - pt.y * angle_sin,
            pt.x * angle_sin + pt.y * angle_cos,
        );
    }

    pub fn expect_near(p0: Point, p1: Point) !void {
        try std.testing.expect(p1.distance(p0) < 1e-9);
    }

    pub fn affine_expect_near(a0: Affine, a1: Affine) !void {
        for (0..6) |i| {
            try std.testing.expect(@abs(a0.a[i] - a1.a[i]) < 1e-9);
        }
    }
};
