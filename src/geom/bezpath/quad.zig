const std = @import("std");
const math = std.math;
const mod = @import("../module.zig");
const Point = mod.Point;
const CubicBez = mod.CubicBez;
const util = mod.util;

/// A single quadratic Bézier segment.
pub const QuadBez = struct {
    p0: Point,
    p1: Point,
    p2: Point,

    pub fn new(p0: Point, p1: Point, p2: Point) @This() {
        return @This(){ .p0 = p0, .p1 = p1, .p2 = p2 };
    }

    /// Raise the order by 1.
    ///
    /// Returns a cubic Bézier segment that exactly represents this quadratic.
    pub fn raise(self: @This()) CubicBez {
        return CubicBez.new(
            self.p0,
            self.p0.sum(self.p1.sub(self.p0).mul(2.0 / 3.0)),
            self.p2.sum(self.p1.sub(self.p2).mul(2.0 / 3.0)),
            self.p2,
        );
    }

    /// Estimate the number of subdivisions for flattening.
    pub fn estimateSubdiv(self: *const @This(), sqrt_tol: f64) FlattenParams {
        // Determine transformation to $y = x^2$ parabola.
        const d01 = self.p1.sub(self.p0).toVec2();
        const d12 = self.p2.sub(self.p1).toVec2();
        const dd = d01.sub(d12);
        const cross = self.p2.sub(self.p0).toVec2().cross(dd);
        const x0 = d01.dot(dd) / cross;
        const x2 = d12.dot(dd) / cross;
        const scale = @abs(cross / (dd.hypot() * (x2 - x0)));

        // Compute number of subdivisions needed.
        const a0 = util.approx_parabola_integral(x0);
        const a2 = util.approx_parabola_integral(x2);
        const val = val: {
            if (math.isFinite(scale)) {
                const da = @abs(a2 - a0);
                const sqrt_scale = @sqrt(scale);
                if (x0 * x2 > 0.0) {
                    break :val da * sqrt_scale;
                } else {
                    const xmin = sqrt_tol / sqrt_scale;
                    break :val sqrt_tol * da / util.approx_parabola_integral(xmin);
                }
            }
            break :val 0.0;
        };
        const k0 = util.approx_parabola_inv_integral(a0);
        const k2 = util.approx_parabola_inv_integral(a2);

        const uscale = 1 / (k2 - k0);

        return FlattenParams{
            .a0 = a0,
            .a2 = a2,
            .k0 = k0,
            .uscale = uscale,
            .val = val,
        };
    }

    // Maps a value from 0..1 to 0..1.
    pub fn determineSubdivT(
        params: *const FlattenParams,
        x: f64,
    ) f64 {
        const a = params.a0 + (params.a2 - params.a0) * x;
        const u = util.approx_parabola_inv_integral(a);
        return (u - params.k0) * params.uscale;
    }

    /// Is this quadratic Bezier curve finite?
    pub fn isFinite(self: *const @This()) bool {
        return self.p0.isFinite() and self.p1.isFinite() and self.p2.isFinite();
    }

    /// Is this quadratic Bezier curve NaN?
    pub fn isNan(self: *const @This()) bool {
        return self.p0.isNan() or self.p1.isNan() or self.p2.isNan();
    }

    pub fn eval(self: *const @This(), t: f64) Point {
        const mt = 1.0 - t;
        const v0 = self.p0.toVec2().mul(mt * mt);
        const v2 = self.p2.toVec2().mul(t);
        const v1 = self.p1.toVec2().mul(mt * 2.0).sum(v2).mul(t);
        return v0.sum(v1).toPoint();
    }

    pub fn subsegment(self: *const @This(), t0: f64, t1: f64) QuadBez {
        const p0 = self.eval(t0);
        const p2 = self.eval(t1);
        const p1 = p0.sum(self.p1.sub(self.p0).lerp(self.p2.sub(self.p1), t0).mul(t1 - t0));
        return QuadBez{ .p0 = p0, .p1 = p1, .p2 = p2 };
    }

    pub fn subdivide(self: *const @This()) struct { QuadBez, QuadBez } {
        const pm = self.eval(0.5);
        return .{
            QuadBez.new(self.p0, self.p0.midpoint(self.p1), pm),
            QuadBez.new(pm, self.p1.midpoint(self.p2), self.p2),
        };
    }

    pub fn start(self: *const @This()) Point {
        return self.p0;
    }

    pub fn end(self: *const @This()) Point {
        return self.p2;
    }

    /// An iterator for quadratic beziers.
    pub const QuadBezIter = struct {
        quad: QuadBez,
        ix: usize,
    };

    pub const FlattenParams = struct {
        a0: f64,
        a2: f64,
        k0: f64,
        uscale: f64,
        val: f64,
    };
};
