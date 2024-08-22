const std = @import("std");
const math = std.math;
const mod = @import("../module.zig");
const Point = mod.Point;

/// A single quadratic Bézier segment.
pub const QuadBez = struct {
    start: Point,
    p0: Point,
    end: Point,

    pub fn new(start: Point, p0: Point, end: Point) QuadBez {
        return QuadBez{ .start = start, .p0 = p0, .end = end };
    }
};
