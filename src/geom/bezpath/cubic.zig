const std = @import("std");
const math = std.math;
const mod = @import("../module.zig");
const Point = mod.Point;

pub const CubicBez = struct {
    start: Point,
    p0: Point,
    p1: Point,
    end: Point,

    pub fn new(start: Point, p0: Point, p1: Point, end: Point) CubicBez {
        return CubicBez{ .start = start, .p0 = p0, .p1 = p1, .end = end };
    }
};

pub const CuspType = enum {
    /// Cusp is a loop.
    Loop,
    /// Cusp has two closely spaced inflection points.
    DoubleInflection,
};
