const std = @import("std");
const math = std.math;
const Point = @import("point.zig").Point;

pub const CubicBez = struct {
    p0: Point,
    p1: Point,
    p2: Point,
    p3: Point,
};

pub const CuspType = enum {
    /// Cusp is a loop.
    Loop,
    /// Cusp has two closely spaced inflection points.
    DoubleInflection,
};