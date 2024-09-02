const std = @import("std");
const math = std.math;
const mod = @import("../module.zig");
const Point = mod.Point;

pub const CubicBez = struct {
    p0: Point,
    p1: Point,
    p2: Point,
    p3: Point,

    pub fn new(p0: Point, p1: Point, p2: Point, p3: Point) @This() {
        return @This(){ .p0 = p0, .p1 = p1, .p2 = p2, .p3 = p3 };
    }

    /// Classification result for cusp detection.
    pub const CuspType = enum {
        /// Cusp is a loop.
        Loop,
        /// Cusp has two closely spaced inflection points.
        DoubleInflection,
    };
};
