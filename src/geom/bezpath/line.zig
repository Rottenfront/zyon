const std = @import("std");
const math = std.math;
const mod = @import("../module.zig");
const Point = mod.Point;

pub const Line = struct {
    start: Point,
    end: Point,

    pub fn new(start: Point, end: Point) Line {
        return Line{ .start = start, .end = end };
    }
};
