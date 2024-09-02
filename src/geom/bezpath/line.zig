const std = @import("std");
const math = std.math;
const mod = @import("../module.zig");
const Point = mod.Point;

pub const Line = struct {
    p0: Point,
    p1: Point,

    pub fn new(p0: Point, p1: Point) @This() {
        return @This(){ .p0 = p0, .p1 = p1 };
    }
};
