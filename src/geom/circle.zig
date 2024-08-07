const std = @import("std");
const math = std.math;
const Point = @import("point.zig").Point;

pub const Circle = struct {
    center: Point,
    radius: f64,
};