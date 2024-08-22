const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Point = mod.Point;

pub const Circle = struct {
    center: Point,
    radius: f64,
};
