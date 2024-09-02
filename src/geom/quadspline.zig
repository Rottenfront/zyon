const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Point = mod.Point;

pub const QuadSpline = struct {
    data: std.ArrayList(Point),

    pub const ToQuadBez = struct {
        idx: usize,
        points: *std.ArrayList(Point),
    };
};
