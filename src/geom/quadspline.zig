const std = @import("std");
const Point = @import("point.zig").Point;

pub const QuadSpline = struct {
    data: std.ArrayList(Point),
};

pub const ToQuadBez = struct {
    idx: usize,
    points: *std.ArrayList(Point),
};