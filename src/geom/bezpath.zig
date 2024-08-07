const std = @import("std");
const Point = @import("point.zig").Point;
const util = @import("util.zig");

pub const BezPath = struct {
    data: std.ArrayList(PathEl),
};

pub const PathEl = union(enum) {
    move_to: Point,
    line_to: Point,
    quad_to: struct { p0: Point, end: Point },
    curve_to: struct { p0: Point, p1: Point, end: Point },
    close_path,
};