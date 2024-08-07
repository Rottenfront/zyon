const Point = @import("point.zig").Point;

/// A single quadratic Bézier segment.
pub const QuadBez = struct {
    p0: Point,
    p1: Point,
    p2: Point,
};