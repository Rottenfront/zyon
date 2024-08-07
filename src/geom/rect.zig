const math = @import("std").math;
const Affine = @import("affine.zig").Affine;
const Point = @import("point.zig").Point;

pub const Rect = struct {
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,

    pub fn new(x0: f64, y0: f64, x1: f64, y1: f64) Rect {
        return Rect{
            .x0 = x0,
            .y0 = y0,
            .x1 = x1,
            .y1 = y1,
        };
    }

    pub fn from_points(p0: Point, p1: Point) Rect {
        return Rect{
            .x0 = p0.x,
            .y0 = p0.y,
            .x1 = p1.x,
            .y1 = p1.y,
        };
    }

    pub fn width(self: Rect) f64 {
        return self.x1 - self.x0;
    }

    pub fn height(self: Rect) f64 {
        return self.y1 - self.y0;
    }

    pub fn rectUnion(self: *Rect, other: Rect) Rect {
        return Rect.new(
            self.x0.min(other.x0),
            self.y0.min(other.y0),
            self.x1.max(other.x1),
            self.y1.max(other.y1),
        );
    }
};
