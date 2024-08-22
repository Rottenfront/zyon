const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;
const Size = mod.Size;
const TranslateScale = mod.TranslateScale;
const Vec2 = mod.Vec2;

/// A 2D point.
///
/// This type represents a point in 2D space. It has the same layout as `Vec2`, but
/// its meaning is different: `Vec2` represents a change in location (for example velocity).
///
/// In general, library overloads math operators where it makes sense, for example implementing
/// `Affine * Point` as the point under the affine transformation. However `Point.sum(Point)` and
/// `Point.mul(f64)` are not implemented, because the operations do not make geometric sense. If you
/// need to apply these operations, then 1) check what you're doing makes geometric sense, then 2)
/// use [`Point.toVec2`] to convert the point to a `Vec2`.
pub const Point = struct {
    /// The x-coordinate.
    x: f64,
    /// The y-coordinate.
    y: f64,

    /// The point (0, 0).
    pub const ZERO: Point = Point{ .x = 0.0, .y = 0.0 };

    /// The point at the origin; (0, 0).
    pub const ORIGIN: Point = Point{ .x = 0.0, .y = 0.0 };

    /// Create a new `Point` with the provided `x` and `y` coordinates.
    pub fn new(x: f64, y: f64) Point {
        return Point{ .x = x, .y = y };
    }

    /// Convert this point into a `Vec2`.
    pub fn toVec2(self: Point) Vec2 {
        return Vec2.new(self.x, self.y);
    }

    /// Linearly interpolate between two points.
    pub fn lerp(self: Point, other: Point, t: f64) Point {
        return self.toVec2().lerp(other.toVec2(), t).toPoint();
    }

    /// Determine the midpoint of two points.
    pub fn midpoint(self: Point, other: Point) Point {
        return Point.new(0.5 * (self.x + other.x), 0.5 * (self.y + other.y));
    }

    /// Euclidean distance.
    pub fn distance(self: Point, other: Point) f64 {
        return self.toVec2().sub(other.toVec2()).hypot();
    }

    /// Squared Euclidean distance.
    pub fn distanceSquared(self: Point, other: Point) f64 {
        return self.toVec2().sub(other.toVec2()).hypot2();
    }

    /// Returns a new `Point`,
    /// with `x` and `y` rounded to the nearest integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Point::new(3.3, 3.6).round();
    /// const b = Point::new(3.0, -3.1).round();
    /// try expectEqual(a.x, 3.0);
    /// try expectEqual(a.y, 4.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -3.0);
    /// ```
    pub fn round(self: Point) Point {
        return Point.new(self.x.round(), self.y.round());
    }

    /// Returns a new `Point`,
    /// with `x` and `y` rounded up to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Point::new(3.3, 3.6).ceil();
    /// const b = Point::new(3.0, -3.1).ceil();
    /// try expectEqual(a.x, 4.0);
    /// try expectEqual(a.y, 4.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -3.0);
    /// ```
    pub fn ceil(self: Point) Point {
        return Point.new(self.x.ceil(), self.y.ceil());
    }

    /// Returns a new `Point`,
    /// with `x` and `y` rounded down to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Point::new(3.3, 3.6).floor();
    /// const b = Point::new(3.0, -3.1).floor();
    /// try expectEqual(a.x, 3.0);
    /// try expectEqual(a.y, 3.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -4.0);
    /// ```
    pub fn floor(self: Point) Point {
        return Point.new(self.x.floor(), self.y.floor());
    }

    /// Returns a new `Point`,
    /// with `x` and `y` rounded towards zero to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Point::new(3.3, 3.6).trunc();
    /// const b = Point::new(3.0, -3.1).trunc();
    /// try expectEqual(a.x, 3.0);
    /// try expectEqual(a.y, 3.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -3.0);
    /// ```
    pub fn trunc(self: Point) Point {
        return Point.new(math.trunc(self.y), math.trunc(self.y));
    }

    /// Is this point finite?
    pub fn isFinite(self: *const Point) bool {
        return math.isFinite(self.x) and math.isFinite(self.y);
    }

    /// Is this point NaN?
    pub fn isNan(self: *const Point) bool {
        return math.isNan(self.x) or math.isNan(self.y);
    }

    pub fn sum(self: Point, other: Point) Point {
        return Point{ .x = self.x + other.x, .y = self.y + other.y };
    }


    /// Apply an affine to a point
    pub fn applyAffine(self: Point, transform: Affine) Point {
        return Point{
            .x = self.x * transform.a[0] + self.y * transform.a[2] + transform.a[4],
            .y = self.x * transform.a[1] + self.y * transform.a[3] + transform.a[5],
        };
    }

    /// Apply an affine to a point
    pub fn applyTranslateScale(self: Point, transform: TranslateScale) Point {
        return Point{
            .x = self.x * transform.scale + transform.translation.x,
            .y = self.y * transform.scale + transform.translation.y,
        };
    }
};
