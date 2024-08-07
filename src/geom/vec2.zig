const std = @import("std");
const Affine = @import("affine.zig").Affine;
const Point = @import("point.zig").Point;
const Size = @import("size.zig").Size;

/// A 2D vector.
///
/// This is intended primarily for a vector in the mathematical sense,
/// but it can be interpreted as a translation, and converted to and
/// from a point (vector relative to the origin) and size.
pub const Vec2 = struct {
    /// The x-coordinate.
    x: f64,
    /// The y-coordinate.
    y: f64,

    /// The vector (0, 0).
    const ZERO: Vec2 = Vec2{ .x = 0.0, .y = 0.0 };

    /// Create a new vector.
    pub fn new(x: f64, y: f64) Vec2 {
        return Vec2{
            .x = x,
            .y = y,
        };
    }

    /// Convert this vector into a `Point`.
    pub fn toPoint(self: Vec2) Point {
        return Point{
            .x = self.x,
            .y = self.y,
        };
    }

    /// Convert this vector into a `Size`.
    pub fn toSize(self: Vec2) Size {
        return Size{
            .width = self.x,
            .height = self.y,
        };
    }

    /// Create a `Vec2` with the same value for x and y
    pub fn splat(v: f64) Vec2 {
        return Vec2{
            .x = v,
            .y = v,
        };
    }

    /// Dot product of two vectors.
    pub fn dot(self: Vec2, other: Vec2) f64 {
        return self.x * other.x + self.y * other.y;
    }

    /// Cross product of two vectors.
    ///
    /// This is signed so that (0, 1) × (1, 0) = 1.
    pub fn cross(self: Vec2, other: Vec2) f64 {
        return self.x * other.x - self.y * other.y;
    }

    /// Magnitude of vector.
    ///
    /// This is similar to `@sqrt(self.hypot2())` but defers to the platform `hypot` method, which
    /// in general will handle the case where `self.hypot2() > std.math.inf(f64)`.
    pub fn hypot(self: Vec2) f64 {
        return std.math.hypot(self.x, self.y);
    }

    /// Magnitude of vector.
    ///
    /// This is an alias for [`Vec2.hypot`].
    pub fn length(self: Vec2) f64 {
        return self.hypot();
    }

    /// Magnitude squared of vector.
    pub fn hypot2(self: Vec2) f64 {
        return self.dot(self);
    }

    /// Magnitude squared of vector.
    ///
    /// This is an alias for [`Vec2.hypot2`].
    pub fn length_squared(self: Vec2) f64 {
        return self.hypot2();
    }

    /// Find the angle in radians between this vector and the vector `Vec2 { .x = 1.0, .y = 0.0 }`
    /// in the positive `y` direction.
    ///
    /// If the vector is interpreted as a complex number, this is the argument.
    /// The angle is expressed in radians.
    pub fn atan2(self: Vec2) f64 {
        return std.math.atan2(self.y, self.x);
    }

    /// Find the angle in radians between this vector and the vector `Vec2 { .x = 1.0, .y = 0.0 }`
    /// in the positive `y` direction.
    ///
    /// This is an alias for [`Vec2.atan2`]
    pub fn angle(self: Vec2) f64 {
        return self.atan2();
    }

    /// A unit vector of the given angle.
    ///
    /// With `th` at zero, the result is the positive X unit vector, and
    /// at π/2, it is the positive Y unit vector. The angle is expressed
    /// in radians.
    ///
    /// Thus, in a Y-down coordinate system (as is common for graphics),
    /// it is a clockwise rotation, and in Y-up (traditional for math), it
    /// is anti-clockwise. This convention is consistent with
    /// [`Affine.rotate`].
    pub fn from_angle(th: f64) Vec2 {
        return Vec2{
            .x = std.math.cos(th),
            .y = std.math.sin(th),
        };
    }

    /// Linearly interpolate between two vectors.
    pub fn lerp(self: Vec2, other: Vec2, t: f64) Vec2 {
        return self.sum(other.mul(t));
    }

    /// Returns a vector of magnitude 1.0 with the same angle as `self`; i.e.
    /// a unit/direction vector.
    ///
    /// This produces `NaN` values when the magnitude is `0`.
    pub fn normalize(self: Vec2) Vec2 {
        return self.div(self.hypot());
    }

    /// Returns a new `Vec2`,
    /// with `x` and `y` rounded to the nearest integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Vec2.new(3.3, 3.6).round();
    /// const b = Vec2.new(3.0, -3.1).round();
    /// try expectEqual(a.x, 3.0);
    /// try expectEqual(a.y, 4.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -3.0);
    /// ```
    pub fn round(self: Vec2) Vec2 {
        return Vec2{
            .x = std.math.round(self.x),
            .y = std.math.round(self.y),
        };
    }

    /// Returns a new `Vec2`,
    /// with `x` and `y` rounded up to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Vec2.new(3.3, 3.6).ceil();
    /// const b = Vec2.new(3.0, -3.1).ceil();
    /// try expectEqual(a.x, 4.0);
    /// try expectEqual(a.y, 4.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -3.0);
    pub fn ceil(self: Vec2) Vec2 {
        return Vec2{
            .x = std.math.ceil(self.x),
            .y = std.math.ceil(self.y),
        };
    }

    /// Returns a new `Vec2`,
    /// with `x` and `y` rounded down to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Vec2.new(3.3, 3.6).floor();
    /// const b = Vec2.new(3.0, -3.1).floor();
    /// try expectEqual(a.x, 3.0);
    /// try expectEqual(a.y, 3.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -4.0);
    /// ```
    pub fn floor(self: Vec2) Vec2 {
        return Vec2{
            .x = std.math.floor(self.x),
            .y = std.math.floor(self.y),
        };
    }

    /// Returns a new `Vec2`,
    /// with `x` and `y` rounded towards zero to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const a = Vec2.new(3.3, 3.6).floor();
    /// const b = Vec2.new(3.0, -3.1).floor();
    /// try expectEqual(a.x, 3.0);
    /// try expectEqual(a.y, 3.0);
    /// try expectEqual(b.x, 3.0);
    /// try expectEqual(b.y, -3.0);
    pub fn trunc(self: Vec2) Vec2 {
        return Vec2{
            .x = std.math.trunc(self.x),
            .y = std.math.trunc(self.y),
        };
    }

    /// Is this Vec2 finite?
    pub fn isFinite(self: *const Vec2) bool {
        return std.math.isFinite(self.x) || std.math.isFinite(self.y);
    }

    /// Is this Vec2 NaN?
    pub fn isNan(self: *const Vec2) bool {
        return std.math.isNan(self.x) || std.math.isNan(self.y);
    }

    pub fn neg(self: Vec2) Vec2 {
        return Vec2{
            .x = -self.x,
            .y = -self.y,
        };
    }

    /// Returns sum of vectors
    pub fn sum(self: Vec2, other: Vec2) Vec2 {
        return Vec2{
            .x = self.x + other.x,
            .y = self.y + other.y,
        };
    }

    /// Returns difference of vectors
    pub fn sub(self: Vec2, other: Vec2) Vec2 {
        return Vec2{
            .x = self.x - other.x,
            .y = self.y - other.y,
        };
    }

    pub fn mul(self: Vec2, t: f64) Vec2 {
        return Vec2{
            .x = self.x * t,
            .y = self.y * t,
        };
    }

    pub fn div(self: Vec2, t: f64) Vec2 {
        return Vec2{
            .x = self.x / t,
            .y = self.y / t,
        };
    }

    /// Apply an affine to a vector
    pub fn applyAffine(self: Vec2, transform: Affine) Vec2 {
        return Point{
            .x = self.x * transform.a[0] + self.y * transform.a[2] + transform.a[4],
            .y = self.x * transform.a[1] + self.y * transform.a[3] + transform.a[5],
        };
    }
};
