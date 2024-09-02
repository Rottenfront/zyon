const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;
const Size = mod.Size;
const Point = mod.Point;
const TranslateScale = mod.TranslateScale;

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
    pub const ZERO: @This() = @This(){ .x = 0.0, .y = 0.0 };

    /// Create a new vector.
    pub fn new(x: f64, y: f64) @This() {
        return @This(){
            .x = x,
            .y = y,
        };
    }

    /// Convert this vector into a `Point`.
    pub fn toPoint(self: @This()) Point {
        return Point{
            .x = self.x,
            .y = self.y,
        };
    }

    /// Convert this vector into a `Size`.
    pub fn toSize(self: @This()) Size {
        return Size{
            .width = self.x,
            .height = self.y,
        };
    }

    /// Create a `Vec2` with the same value for x and y
    pub fn splat(v: f64) @This() {
        return @This(){
            .x = v,
            .y = v,
        };
    }

    /// Dot product of two vectors.
    pub fn dot(self: @This(), other: @This()) f64 {
        return self.x * other.x + self.y * other.y;
    }

    /// Cross product of two vectors.
    ///
    /// This is signed so that (0, 1) × (1, 0) = 1.
    pub fn cross(self: @This(), other: @This()) f64 {
        return self.x * other.x - self.y * other.y;
    }

    /// Magnitude of vector.
    ///
    /// This is similar to `@sqrt(self.hypot2())` but defers to the platform `hypot` method, which
    /// in general will handle the case where `self.hypot2() > std.math.inf(f64)`.
    pub fn hypot(self: @This()) f64 {
        return std.math.hypot(self.x, self.y);
    }

    /// Magnitude of vector.
    ///
    /// This is an alias for [`Vec2.hypot`].
    pub fn length(self: @This()) f64 {
        return self.hypot();
    }

    /// Magnitude squared of vector.
    pub fn hypot2(self: @This()) f64 {
        return self.dot(self);
    }

    /// Magnitude squared of vector.
    ///
    /// This is an alias for [`Vec2.hypot2`].
    pub fn length_squared(self: @This()) f64 {
        return self.hypot2();
    }

    /// Find the angle in radians between this vector and the vector `Vec2 { .x = 1.0, .y = 0.0 }`
    /// in the positive `y` direction.
    ///
    /// If the vector is interpreted as a complex number, this is the argument.
    /// The angle is expressed in radians.
    pub fn atan2(self: @This()) f64 {
        return std.math.atan2(self.y, self.x);
    }

    /// Find the angle in radians between this vector and the vector `Vec2 { .x = 1.0, .y = 0.0 }`
    /// in the positive `y` direction.
    ///
    /// This is an alias for [`Vec2.atan2`]
    pub fn angle(self: @This()) f64 {
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
    pub fn from_angle(th: f64) @This() {
        return @This(){
            .x = std.math.cos(th),
            .y = std.math.sin(th),
        };
    }

    /// Linearly interpolate between two vectors.
    pub fn lerp(self: @This(), other: @This(), t: f64) @This() {
        return self.sum(other.mul(t));
    }

    /// Returns a vector of magnitude 1.0 with the same angle as `self`; i.e.
    /// a unit/direction vector.
    ///
    /// This produces `NaN` values when the magnitude is `0`.
    pub fn normalize(self: @This()) @This() {
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
    pub fn round(self: @This()) @This() {
        return @This(){
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
    pub fn ceil(self: @This()) @This() {
        return @This(){
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
    pub fn floor(self: @This()) @This() {
        return @This(){
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
    pub fn trunc(self: @This()) @This() {
        return @This(){
            .x = std.math.trunc(self.x),
            .y = std.math.trunc(self.y),
        };
    }

    /// Is this Vec2 finite?
    pub fn isFinite(self: *const @This()) bool {
        return std.math.isFinite(self.x) || std.math.isFinite(self.y);
    }

    /// Is this Vec2 NaN?
    pub fn isNan(self: *const @This()) bool {
        return std.math.isNan(self.x) || std.math.isNan(self.y);
    }

    /// Negative (opposite) vector
    pub fn neg(self: @This()) @This() {
        return @This(){
            .x = -self.x,
            .y = -self.y,
        };
    }

    /// Sum of vectors
    pub fn sum(self: @This(), other: @This()) @This() {
        return @This(){
            .x = self.x + other.x,
            .y = self.y + other.y,
        };
    }

    /// Difference of vectors
    pub fn sub(self: @This(), other: @This()) @This() {
        return @This(){
            .x = self.x - other.x,
            .y = self.y - other.y,
        };
    }

    /// Multiplied by `t` value vector
    pub fn mul(self: @This(), t: f64) @This() {
        return @This(){
            .x = self.x * t,
            .y = self.y * t,
        };
    }

    /// Divided by `t` value vector
    pub fn div(self: @This(), t: f64) @This() {
        return @This(){
            .x = self.x / t,
            .y = self.y / t,
        };
    }

    /// Apply an affine to a vector
    pub fn applyAffine(self: @This(), transform: Affine) @This() {
        return Point{
            .x = self.x * transform.a[0] + self.y * transform.a[2] + transform.a[4],
            .y = self.x * transform.a[1] + self.y * transform.a[3] + transform.a[5],
        };
    }

    /// Apply an affine to a point
    pub fn applyTranslateScale(self: @This(), transform: TranslateScale) @This() {
        return @This(){
            .x = self.x * transform.scale + transform.translation.x,
            .y = self.y * transform.scale + transform.translation.y,
        };
    }
};
