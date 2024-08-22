const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Rect = mod.Rect;
const RoundRect = mod.RoundRect;
const RoundRectRadii = mod.RoundRectRadii;
const Vec2 = mod.Vec2;

/// A 2D size.
pub const Size = struct {
    /// The width.
    width: f64,
    /// The height.
    height: f64,

    /// A size with zero width or height
    pub const ZERO: Size = Size.new(0.0, 0.0);

    /// Create a new `Size` with the provided `width` and `height`.
    pub fn new(width: f64, height: f64) Size {
        return Size{ .width = width, .height = height };
    }

    /// Returns the max of `width` and `height`.
    ///
    /// # Examples
    ///
    /// ```
    /// const size = Size.new(-10.5, 42.0);
    /// try expectEqual(size.max_side(), 42.0);
    /// ```
    pub fn max_side(self: *const Size) f64 {
        return @max(self.width, self.height);
    }

    /// Returns the min of `width` and `height`.
    ///
    /// # Examples
    ///
    /// ```
    /// const size = Size.new(-10.5, 42.0);
    /// try expectEqual(size.max_side(), -10.5);
    /// ```
    pub fn min_size(self: *const Size) f64 {
        return @min(self.width, self.height);
    }

    /// The area covered by this size.
    pub fn area(self: *const Size) f64 {
        return self.width * self.height;
    }

    /// Whether this size has zero area.
    ///
    /// Note: a size with negative area is not considered empty.
    pub fn isEmpty(self: *const Size) f64 {
        return self.area() == 0.0;
    }

    /// Returns a new size bounded by `min` and `max.`
    ///
    /// # Examples
    ///
    /// ```
    /// use kurbo::Size;
    ///
    /// const this = Size.new(0.0, 100.0);
    /// const min = Size.new(10.0, 10.0);
    /// const max = Size.new(50.0, 50.0);
    /// try expectEqual(this.clamp(min, max), Size.new(10.0, 50.0))
    /// ```
    pub fn clamp(self: Size, min: Size, max: Size) Size {
        const width = @min(@max(self.width, min.width), max.width);
        const height = @min(@max(self.height, min.height), max.height);
        return Size.new(width, height);
    }

    /// Convert this size into a [`Vec2`], with `width` mapped to `x` and `height`
    /// mapped to `y`.
    pub fn toVec2(self: Size) Size {
        return Vec2.new(self.width, self.height);
    }

    /// Returns a new `Size`,
    /// with `width` and `height` rounded to the nearest integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const size_pos = Size.new(3.3, 3.6).round();
    /// try expectEqual(size_pos.width, 3.0);
    /// try expectEqual(size_pos.height, 4.0);
    /// const size_neg = Size.new(-3.3, -3.6).round();
    /// try expectEqual(size_neg.width, -3.0);
    /// try expectEqual(size_neg.height, -4.0);
    /// ```
    pub fn round(self: Size) Size {
        return Size.new(math.round(self.width), math.round(self.height));
    }

    /// Returns a new `Size`,
    /// with `width` and `height` rounded up to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const size_pos = Size.new(3.3, 3.6).ceil();
    /// try expectEqual(size_pos.width, 4.0);
    /// try expectEqual(size_pos.height, 4.0);
    /// const size_neg = Size.new(-3.3, -3.6).ceil();
    /// try expectEqual(size_neg.width, -3.0);
    /// try expectEqual(size_neg.height, -3.0);
    /// ```
    pub fn ceil(self: Size) Size {
        return Size.new(math.ceil(self.width), math.ceil(self.height));
    }

    /// Returns a new `Size`,
    /// with `width` and `height` rounded down to the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const size_pos = Size.new(3.3, 3.6).floor();
    /// try expectEqual(size_pos.width, 3.0);
    /// try expectEqual(size_pos.height, 3.0);
    /// const size_neg = Size.new(-3.3, -3.6).floor();
    /// try expectEqual(size_neg.width, -4.0);
    /// try expectEqual(size_neg.height, -4.0);
    /// ```
    pub fn floor(self: Size) Size {
        return Size.new(math.floor(self.width), math.floor(self.height));
    }

    /// Returns a new `Size`,
    /// with `width` and `height` rounded down towards zero the nearest integer,
    /// unless they are already an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// const size_pos = Size.new(3.3, 3.6).trunc();
    /// try expectEqual(size_pos.width, 3.0);
    /// try expectEqual(size_pos.height, 3.0);
    /// const size_neg = Size.new(-3.3, -3.6).trunc();
    /// try expectEqual(size_neg.width, -3.0);
    /// try expectEqual(size_neg.height, -3.0);
    /// ```
    pub fn trunc(self: Size) Size {
        return Size.new(math.trunc(self.width), math.trunc(self.height));
    }

    /// Returns the aspect ratio of a rectangle with the given size.
    ///
    /// If the width is `0`, the output will be `sign(self.height) * infinity`. If The width and
    /// height are `0`, then the output will be `NaN`.
    pub fn aspectRatio(self: *const Size) f64 {
        return self.height / self.width;
    }

    /// Convert this `Size` into a [`Rect`] with origin `(0.0, 0.0)`.
    pub fn toRect(self: Size) Rect {
        return Rect.new(0.0, 0.0, self.width, self.height);
    }

    /// Convert this `Size` into a [`RoundedRect`] with origin `(0.0, 0.0)` and
    /// the provided corner radius.
    pub fn toRoundRect(self: Size, radii: RoundRectRadii) RoundRect {
        return self.toRect().toRoundRect(radii);
    }

    /// Is this size finite?
    pub fn isFinite(self: *const Size) bool {
        return math.isFinite(self.width) and math.isFinite(self.height);
    }

    /// Is this size NaN?
    pub fn isNan(self: *const Size) bool {
        return math.isNan(self.width) or math.isNan(self.height);
    }
};

test "aspect ratio" {
    const s = Size.new(1.0, 1.0);
    try std.testing.expect(@abs(s.aspectRatio() - 1.0) < 1e-6);
}
