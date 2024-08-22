const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;
const Point = mod.Point;
const Vec2 = mod.Vec2;

/// A transformation consisting of a uniform scaling followed by a translation.
///
/// If the translation is `(x, y)` and the scale is `s`, then this
/// transformation represents this augmented matrix:
///
/// ```text
/// | s 0 x |
/// | 0 s y |
/// | 0 0 1 |
/// ```
///
/// See [`Affine`] for more details about the
/// equivalence with augmented matrices.
///
/// Various multiplication ops are defined, and these are all defined
/// to be consistent with matrix multiplication. Therefore,
/// `TranslateScale * Point` is defined but not the other way around.
///
/// Also note that multiplication is not commutative.
///
/// This transformation is less powerful than `Affine`, but can be applied
/// to more primitives, especially including [`Rect`].
pub const TranslateScale = struct {
    /// The translation component of this transformation
    translation: Vec2,
    /// The scale component of this transformation
    scale: f64,

    pub fn new(translation: Vec2, s: f64) TranslateScale {
        return TranslateScale{ .translation = translation, .scale = s };
    }

    pub fn scale(s: f64) TranslateScale {
        return TranslateScale{ .translation = Vec2.ZERO, .scale = s };
    }

    pub fn translate(translation: Vec2) TranslateScale {
        return TranslateScale{ .translation = translation, .scale = 1.0 };
    }

    /// Create a transform that scales about a point other than the origin.
    ///
    /// # Examples
    ///
    /// ```
    /// const center = Point.new(1.0, 1.0);
    /// const ts = TranslateScale.fromScaleAbout(2.0, center);
    /// // Should keep the point (1.0, 1.0) stationary
    /// assert_near(center.applyTranslateScale(ts), center);
    /// // (2.0, 2.0) -> (3.0, 3.0)
    /// assert_near(Point.new(2.0, 2.0).applyTranslateScale(ts), Point.new(3.0, 3.0));
    /// ```
    pub fn fromScaleSbout(s: f64, focus: Point) TranslateScale {
        // We need to create a transform that is equivalent to translating `focus`
        // to the origin, followed by a normal scale, followed by reversing the translation.
        // We need to find the (translation ∘ scale) that matches this.
        const focus_vector = focus.toVec2();
        const translation = focus_vector.sub(focus_vector.mul(s));
        return TranslateScale{ .translation = translation, .scale = s };
    }

    /// Compute the inverse transform.
    ///
    /// Multiplying a transform with its inverse (either on the
    /// left or right) results in the identity transform
    /// (modulo floating point rounding errors).
    ///
    /// Produces NaN values when scale is zero.
    pub fn inverse(self: TranslateScale) TranslateScale {
        const scale_recip = 1.0 / self.scale;
        return TranslateScale{
            .translation = self.translation.mul(-scale_recip),
            .scale = scale_recip,
        };
    }

    /// Is this translate/scale finite?
    pub fn isFinite(self: *TranslateScale) bool {
        return self.translation.isFinite() and math.isFinite(self.scale);
    }

    /// Is this translate/scale NaN?
    pub fn isNan(self: *TranslateScale) bool {
        return self.translation.isNan() or math.isNan(self.scale);
    }

    /// Apply self on other `TranslateScale`
    pub fn apply(self: TranslateScale, other: TranslateScale) TranslateScale {
        return TranslateScale{ .translation = self.translation.sum(other.translation.mul(self.scale)), .scale = self.scale * other.scale };
    }
};

//    use crate::{Affine, Point, TranslateScale, Vec2};
//
//    fn assert_near(p0: Point, p1: Point) {
//        assert!((p1 - p0).hypot() < 1e-9, "{p0:?} != {p1:?}");
//    }
//
//    #[test]
//    fn translate_scale() {
//        let p = Point::new(3.0, 4.0);
//        let ts = TranslateScale::new(Vec2::new(5.0, 6.0), 2.0);
//
//        assert_near(ts * p, Point::new(11.0, 14.0));
//    }
//
//    #[test]
//    fn conversions() {
//        let p = Point::new(3.0, 4.0);
//        let s = 2.0;
//        let t = Vec2::new(5.0, 6.0);
//        let ts = TranslateScale::new(t, s);
//
//        // Test that conversion to affine is consistent.
//        let a: Affine = ts.into();
//        assert_near(ts * p, a * p);
//
//        assert_near((s * p.to_vec2()).to_point(), TranslateScale::scale(s) * p);
//        assert_near(p + t, TranslateScale::translate(t) * p);
//    }
//
//    #[test]
//    fn inverse() {
//        let p = Point::new(3.0, 4.0);
//        let ts = TranslateScale::new(Vec2::new(5.0, 6.0), 2.0);
//
//        assert_near(p, (ts * ts.inverse()) * p);
//        assert_near(p, (ts.inverse() * ts) * p);
//    }

fn assert_near(p0: Point, p1: Point) !void {
    try std.testing.expect(p1.distance(p0) < 1e-9);
}

test "translate/scale" {
    const p = Point.new(3.0, 4.0);
    const ts = TranslateScale.new(Vec2.new(5.0, 6.0), 2.0);

    assert_near(p.applyTranslateScale(ts), Point.new(11.0, 14.0));
}
