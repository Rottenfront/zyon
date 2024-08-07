const std = @import("std");
const math = std.math;
const Point = @import("point.zig").Point;
const Rect = @import("rect.zig").Rect;
const Vec2 = @import("vec2.zig").Vec2;

/// A 2D affine transform.
pub const Affine = struct {
    a: [6]f64,

    /// The identity transform.
    pub const IDENTITY: Affine = Affine{ .a = [_]f64{ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 } };

    /// A transform that is flipped on the y-axis. Useful for converting between
    /// y-up and y-down spaces.
    pub const FLIP_Y: Affine = Affine{ .a = [_]f64{ 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 } };

    /// A transform that is flipped on the x-axis.
    pub const FLIP_X: Affine = Affine{ .a = [_]f64{ -1.0, 0.0, 0.0, 1.0, 0.0, 0.0 } };

    /// Construct an affine transform from coefficients.
    ///
    /// If the coefficients are `(a, b, c, d, e, f)`, then the resulting
    /// transformation represents this augmented matrix:
    ///
    /// ```text
    /// | a c e |
    /// | b d f |
    /// | 0 0 1 |
    /// ```
    ///
    /// Note that this convention is transposed from PostScript and
    /// Direct2D, but is consistent with the
    /// [Wikipedia](https://en.wikipedia.org/wiki/Affine_transformation)
    /// formulation of affine transformation as augmented matrix.
    pub fn new(a: [6]f64) Affine {
        return Affine{ .a = a };
    }

    /// An affine transform representing uniform scaling.
    pub fn scale(s: f64) Affine {
        return Affine{ .a = [_]f64{ s, 0.0, 0.0, s, 0.0, 0.0 } };
    }

    /// An affine transform representing non-uniform scaling
    /// with different scale values for x and y
    pub fn scale_non_uniform(scale_x: f64, scale_y: f64) Affine {
        return Affine{ .a = [_]f64{ scale_x, 0.0, 0.0, scale_y, 0.0, 0.0 } };
    }

    /// An affine transform representing rotation.
    ///
    /// The convention for rotation is that a positive angle rotates a
    /// positive X direction into positive Y. Thus, in a Y-down coordinate
    /// system (as is common for graphics), it is a clockwise rotation, and
    /// in Y-up (traditional for math), it is anti-clockwise.
    ///
    /// The angle, `th`, is expressed in radians.
    pub fn rotate(th: f64) Affine {
        const s = math.sin(th);
        const c = math.cos(th);
        return Affine{ .a = [_]f64{ c, s, -s, c, 0.0, 0.0 } };
    }

    /// An affine transform representing a rotation of `th` radians about `center`.
    ///
    /// See [`Affine.rotate()`] for more info.
    pub fn rotate_about(th: f64, center: Point) Affine {
        const vector = center.toVec2();
        return Affine.translate(vector.neg()).then_rotate(th).then_translate(vector);
    }

    /// An affine transform representing translation.
    pub fn translate(v: Vec2) Affine {
        return Affine{ .a = [_]f64{ 1.0, 0.0, 0.0, 1.0, v.x, v.y } };
    }

    /// An affine transformation representing a skew.
    ///
    /// The `skew_x` and `skew_y` parameters represent skew factors for the
    /// horizontal and vertical directions, respectively.
    ///
    /// This is commonly used to generate a faux oblique transform for
    /// font rendering. In this case, you can slant the glyph 20 degrees
    /// clockwise in the horizontal direction (assuming a Y-up coordinate
    /// system):
    ///
    /// ```
    /// const obliqueTransform = Affine.skew(math.tan(math.degreesToRadians(20.0)), 0.0);
    /// ```
    pub fn skew(skew_x: f64, skew_y: f64) Affine {
        return Affine{ .a = [_]f64{ 1.0, skew_y, skew_x, 1.0, 0.0, 0.0 } };
    }

    /// Create an affine transform that represents reflection about the line
    /// `point + direction * t, t in (-infty, infty)`
    pub fn reflect(point: Point, direction: Vec2) Affine {
        const n = Vec2.new(direction.y, -direction.x).normalize();

        // Compute Householder reflection matrix
        const x2 = n.x * n.x;
        const xy = n.x * n.y;
        const y2 = n.y * n.y;
        // Here we also add in the post translation, because it doesn't require any further calc.
        const affine = Affine{ .a = [_]f64{
            1.0 - 2.0 * x2,
            -2.0 * xy,
            -2.0 * xy,
            1.0 - 2.0 * y2,
            point.x,
            point.y,
        } };
        return affine.pre_translate(point.toVec2().neg());
    }

    /// A rotation by `th` followed by `self`.
    ///
    /// Equivalent to `self.apply(Affine.rotate(th))`
    pub fn pre_rotate(self: Affine, th: f64) Affine {
        return self.apply(Affine.rotate(th));
    }

    /// A rotation by `th` about `center` followed by `self`.
    ///
    /// Equivalent to `self.apply(Affine.rotate_about(th))`
    pub fn pre_rotate_about(self: Affine, th: f64, center: Point) Affine {
        return self.apply(Affine.rotate_about(th, center));
    }

    /// A scale by `scale` followed by `self`.
    ///
    /// Equivalent to `self.apply(Affine.scale(scale))`
    pub fn pre_scale(self: Affine, s: f64) Affine {
        return self.apply(Affine.scale(s));
    }

    /// A scale by `(scale_x, scale_y)` followed by `self`.
    ///
    /// Equivalent to `self.apply(Affine.scale_non_uniform(scale_x, scale_y))`
    pub fn pre_scale_non_uniform(self: Affine, scale_x: f64, scale_y: f64) Affine {
        return self.apply(Affine.scale_non_uniform(scale_x, scale_y));
    }

    /// A translation of `trans` followed by `self`.
    ///
    /// Equivalent to `Affine.translate(trans).apply(self)`
    pub fn pre_translate(self: Affine, v: Vec2) Affine {
        return self.apply(Affine.translate(v));
    }

    /// `self` followed by a rotation of `th`.
    ///
    /// Equivalent to `Affine::rotate(th).apply(self)`
    pub fn then_rotate(self: Affine, th: f64) Affine {
        return Affine.rotate(th).apply(self);
    }

    /// `self` followed by a rotation of `th` about `center`.
    ///
    /// Equivalent to `Affine::rotate_about(th, center).apply(self)`
    pub fn then_rotate_about(self: Affine, th: f64, center: Point) Affine {
        return Affine.rotate_about(th, center).apply(self);
    }

    /// `self` followed by a scale of `scale`.
    ///
    /// Equivalent to `Affine::scale(scale).apply(self)`
    pub fn then_scale(self: Affine, s: f64) Affine {
        return Affine.scale(s).apply(self);
    }

    /// `self` followed by a scale of `(scale_x, scale_y)`.
    ///
    /// Equivalent to `Affine::scale_non_uniform(scale_x, scale_y).apply(self)`
    pub fn then_scale_non_uniform(self: Affine, scale_x: f64, scale_y: f64) Affine {
        return Affine.scale_non_uniform(scale_x, scale_y).apply(self);
    }

    /// `self` followed by a translation of `trans`.
    ///
    /// Equivalent to `Affine::translate(trans).apply(self)`
    pub fn then_translate(self: Affine, v: Vec2) Affine {
        var affine = self;
        affine.a[4] += v.x;
        affine.a[5] += v.y;
        return affine;
    }

    /// Creates an affine transformation that takes the unit square to the given rectangle.
    ///
    /// Useful when you want to draw into the unit square but have your output fill any rectangle.
    /// In this case push the `Affine` onto the transform stack.
    pub fn mapUnitSquare(rect: Rect) Affine {
        return Affine{ .a = [_]f64{ rect.width(), 0.0, 0.0, rect.height(), rect.x0, rect.y0 } };
    }

    /// Get the coefficients of the transform.
    pub fn asCoeffs(self: Affine) [6]f64 {
        return self.a;
    }

    /// Compute the determinant of this transform.
    pub fn determinant(self: Affine) f64 {
        return self.a[0] * self.a[3] - self.a[1] * self.a[2];
    }

    /// Compute the inverse transform.
    ///
    /// Produces NaN values when the determinant is zero.
    pub fn inverse(self: Affine) Affine {
        const invDet = 1.0 / self.determinant();
        return Affine{ .a = [_]f64{
            invDet * self.a[3],
            -invDet * self.a[1],
            -invDet * self.a[2],
            invDet * self.a[0],
            invDet * (self.a[2] * self.a[5] - self.a[3] * self.a[4]),
            invDet * (self.a[1] * self.a[4] - self.a[0] * self.a[5]),
        } };
    }

    /// Compute the bounding box of a transformed rectangle.
    ///
    /// Returns the minimal `Rect` that encloses the given `Rect` after affine transformation.
    /// If the transform is axis-aligned, then this bounding box is "tight", in other words the
    /// returned `Rect` is the transformed rectangle.
    ///
    /// The returned rectangle always has non-negative width and height.
    pub fn transformRectBbox(self: Affine, rect: Rect) Rect {
        const p00 = Point.new(rect.x0, rect.y0).apply(self);
        const p01 = Point.new(rect.x0, rect.y1).apply(self);
        const p10 = Point.new(rect.x1, rect.y0).apply(self);
        const p11 = Point.new(rect.x1, rect.y1).apply(self);
        return Rect.from_points(p00, p01).rectUnion(Rect.from_points(p10, p11));
    }

    /// Is this map finite?
    pub fn isFinite(self: *const Affine) bool {
        return math.isFinite(self.a[0]) and math.isFinite(self.a[1]) and math.isFinite(self.a[2]) and math.isFinite(self.a[3]) and math.isFinite(self.a[4]) and math.isFinite(self.a[5]);
    }

    /// Is this map NaN?
    pub fn isNan(self: *const Affine) bool {
        return math.isNan(self.a[0]) or math.isNan(self.a[1]) or math.isNan(self.a[2]) or math.isNan(self.a[3]) or math.isNan(self.a[4]) or math.isNan(self.a[5]);
    }

    /// Compute the singular value decomposition of the linear transformation (ignoring the
    /// translation).
    ///
    /// All non-degenerate linear transformations can be represented as
    ///
    ///  1. a rotation about the origin.
    ///  2. a scaling along the x and y axes
    ///  3. another rotation about the origin
    ///
    /// composed together. Decomposing a 2x2 matrix in this way is called a "singular value
    /// decomposition" and is written `U Σ V^T`, where U and V^T are orthogonal (rotations) and Σ
    /// is a diagonal matrix (a scaling).
    ///
    /// Since currently this function is used to calculate ellipse radii and rotation from an
    /// affine map on the unit circle, we don't calculate V^T, since a rotation of the unit (or
    /// any) circle about its center always results in the same circle. This is the reason that an
    /// ellipse mapped using an affine map is always an ellipse.
    ///
    /// Will return NaNs if the matrix (or equivalently the linear map) is singular.
    ///
    /// First part of the return tuple is the scaling, second part is the angle of rotation (in
    /// radians)
    pub fn svd(self: *const Affine) struct { Vec2, f64 } {
        const a = self.a[0];
        const a2 = a * a;
        const b = self.a[1];
        const b2 = b * b;
        const c = self.a[2];
        const c2 = c * c;
        const d = self.a[3];
        const d2 = d * d;
        const ab = a * b;
        const cd = c * d;
        const angle = 0.5 * (2.0 * (ab + cd)).atan2(a2 - b2 + c2 - d2);
        const s1 = a2 + b2 + c2 + d2;
        const s2 = math.sqrt(math.powi(a2 - b2 + c2 - d2, 2) + 4.0 * math.powi(ab + cd, 2));
        return .{
            Vec2{
                .x = math.sqrt(0.5 * (s1 + s2)),
                .y = math.sqrt(0.5 * (s1 - s2)),
            },
            angle,
        };
    }

    /// Returns the translation part of this affine map (`(self.0[4], self.0[5])`).
    pub fn translation(self: Affine) Vec2 {
        return Vec2{ .x = self.a[4], .y = self.a[5] };
    }

    /// Replaces the translation portion of this affine map
    ///
    /// The translation can be seen as being applied after the linear part of the map.
    pub fn with_translation(self: Affine, trans: Vec2) Affine {
        var affine = self;
        affine.a[4] = trans.x;
        affine.a[5] = trans.y;
        return affine;
    }

    /// Apply self on other `Affine`
    pub fn apply(self: Affine, other: Affine) Affine {
        return Affine{ .a = [_]f64{
            self.a[0] * other.a[0] + self.a[2] * other.a[1],
            self.a[1] * other.a[0] + self.a[3] * other.a[1],
            self.a[0] * other.a[2] + self.a[2] * other.a[3],
            self.a[1] * other.a[2] + self.a[3] * other.a[3],
            self.a[0] * other.a[4] + self.a[2] * other.a[5] + self.a[4],
            self.a[1] * other.a[4] + self.a[3] * other.a[5] + self.a[5],
        } };
    }
};

fn assert_near(p0: Point, p1: Point) !void {
    try std.testing.expect(p1.distance(p0) < 1e-9);
}

fn affine_assert_near(a0: Affine, a1: Affine) !void {
    for (0..6) |i| {
        try std.testing.expect(@abs(a0.a[i] - a1.a[i]) < 1e-9);
    }
}

test "affine basic" {
    const p = Point.new(3.0, 4.0);

    assert_near(Affine.IDENTITY * p, p);
    assert_near(Affine.scale(2.0) * p, Point.new(6.0, 8.0));
    assert_near(Affine.rotate(0.0) * p, p);
    assert_near(Affine.rotate(math.pi / 2.0) * p, Point.new(-4.0, 3.0));
    assert_near(Affine.translate(Vec2.new(5.0, 6.0)) * p, Point.new(8.0, 10.0));
    assert_near(Affine.skew(0.0, 0.0) * p, p);
    assert_near(Affine.skew(2.0, 4.0) * p, Point.new(11.0, 16.0));
}

test "affine apply" {
    const a1 = Affine.new([_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
    const a2 = Affine.new([_]f64{ 0.1, 1.2, 2.3, 3.4, 4.5, 5.6 });

    const px = Point.new(1.0, 0.0);
    const py = Point.new(0.0, 1.0);
    const pxy = Point.new(1.0, 1.0);

    assert_near(px.applyAffine(a2).applyAffine(a1), px.applyAffine(a1).applyAffine(a2));
    assert_near(py.applyAffine(a2).applyAffine(a1), py.applyAffine(a1).applyAffine(a2));
    assert_near(pxy.applyAffine(a2).applyAffine(a1), pxy.applyAffine(a1).applyAffine(a2));
}

test "affine invariant" {
    const a = Affine.new([_]f64{ 0.1, 1.2, 2.3, 3.4, 4.5, 5.6 });
    const a_inv = a.inverse();

    const px = Point.new(1.0, 0.0);
    const py = Point.new(0.0, 1.0);
    const pxy = Point.new(1.0, 1.0);
    assert_near(px.applyAffine(a_inv).applyAffine(a), px);
    assert_near(py.applyAffine(a_inv).applyAffine(a), py);
    assert_near(pxy.applyAffine(a_inv).applyAffine(a), pxy);
    assert_near(px.applyAffine(a).applyAffine(a_inv), px);
    assert_near(py.applyAffine(a).applyAffine(a_inv), py);
    assert_near(pxy.applyAffine(a).applyAffine(a_inv), pxy);
}

test "reflection" {
    affine_assert_near(
        Affine.reflect(Point.ZERO, Vec2.new(1.0, 0.0)),
        Affine.new([_]f64{ 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 }),
    );
    affine_assert_near(
        Affine.reflect(Point.ZERO, Vec2.new(0.0, 1.0)),
        Affine.new([_]f64{ -1.0, 0.0, 0.0, 1.0, 0.0, 0.0 }),
    );
    // y = x
    affine_assert_near(
        Affine.reflect(Point.ZERO, Vec2.new(1.0, 1.0)),
        Affine.new([_]f64{ 0.0, 1.0, 1.0, 0.0, 0.0, 0.0 }),
    );

    // no translate
    const point = Point.new(0.0, 0.0);
    const vec = Vec2.new(1.0, 1.0);
    const map = Affine.reflect(point, vec);
    assert_near(Point.new(0.0, 0.0).applyAffine(map), Point.new(0.0, 0.0));
    assert_near(Point.new(1.0, 1.0).applyAffine(map), Point.new(1.0, 1.0));
    assert_near(Point.new(1.0, 2.0).applyAffine(map), Point.new(2.0, 1.0));

    // with translate
    const point1 = Point.new(1.0, 0.0);
    const vec1 = Vec2.new(1.0, 1.0);
    const map1 = Affine.reflect(point1, vec1);
    assert_near(Point.new(1.0, 0.0).applyAffine(map1), Point.new(1.0, 0.0));
    assert_near(Point.new(2.0, 1.0).applyAffine(map1), Point.new(2.0, 1.0));
    assert_near(Point.new(2.0, 2.0).applyAffine(map1), Point.new(3.0, 1.0));
}
