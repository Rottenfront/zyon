const std = @import("std");
const math = std.math;

pub const ArrayAccessError = error{OutOfBounds};

/// A 2D affine transform.
pub const Affine = struct {
    a: [6]f64,

    /// The identity transform.
    pub const IDENTITY: @This() = @This().scale(1.0);

    /// A transform that is flipped on the y-axis. Useful for converting between
    /// y-up and y-down spaces.
    pub const FLIP_Y: @This() = @This().new([_]f64{ 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 });

    /// A transform that is flipped on the x-axis.
    pub const FLIP_X: @This() = @This().new([_]f64{ -1.0, 0.0, 0.0, 1.0, 0.0, 0.0 });

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
    pub fn new(a: [6]f64) @This() {
        return @This(){ .a = a };
    }

    /// An affine transform representing uniform scaling.
    pub fn scale(s: f64) @This() {
        return @This(){ .a = [_]f64{ s, 0.0, 0.0, s, 0.0, 0.0 } };
    }

    /// An affine transform representing non-uniform scaling
    /// with different scale values for x and y
    pub fn scale_non_uniform(scale_x: f64, scale_y: f64) @This() {
        return @This(){ .a = [_]f64{ scale_x, 0.0, 0.0, scale_y, 0.0, 0.0 } };
    }

    /// An affine transform representing rotation.
    ///
    /// The convention for rotation is that a positive angle rotates a
    /// positive X direction into positive Y. Thus, in a Y-down coordinate
    /// system (as is common for graphics), it is a clockwise rotation, and
    /// in Y-up (traditional for math), it is anti-clockwise.
    ///
    /// The angle, `th`, is expressed in radians.
    pub fn rotate(th: f64) @This() {
        const s = math.sin(th);
        const c = math.cos(th);
        return @This(){ .a = [_]f64{ c, s, -s, c, 0.0, 0.0 } };
    }

    /// An affine transform representing a rotation of `th` radians about `center`.
    ///
    /// See [`Affine.rotate()`] for more info.
    pub fn rotate_about(th: f64, center: Point) @This() {
        const vector = center.toVec2();
        return @This().translate(vector.neg()).then_rotate(th).then_translate(vector);
    }

    /// An affine transform representing translation.
    pub fn translate(v: Vec2) @This() {
        return @This(){ .a = [_]f64{ 1.0, 0.0, 0.0, 1.0, v.x, v.y } };
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
    pub fn skew(skew_x: f64, skew_y: f64) @This() {
        return @This(){ .a = [_]f64{ 1.0, skew_y, skew_x, 1.0, 0.0, 0.0 } };
    }

    /// Create an affine transform that represents reflection about the line
    /// `point + direction * t, t in (-infty, infty)`
    pub fn reflect(point: Point, direction: Vec2) @This() {
        const n = Vec2.new(direction.y, -direction.x).normalize();

        // Compute Householder reflection matrix
        const x2 = n.x * n.x;
        const xy = n.x * n.y;
        const y2 = n.y * n.y;
        // Here we also add in the post translation,
        // because it doesn't require any further calc.
        const affine = @This(){ .a = [_]f64{
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
    pub fn pre_rotate(self: @This(), th: f64) @This() {
        return self.apply(@This().rotate(th));
    }

    /// A rotation by `th` about `center` followed by `self`.
    pub fn pre_rotate_about(self: @This(), th: f64, center: Point) @This() {
        return self.apply(@This().rotate_about(th, center));
    }

    /// A scale by `scale` followed by `self`.
    pub fn pre_scale(self: @This(), s: f64) @This() {
        return self.apply(@This().scale(s));
    }

    /// A scale by `(scale_x, scale_y)` followed by `self`.
    pub fn pre_scale_non_uniform(self: @This(), scale_x: f64, scale_y: f64) @This() {
        return self.apply(@This().scale_non_uniform(scale_x, scale_y));
    }

    /// A translation of `trans` followed by `self`.
    pub fn pre_translate(self: @This(), v: Vec2) @This() {
        return self.apply(@This().translate(v));
    }

    /// `self` followed by a rotation of `th`.
    pub fn then_rotate(self: @This(), th: f64) @This() {
        return @This().rotate(th).apply(self);
    }

    /// `self` followed by a rotation of `th` about `center`.
    pub fn then_rotate_about(self: @This(), th: f64, center: Point) @This() {
        return @This().rotate_about(th, center).apply(self);
    }

    /// `self` followed by a scale of `scale`.
    pub fn then_scale(self: @This(), s: f64) @This() {
        return @This().scale(s).apply(self);
    }

    /// `self` followed by a scale of `(scale_x, scale_y)`.
    pub fn then_scale_non_uniform(self: @This(), scale_x: f64, scale_y: f64) @This() {
        return @This().scale_non_uniform(scale_x, scale_y).apply(self);
    }

    /// `self` followed by a translation of `trans`.
    pub fn then_translate(self: @This(), v: Vec2) @This() {
        var affine = self;
        affine.a[4] += v.x;
        affine.a[5] += v.y;
        return affine;
    }

    /// Creates an affine transformation that takes the unit square to the given rectangle.
    ///
    /// Useful when you want to draw into the unit square but have your output fill any rectangle.
    /// In this case push the `Affine` onto the transform stack.
    pub fn mapUnitSquare(rect: Rect) @This() {
        return @This(){ .a = [_]f64{ rect.width(), 0.0, 0.0, rect.height(), rect.x0, rect.y0 } };
    }

    /// Get the coefficients of the transform.
    pub fn asCoeffs(self: @This()) [6]f64 {
        return self.a;
    }

    /// Compute the determinant of this transform.
    pub fn determinant(self: @This()) f64 {
        return self.a[0] * self.a[3] - self.a[1] * self.a[2];
    }

    /// Compute the inverse transform.
    ///
    /// Produces NaN values when the determinant is zero.
    pub fn inverse(self: @This()) @This() {
        const invDet = 1.0 / self.determinant();
        return @This(){ .a = [_]f64{
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
    pub fn transformRectBbox(self: @This(), rect: Rect) Rect {
        const p00 = Point.new(rect.x0, rect.y0).apply(self);
        const p01 = Point.new(rect.x0, rect.y1).apply(self);
        const p10 = Point.new(rect.x1, rect.y0).apply(self);
        const p11 = Point.new(rect.x1, rect.y1).apply(self);
        return Rect.from_points(p00, p01).rectUnion(Rect.from_points(p10, p11));
    }

    /// Is this map finite?
    pub fn isFinite(self: *const @This()) bool {
        return math.isFinite(self.a[0]) and math.isFinite(self.a[1]) and math.isFinite(self.a[2]) and math.isFinite(self.a[3]) and math.isFinite(self.a[4]) and math.isFinite(self.a[5]);
    }

    /// Is this map NaN?
    pub fn isNan(self: *const @This()) bool {
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
    pub fn svd(self: *const @This()) struct { Vec2, f64 } {
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
    pub fn translation(self: @This()) Vec2 {
        return Vec2{ .x = self.a[4], .y = self.a[5] };
    }

    /// Replaces the translation portion of this affine map
    ///
    /// The translation can be seen as being applied after the linear part of the map.
    pub fn with_translation(self: @This(), trans: Vec2) @This() {
        var affine = self;
        affine.a[4] = trans.x;
        affine.a[5] = trans.y;
        return affine;
    }

    /// Apply self on other `@This()`
    pub fn apply(self: @This(), other: @This()) @This() {
        return @This(){ .a = [_]f64{
            self.a[0] * other.a[0] + self.a[2] * other.a[1],
            self.a[1] * other.a[0] + self.a[3] * other.a[1],
            self.a[0] * other.a[2] + self.a[2] * other.a[3],
            self.a[1] * other.a[2] + self.a[3] * other.a[3],
            self.a[0] * other.a[4] + self.a[2] * other.a[5] + self.a[4],
            self.a[1] * other.a[4] + self.a[3] * other.a[5] + self.a[5],
        } };
    }
};

/// A single arc segment.
pub const Arc = struct {
    /// The arc's centre point.
    center: Point,
    /// The arc's radii, where the vector's x-component is the radius in the
    /// positive x direction after applying `x_rotation`.
    radii: Vec2,
    /// The start angle in radians.
    start_angle: f64,
    /// The angle between the start and end of the arc, in radians.
    sweep_angle: f64,
    /// How much the arc is rotated, in radians.
    x_rotation: f64,

    /// Create a new `Arc`
    pub fn new(center: Point, radii: Vec2, start_angle: f64, sweep_angle: f64, x_rotation: f64) @This() {
        return @This(){
            .center = center,
            .radii = radii,
            .start_angle = start_angle,
            .sweep_angle = sweep_angle,
            .x_rotation = x_rotation,
        };
    }

    /// Create an iterator generating Bezier path elements.
    ///
    /// The generated elements can be appended to an existing bezier path.
    pub fn appendIter(self: *const @This(), tolerance: f64) AppendIter {
        const sign = math.sign(self.sweep_angle);
        const scaled_err = @max(self.radii.x, self.radii.y) / tolerance;
        // Number of subdivisions per ellipse based on error tolerance.
        // Note: this may slightly underestimate the error for quadrants.
        const n_err = @max(math.pow(1.1163 * scaled_err, 1.0 / 6.0), 3.999_999);
        const n_f64 = math.ceil(n_err * @abs(self.sweep_angle) * (1.0 / (2.0 * math.pi)));
        const angle_step = self.sweep_angle / n_f64;
        const n: usize = math.round(n_f64);
        const arm_len = (4.0 / 3.0) * math.tan(@abs(0.25 * angle_step)) * sign;
        const angle0 = self.start_angle;
        const p0 = util.sample_ellipse(self.radii, self.x_rotation, angle0);

        return AppendIter{
            .idx = 0,

            .center = self.center,
            .radii = self.radii,
            .x_rotation = self.x_rotation,
            .n = n,
            .arm_len = arm_len,
            .angle_step = angle_step,

            .p0 = p0,
            .angle0 = angle0,
        };
    }

    /// Converts an Arc into a list of cubic bezier segments.
    pub fn to_cubic_beziers(
        self: @This(),
        tolerance: f64,
        allocator: ?std.mem.Allocator,
    ) !std.ArrayList(BezPath.Element) {
        var path = self.appendIter(tolerance);
        var result = std.ArrayList(BezPath.Element).init(alloc: {
            if (allocator != null) {
                break :alloc allocator;
            }
            break :alloc std.heap.page_allocator;
        });

        var current_elem = path.next();
        while (current_elem != null) {
            result.append(current_elem);
            current_elem = path.next();
        }

        return result;
    }

    pub const AppendIter = struct {
        idx: isize,

        center: Point,
        radii: Vec2,
        x_rotation: f64,
        n: usize,
        arm_len: f64,
        angle_step: f64,

        p0: Vec2,
        angle0: f64,

        pub fn next(self: *@This()) ?BezPath.Element {
            if (self.idx == -1) {
                self.idx += 1;
                return BezPath.Element{
                    .move_to = struct {
                        .p = self.p0,
                    },
                };
            }

            if (self.idx >= self.n) {
                return null;
            }

            const angle1 = self.angle0 + self.angle_step;
            const p0 = self.p0;
            const p1 = p0 + util.sample_ellipse(
                self.radii,
                self.x_rotation,
                self.angle0 + math.pi / 2.0,
            ).mul(self.arm_len);
            const p3 = util.sample_ellipse(self.radii, self.x_rotation, angle1);
            const p2 = p3 - util.sample_ellipse(
                self.radii,
                self.x_rotation,
                self.angle1 + math.pi / 2.0,
            ).mul(self.arm_len);

            self.angle0 = angle1;
            self.p0 = p3;
            self.idx += 1;

            return BezPath.Element{ .curve_to = struct {
                .p0 = self.center.sum(p1),
                .p1 = self.center.sum(p2),
                .end = self.center.sum(p3),
            } };
        }
    };

    pub fn pathElements(self: *const @This(), tolerance: f64) AppendIter {
        var iter = self.appendIter(tolerance);
        iter.idx = -1;
        return iter;
    }

    /// Note: shape isn't closed so area is not well defined.
    pub fn area(self: *const @This()) f64 {
        return math.pi * self.radii.x * self.radii.y;
    }

    /// The perimeter of the arc.
    ///
    /// For now we just approximate by using the bezier curve representation.
    pub fn perimeter(self: *const @This(), accuracy: f64) f64 {
        return self.pathElements(0.1).perimeter(accuracy);
    }
};

/// A Bézier path.
///
/// These docs assume basic familiarity with Bézier curves; for an introduction,
/// see Pomax's wonderful [A Primer on Bézier Curves].
///
/// This path can contain lines, quadratics ([`QuadBez`]) and cubics
/// ([`CubicBez`]), and may contain multiple subpaths.
///
/// # Elements and Segments
///
/// A Bézier path can be represented in terms of either 'elements' ([`PathEl`])
/// or 'segments' ([`PathSeg`]). Elements map closely to how Béziers are
/// generally used in PostScript-style drawing APIs; they can be thought of as
/// instructions for drawing the path. Segments more directly describe the
/// path itself, with each segment being an independent line or curve.
///
/// These different representations are useful in different contexts.
/// For tasks like drawing, elements are a natural fit, but when doing
/// hit-testing or subdividing, we need to have access to the segments.
///
/// Conceptually, a `BezPath` contains zero or more subpaths. Each subpath
/// *always* begins with a `MoveTo`, then has zero or more `LineTo`, `QuadTo`,
/// and `CurveTo` elements, and optionally ends with a `ClosePath`.
///
/// Internally, a `BezPath` is a list of [`PathEl`]s.
///
/// # Advanced functionality
///
/// In addition to the basic API, there are several useful pieces of advanced
/// functionality available on `BezPath`:
///
/// - [`flatten`] does Bézier flattening, converting a curve to a series of
///   line segments
/// - [`intersect_line`] computes intersections of a path with a line, useful
///   for things like subdividing
///
/// [A Primer on Bézier Curves]: https://pomax.github.io/bezierinfo/
pub const BezPath = struct {
    pub const PathList: type = std.ArrayList(Element);
    pub const allocator = std.heap.page_allocator;

    data: PathList,

    pub fn new() @This() {
        return @This(){
            .data = PathList.init(@This().allocator),
        };
    }

    pub fn pop(self: *@This()) ?Element {
        return self.data.popOrNull();
    }

    pub fn push(self: *@This(), el: Element) @This().allocator.Error!void {
        try self.data.append(el);
        std.debug.assert(match: {
            switch (self.data.items[0]) {
                .move_to => |_| break :match true,
                else => break :match false,
            }
        });
    }

    pub fn move_to(self: *@This(), p: Point) @This().allocator.Error!void {
        return self.push(Element{ .move_to = .{ .p = p } });
    }

    pub fn line_to(self: *@This(), end: Point) @This().allocator.Error!void {
        return self.push(Element{ .line_to = .{ .end = end } });
    }

    pub fn quad_to(
        self: *@This(),
        p0: Point,
        end: Point,
    ) @This().allocator.Error!void {
        return self.push(Element{ .quad_to = .{ .p0 = p0, .end = end } });
    }

    pub fn curve_to(
        self: *@This(),
        p0: Point,
        p1: Point,
        end: Point,
    ) @This().allocator.Error!void {
        return self.push(Element{ .curve_to = .{
            .p0 = p0,
            .p1 = p1,
            .end = end,
        } });
    }

    pub fn close_path(self: *@This()) @This().allocator.Error!void {
        return self.push(Element{.close_path});
    }

    /// The element of a Bézier path.
    ///
    /// A valid path has `MoveTo` at the beginning of each subpath.
    pub const Element = union(enum) {
        /// Move directly to the point without drawing anything, starting a new
        /// subpath.
        move_to: struct { p: Point },
        /// Draw a line from the current location to the point.
        line_to: struct { end: Point },
        /// Draw a quadratic bezier using the current location and the two points.
        quad_to: struct { p0: Point, end: Point },
        /// Draw a cubic bezier using the current location and the three points.
        curve_to: struct { p0: Point, p1: Point, end: Point },
        /// Close off the path.
        close_path,
    };

    /// A segment of a Bézier path.
    pub const Segment = union(enum) {
        /// A line segment.
        line: Line,
        /// A quadratic bezier segment.
        quad: QuadBez,
        /// A cubic bezier segment.
        curve: CubicBez,
    };

    /// An intersection of a [`Line`] and a [`PathSeg`].
    ///
    /// This can be generated with the [`PathSeg::intersect_line`] method.
    pub const LineIntersection = struct {
        /// The 'time' that the intersection occurs, on the line.
        ///
        /// This value is in the range 0..1.
        line_t: f64,

        /// The 'time' that the intersection occurs, on the path segment.
        ///
        /// This value is nominally in the range 0..1, although it may slightly exceed
        /// that range at the boundaries of segments.
        segment_t: f64,
    };

    /// The minimum distance between two Bézier curves.
    pub const MinDistance = struct {
        /// The shortest distance between any two points on the two curves.
        distance: f64,
        /// The position of the nearest point on the first curve, as a parameter.
        ///
        /// To resolve this to a [`Point`], use [`ParamCurve.eval`].
        t1: f64,
        /// The position of the nearest point on the second curve, as a parameter.
        ///
        /// To resolve this to a [`Point`], use [`ParamCurve.eval`].
        t2: f64,
    };
};

pub const CubicBez = struct {
    start: Point,
    p0: Point,
    p1: Point,
    end: Point,

    pub fn new(start: Point, p0: Point, p1: Point, end: Point) @This() {
        return @This(){ .start = start, .p0 = p0, .p1 = p1, .end = end };
    }

    /// Classification result for cusp detection.
    pub const CuspType = enum {
        /// Cusp is a loop.
        Loop,
        /// Cusp has two closely spaced inflection points.
        DoubleInflection,
    };
};

pub const Line = struct {
    start: Point,
    end: Point,

    pub fn new(start: Point, end: Point) @This() {
        return @This(){ .start = start, .end = end };
    }
};

/// A single quadratic Bézier segment.
pub const QuadBez = struct {
    start: Point,
    p0: Point,
    end: Point,

    pub fn new(start: Point, p0: Point, end: Point) @This() {
        return @This(){ .start = start, .p0 = p0, .end = end };
    }

    /// Raise the order by 1.
    ///
    /// Returns a cubic Bézier segment that exactly represents this quadratic.
    pub fn raise(self: @This()) CubicBez {
        return CubicBez.new(
            self.start,
            self.start.sum(self.p0.sub(self.start).mul(2.0 / 3.0)),
            self.end.sum(self.p0.sub(self.end).mul(2.0 / 3.0)),
            self.end,
        );
    }

    /// Estimate the number of subdivisions for flattening.
    pub fn estimateSubdiv(self: *const @This(), sqrt_tol: f64) FlattenParams {
        // Determine transformation to $y = x^2$ parabola.
        const dstart_p0 = self.p0.sub(self.start).toVec2();
        const dp0_end = self.end.sub(self.p0).toVec2();
        const dd = dstart_p0.sub(dp0_end);
        const cross = self.end.sub(self.start).toVec2().cross(dd);
        const xstart = dstart_p0.dot(dd) / cross;
        const xend = dp0_end.dot(dd) / cross;
        const scale = @abs(cross / (dd.hypot() * (xend - xstart)));

        // Compute number of subdivisions needed.
        const astart = util.approx_parabola_integral(xstart);
        const aend = util.approx_parabola_integral(xend);
        const val = val: {
            if (math.isFinite(scale)) {
                const da = @abs(aend - astart);
                const sqrt_scale = @sqrt(scale);
                if (xstart * xend > 0.0) {
                    break :val da * sqrt_scale;
                } else {
                    const xmin = sqrt_tol / sqrt_scale;
                    break :val sqrt_tol * da / util.approx_parabola_integral(xmin);
                }
            }
            break :val 0.0;
        };
        const ustart = util.approx_parabola_inv_integral(astart);
        const uend = util.approx_parabola_inv_integral(aend);

        const uscale = 1 / (uend - ustart);

        return FlattenParams{
            .astart = astart,
            .aend = aend,
            .ustart = ustart,
            .uscale = uscale,
            .val = val,
        };
    }

    // Maps a value from 0..1 to 0..1.
    pub fn determineSubdivT(
        params: *const FlattenParams,
        x: f64,
    ) f64 {
        const a = params.astart + (params.aend - params.astart) * x;
        const u = util.approx_parabola_inv_integral(a);
        return (u - params.ustart) * params.uscale;
    }

    /// Is this quadratic Bezier curve finite?
    pub fn isFinite(self: *const @This()) bool {
        return self.start.isFinite() and self.p0.isFinite() and self.end.isFinite();
    }

    /// Is this quadratic Bezier curve NaN?
    pub fn isNan(self: *const @This()) bool {
        return self.start.isNan() or self.p0.isNan() or self.end.isNan();
    }

    pub fn pathElements(self: *const @This()) QuadBezIter {
        return QuadBezIter{ .quad = *self, .ix = 0 };
    }

    pub fn area(self: *const @This()) f64 {
        return 0.0;
    }

    pub fn perimeter(self: *const @This(), accuracy: f64) f64 {
        return self.arclen(accuracy);
    }

    pub fn winding(self: *const @This(), pt: Point) i32 {
        return 0;
    }

    pub fn boundBox(self: *const @This()) Rect {}

    /// An iterator for quadratic beziers.
    pub const QuadBezIter = struct {
        quad: QuadBez,
        ix: usize,
    };

    pub const FlattenParams = struct {
        astart: f64,
        aend: f64,
        ustart: f64,
        uscale: f64,
        val: f64,
    };
};

pub const Circle = struct {
    center: Point,
    radius: f64,
};

pub const Ellipse = struct {
    inner: Affine,
};

/// The offset curve of a cubic Bézier.
///
/// This is a representation of the offset curve of a cubic Bézier segment, for
/// purposes of curve fitting.
///
/// See the module-level documentation for a bit more discussion of the approach,
/// and how this struct is to be used.
pub const CubicOffset = struct {
    /// Source curve.
    c: CubicBez,
    /// Derivative of source curve.
    q: QuadBez,
    /// Offset.
    d: f64,
    // c0 + c1 t + c2 t^2 is the cross product of second and first
    // derivatives of the underlying cubic, multiplied by offset (for
    // computing cusp).
    c0: f64,
    c1: f64,
    c2: f64,
};

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
    pub fn trunc(self: @This()) @This() {
        return @This().new(math.trunc(self.y), math.trunc(self.y));
    }

    /// Is this point finite?
    pub fn isFinite(self: *const @This()) bool {
        return math.isFinite(self.x) and math.isFinite(self.y);
    }

    /// Is this point NaN?
    pub fn isNan(self: *const @This()) bool {
        return math.isNan(self.x) or math.isNan(self.y);
    }

    pub fn sum(self: @This(), other: @This()) @This() {
        return @This(){ .x = self.x + other.x, .y = self.y + other.y };
    }

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

pub const QuadSpline = struct {
    data: std.ArrayList(Point),

    pub const ToQuadBez = struct {
        idx: usize,
        points: *std.ArrayList(Point),
    };
};

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

    pub fn toRoundRect(self: Rect, radii: RoundRect.Radii) RoundRect {
        return RoundRect.new(self, radii);
    }
};

pub const RoundRect = struct {
    rect: Rect,
    radii: Radii,

    pub fn new(rect: Rect, radii: Radii) RoundRect {
        return RoundRect{ .rect = rect, .radii = radii };
    }

    pub const Radii = struct {
        /// The radius of the top-left corner.
        top_left: f64,
        /// The radius of the top-right corner.
        top_right: f64,
        /// The radius of the bottom-right corner.
        bottom_right: f64,
        /// The radius of the bottom-left corner.
        bottom_left: f64,

        pub const ZERO: Radii = Radii.fromSingleRadius(0.0);

        /// Create a new `Radii`. This function takes radius values for
        /// the four corners. The argument order is `top_left`, `top_right`,
        /// `bottom_right`, `bottom_left`, or clockwise starting from `top_left`.
        pub fn new(top_left: f64, top_right: f64, bottom_right: f64, bottom_left: f64) Radii {
            return Radii{
                .top_left = top_left,
                .top_right = top_right,
                .bottom_right = bottom_right,
                .bottom_left = bottom_left,
            };
        }

        /// Create a new `Radii` from a single radius. The `radius`
        /// argument will be set as the radius for all four corners.
        pub fn fromSingleRadius(radius: f64) Radii {
            return Radii{
                .top_left = radius,
                .top_right = radius,
                .bottom_right = radius,
                .bottom_left = radius,
            };
        }
    };
};

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
    pub fn toRoundRect(self: Size, radii: RoundRect.Radii) RoundRect {
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

    pub fn new(translation: Vec2, s: f64) @This() {
        return @This(){ .translation = translation, .scale = s };
    }

    pub fn scale(s: f64) @This() {
        return @This(){ .translation = Vec2.ZERO, .scale = s };
    }

    pub fn translate(translation: Vec2) @This() {
        return @This(){ .translation = translation, .scale = 1.0 };
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
    pub fn fromScaleSbout(s: f64, focus: Point) @This() {
        // We need to create a transform that is equivalent to translating `focus`
        // to the origin, followed by a normal scale, followed by reversing the translation.
        // We need to find the (translation ∘ scale) that matches this.
        const focus_vector = focus.toVec2();
        const translation = focus_vector.sub(focus_vector.mul(s));
        return @This(){ .translation = translation, .scale = s };
    }

    /// Compute the inverse transform.
    ///
    /// Multiplying a transform with its inverse (either on the
    /// left or right) results in the identity transform
    /// (modulo floating point rounding errors).
    ///
    /// Produces NaN values when scale is zero.
    pub fn inverse(self: @This()) @This() {
        const scale_recip = 1.0 / self.scale;
        return @This(){
            .translation = self.translation.mul(-scale_recip),
            .scale = scale_recip,
        };
    }

    /// Is this translate/scale finite?
    pub fn isFinite(self: *@This()) bool {
        return self.translation.isFinite() and math.isFinite(self.scale);
    }

    /// Is this translate/scale NaN?
    pub fn isNan(self: *@This()) bool {
        return self.translation.isNan() or math.isNan(self.scale);
    }

    /// Apply self on other `TranslateScale`
    pub fn apply(self: @This(), other: @This()) @This() {
        return @This(){ .translation = self.translation.sum(other.translation.mul(self.scale)), .scale = self.scale * other.scale };
    }

    pub fn toAffine(self: @This()) Affine {
        return Affine{ .a = [_]f64{
            self.scale,
            0.0,
            0.0,
            self.scale,
            self.translation.x,
            self.translation.y,
        } };
    }
};

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

pub const util = struct {
    /// An approximation to $\int (1 + 4x^2) ^ -0.25 dx$
    ///
    /// This is used for flattening curves.
    pub fn approx_parabola_integral(x: f64) f64 {
        const D: f64 = 0.67;
        return x / (1.0 - D + @sqrt(@sqrt(math.powi(D, 4) + 0.25 * x * x)));
    }

    /// An approximation to the inverse parabola integral.
    pub fn approx_parabola_inv_integral(x: f64) f64 {
        const B: f64 = 0.39;
        return x * (1.0 - B + @sqrt(B * B + 0.25 * x * x));
    }

    /// Take the ellipse radii, how the radii are rotated, and the sweep angle, and return a point on
    /// the ellipse.
    pub fn sample_ellipse(radii: Vec2, x_rotation: f64, angle: f64) Vec2 {
        const angle_sin = math.sin(angle);
        const angle_cos = math.cos(angle);
        const u = radii.x * angle_cos;
        const v = radii.y * angle_sin;
        return rotate_pt(Vec2.new(u, v), x_rotation);
    }

    /// Rotate `pt` about the origin by `angle` radians.
    pub fn rotate_pt(pt: Vec2, angle: f64) Vec2 {
        const angle_sin = math.sin(angle);
        const angle_cos = math.cos(angle);
        return Vec2.new(
            pt.x * angle_cos - pt.y * angle_sin,
            pt.x * angle_sin + pt.y * angle_cos,
        );
    }

    fn expect_near(p0: Point, p1: Point) !void {
        try std.testing.expect(p1.distance(p0) < 1e-9);
    }

    fn affine_expect_near(a0: Affine, a1: Affine) !void {
        for (0..6) |i| {
            try std.testing.expect(@abs(a0.a[i] - a1.a[i]) < 1e-9);
        }
    }
};

test "affine basic" {
    const p = Point.new(3.0, 4.0);

    try util.expect_near(p.applyAffine(Affine.IDENTITY), p);
    try util.expect_near(p.applyAffine(Affine.scale(2.0)), Point.new(6.0, 8.0));
    try util.expect_near(p.applyAffine(Affine.rotate(0.0)), p);
    try util.expect_near(p.applyAffine(Affine.rotate(math.pi / 2.0)), Point.new(-4.0, 3.0));
    try util.expect_near(p.applyAffine(Affine.translate(Vec2.new(5.0, 6.0))), Point.new(8.0, 10.0));
    try util.expect_near(p.applyAffine(Affine.skew(0.0, 0.0)), p);
    try util.expect_near(p.applyAffine(Affine.skew(2.0, 4.0)), Point.new(11.0, 16.0));
}

test "affine apply" {
    const a1 = Affine.new([_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
    const a2 = Affine.new([_]f64{ 0.1, 1.2, 2.3, 3.4, 4.5, 5.6 });

    const px = Point.new(1.0, 0.0);
    const py = Point.new(0.0, 1.0);
    const pxy = Point.new(1.0, 1.0);

    try util.expect_near(px.applyAffine(a2).applyAffine(a1), px.applyAffine(a1.apply(a2)));
    try util.expect_near(py.applyAffine(a2).applyAffine(a1), py.applyAffine(a1.apply(a2)));
    try util.expect_near(pxy.applyAffine(a2).applyAffine(a1), pxy.applyAffine(a1.apply(a2)));
}

test "affine invariant" {
    const a = Affine.new([_]f64{ 0.1, 1.2, 2.3, 3.4, 4.5, 5.6 });
    const a_inv = a.inverse();

    const px = Point.new(1.0, 0.0);
    const py = Point.new(0.0, 1.0);
    const pxy = Point.new(1.0, 1.0);
    try util.expect_near(px.applyAffine(a_inv).applyAffine(a), px);
    try util.expect_near(py.applyAffine(a_inv).applyAffine(a), py);
    try util.expect_near(pxy.applyAffine(a_inv).applyAffine(a), pxy);
    try util.expect_near(px.applyAffine(a).applyAffine(a_inv), px);
    try util.expect_near(py.applyAffine(a).applyAffine(a_inv), py);
    try util.expect_near(pxy.applyAffine(a).applyAffine(a_inv), pxy);
}

test "reflection" {
    try util.affine_expect_near(
        Affine.reflect(Point.ZERO, Vec2.new(1.0, 0.0)),
        Affine.new([_]f64{ 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 }),
    );
    try util.affine_expect_near(
        Affine.reflect(Point.ZERO, Vec2.new(0.0, 1.0)),
        Affine.new([_]f64{ -1.0, 0.0, 0.0, 1.0, 0.0, 0.0 }),
    );
    // y = x
    try util.affine_expect_near(
        Affine.reflect(Point.ZERO, Vec2.new(1.0, 1.0)),
        Affine.new([_]f64{ 0.0, 1.0, 1.0, 0.0, 0.0, 0.0 }),
    );

    // no translate
    const point = Point.new(0.0, 0.0);
    const vec = Vec2.new(1.0, 1.0);
    const map = Affine.reflect(point, vec);
    try util.expect_near(Point.new(0.0, 0.0).applyAffine(map), Point.new(0.0, 0.0));
    try util.expect_near(Point.new(1.0, 1.0).applyAffine(map), Point.new(1.0, 1.0));
    try util.expect_near(Point.new(1.0, 2.0).applyAffine(map), Point.new(2.0, 1.0));

    // with translate
    const point1 = Point.new(1.0, 0.0);
    const vec1 = Vec2.new(1.0, 1.0);
    const map1 = Affine.reflect(point1, vec1);
    try util.expect_near(Point.new(1.0, 0.0).applyAffine(map1), Point.new(1.0, 0.0));
    try util.expect_near(Point.new(2.0, 1.0).applyAffine(map1), Point.new(2.0, 1.0));
    try util.expect_near(Point.new(2.0, 2.0).applyAffine(map1), Point.new(3.0, 1.0));
}

test "aspect ratio" {
    const s = Size.new(1.0, 1.0);
    try std.testing.expect(@abs(s.aspectRatio() - 1.0) < 1e-6);
}

test "translate/scale" {
    const p = Point.new(3.0, 4.0);
    const ts = TranslateScale.new(Vec2.new(5.0, 6.0), 2.0);

    try util.expect_near(p.applyTranslateScale(ts), Point.new(11.0, 14.0));
}

test "conversions" {
    const p = Point.new(3.0, 4.0);
    const s = 2.0;
    const t = Vec2.new(5.0, 6.0);
    const ts = TranslateScale.new(t, s);

    const a = ts.toAffine();

    try util.expect_near(p.applyTranslateScale(ts), p.applyAffine(a));

    try util.expect_near(
        p.toVec2().mul(s).toPoint(),
        p.applyTranslateScale(TranslateScale.scale(s)),
    );
    try util.expect_near(
        p.toVec2().sum(t).toPoint(),
        p.applyTranslateScale(TranslateScale.translate(t)),
    );
}

test "inverse" {
    const p = Point.new(3.0, 4.0);
    const ts = TranslateScale.new(Vec2.new(5.0, 6.0), 2.0);

    try util.expect_near(p, p.applyTranslateScale(ts.apply(ts.inverse())));
    try util.expect_near(p, p.applyTranslateScale(ts.inverse().apply(ts)));
}
