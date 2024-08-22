const std = @import("std");
const math = std.math;
const mod = @import("module.zig");

const Point = mod.Point;

pub const CubicBez = @import("bezpath/cubic.zig").CubicBez;
pub const CuspType = @import("bezpath/cubic.zig").CuspType;
pub const Line = @import("bezpath/line.zig").Line;
pub const QuadBez = @import("bezpath/quad.zig").QuadBez;

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
    pub const PathList: type = std.ArrayList(PathEl);
    pub const allocator = std.heap.page_allocator;

    data: PathList,

    pub fn new() @This() {
        return @This(){
            .data = PathList.init(@This().allocator),
        };
    }

    pub fn pop(self: *@This()) ?PathEl {
        return self.data.popOrNull();
    }

    pub fn push(self: *@This(), el: PathEl) @This().allocator.Error!void {
        try self.data.append(el);
        std.debug.assert(match: {
            switch (self.data.items[0]) {
                .move_to => |_| break :match true,
                else => break :match false,
            }
        });
    }

    pub fn move_to(self: *@This(), p: Point) @This().allocator.Error!void {
        return self.push(PathEl{ .move_to = .{ .p = p } });
    }

    pub fn line_to(self: *@This(), end: Point) @This().allocator.Error!void {
        return self.push(PathEl{ .line_to = .{ .end = end } });
    }

    pub fn quad_to(
        self: *@This(),
        p0: Point,
        end: Point,
    ) @This().allocator.Error!void {
        return self.push(PathEl{ .quad_to = .{ .p0 = p0, .end = end } });
    }

    pub fn curve_to(
        self: *@This(),
        p0: Point,
        p1: Point,
        end: Point,
    ) @This().allocator.Error!void {
        return self.push(PathEl{ .curve_to = .{
            .p0 = p0,
            .p1 = p1,
            .end = end,
        } });
    }

    pub fn close_path(self: *@This()) @This().allocator.Error!void {
        return self.push(PathEl{.close_path});
    }
};

/// The element of a Bézier path.
///
/// A valid path has `MoveTo` at the beginning of each subpath.
pub const PathEl = union(enum) {
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
pub const PathSeg = union(enum) {
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
