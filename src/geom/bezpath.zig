const std = @import("std");
const math = std.math;
const mod = @import("module.zig");

const Point = mod.Point;

pub const CubicBez = @import("bezpath/cubic.zig").CubicBez;
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
    pub const PathList: type = std.ArrayList(Element);
    pub const allocator = std.heap.page_allocator;

    data: PathList,

    /// Create a new path.
    pub fn new() @This() {
        return @This(){
            .data = PathList.init(@This().allocator),
        };
    }

    /// Removes the last [`PathEl`] from the path and returns it, or `None` if the path is empty.
    pub fn pop(self: *@This()) ?Element {
        return self.data.popOrNull();
    }

    /// Push a generic path element onto the path.
    pub fn push(self: *@This(), el: Element) @This().allocator.Error!void {
        try self.data.append(el);
        std.debug.assert(match: {
            switch (self.data.items[0]) {
                .move_to => |_| break :match true,
                else => break :match false,
            }
        });
    }

    /// Push a "move to" element onto the path.
    pub fn move_to(self: *@This(), p: Point) @This().allocator.Error!void {
        return self.push(Element{ .move_to = .{ .p = p } });
    }

    /// Push a "line to" element onto the path.
    ///
    /// Will panic with a debug assert when the path is empty and there is no
    /// "move to" element on the path.
    ///
    /// If `line_to` is called immediately after `close_path` then the current
    /// subpath starts at the initial point of the previous subpath.
    pub fn line_to(self: *@This(), end: Point) @This().allocator.Error!void {
        return self.push(Element{ .line_to = .{ .end = end } });
    }

    /// Push a "quad to" element onto the path.
    ///
    /// Will panic with a debug assert when the path is empty and there is no
    /// "move to" element on the path.
    ///
    /// If `quad_to` is called immediately after `close_path` then the current
    /// subpath starts at the initial point of the previous subpath.
    pub fn quad_to(
        self: *@This(),
        p0: Point,
        end: Point,
    ) @This().allocator.Error!void {
        return self.push(Element{ .quad_to = .{ .p0 = p0, .end = end } });
    }

    /// Push a "curve to" element onto the path.
    ///
    /// Will panic with a debug assert when the path is empty and there is no
    /// "move to" element on the path.
    ///
    /// If `curve_to` is called immediately after `close_path` then the current
    /// subpath starts at the initial point of the previous subpath.
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

    /// Push a "close path" element onto the path.
    ///
    /// Will panic with a debug assert when the path is empty and there is no
    /// "move to" element on the path.
    pub fn close_path(self: *@This()) @This().allocator.Error!void {
        return self.push(Element{.close_path});
    }

    /// Shorten the path, keeping the first `len` elements.
    pub fn truncate(self: *@This(), len: usize) @This().allocator.Error!void {
        if (self.data.items.len > len) {
            return self.data.resize(len);
        }
    }

    pub fn flatten(self: *@This(), tolerance: f64, ST: type, start_context: *ST, callback: fn (*ST, Element) void) void {
        const sqrt_tol = @sqrt(tolerance);
        var last_pt: ?Point = null;
        // var quad_buf = std.ArrayList(struct {
        //     quad: QuadBez,
        //     params: QuadBez.FlattenParams,
        // }).init(allocator);

        for (self.data.items) |el| {
            switch (el) {
                .move_to => |data| {
                    last_pt = data.p;
                    callback(start_context, Element{
                        .move_to = struct { .p = data.p },
                    });
                },
                .line_to => |data| {
                    last_pt = data.p;
                    callback(start_context, Element{
                        .line_to = struct { .p = data.p },
                    });
                },
                .quad_to => |data| {
                    if (last_pt != null) {
                        const q = QuadBez.new(last_pt.?, data.p0, data.end);
                        const params = q.estimateSubdiv(sqrt_tol);
                        const n = @max(@as(usize, math.ceil(0.5 * params.val / sqrt_tol)), 1);
                        const step = 1.0 / @as(f64, n);
                        for (1..n) |i| {
                            const u = @as(f64, i) * step;
                            const t = q.determineSubdivT(&params, u);
                            const p = q.eval(t);
                            callback(start_context, Element{
                                .line_to = struct { .p = p },
                            });
                        }
                        callback(start_context, Element{
                            .line_to = struct { .p = data.end },
                        });
                    }
                    last_pt = data.end;
                },
            }
        }
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
