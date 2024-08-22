const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const CubicBez = mod.CubicBez;
const QuadBez = mod.QuadBez;

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
