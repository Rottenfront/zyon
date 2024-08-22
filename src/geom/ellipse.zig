const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Affine = mod.Affine;

pub const Ellipse = struct {
    inner: Affine,
};
