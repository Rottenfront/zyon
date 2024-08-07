const mod = @import("module.zig");
const Rect = mod.Rect;
const RoundRectRadii = mod.RoundRectRadii;

pub const RoundRect = struct {
    rect: Rect,
    radii: RoundRectRadii,
};