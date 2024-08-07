const mod = @import("module.zig");
const Vec2 = mod.Vec2;

pub const TranslateScale = struct {
    translation: Vec2,
    scale: f64,
};
