const std = @import("std");
const math = std.math;
const mod = @import("module.zig");
const Rect = mod.Rect;

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
