const std = @import("std");
const geom = @import("geom/module.zig");
const Affine = geom.Affine;
const Vec2 = geom.Vec2;

pub fn main() !void {
    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});

    // stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Run `zig build test` to run the tests.\n", .{});

    const aff1 = Affine.scale(2.0);
    const aff2 = Affine.translate(Vec2.new(100.0, 100.0));

    std.debug.print("aff1.apply(aff2).translate_x = {d}\n", .{aff1.apply(aff2).a[5]});
    std.debug.print("aff2.apply(aff1).translate_x = {d}\n", .{aff2.apply(aff1).a[5]});
    std.debug.print("aff2.pre_scale(2.0).translate_x = {d}\n", .{aff2.pre_scale(2.0).a[5]});
    std.debug.print("aff1.then_translate(Vec2.new(100.0, 100.0)).translate_x = {d}\n", .{aff1.then_translate(Vec2.new(100.0, 100.0)).a[5]});

    try bw.flush(); // don't forget to flush!
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
