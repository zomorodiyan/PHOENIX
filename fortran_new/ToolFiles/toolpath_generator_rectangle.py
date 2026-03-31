"""
Toolpath Generator for Rectangle Fill
Generates laser scan toolpath filling a rectangular region in X-Y.
Output format matches .crs files.

Two usage modes (all coordinates in meters, time in seconds):

  Usage 1 — Start corner (default):
    python3 toolpath_generator_rectangle.py \
      --start_x 0.0005 --start_y 0.0005 --start_z 0.0006975 \
      --size_x 0.003 --size_y 0.003 ...

  Usage 2 — Center of rectangle:
    python3 toolpath_generator_rectangle.py \
      --center_x 0.002 --center_y 0.002 --center_z 0.0006975 \
      --size_x 0.003 --size_y 0.003 ...

  If --center_x/--center_y are provided, start_x/start_y are computed as:
    start_x = center_x - size_x / 2
    start_y = center_y - size_y / 2
  --center_z sets start_z directly (Z is constant throughout the toolpath).

Input parameters:
  --start_x         X coordinate of rectangle starting corner (m) [default: 0.0005]
  --start_y         Y coordinate of rectangle starting corner (m) [default: 0.0005]
  --start_z         Z coordinate (constant layer height) (m) [default: 0.0007]
  --center_x        X coordinate of rectangle center (m) [overrides start_x]
  --center_y        Y coordinate of rectangle center (m) [overrides start_y]
  --center_z        Z coordinate (constant layer height) (m) [overrides start_z]
  --size_x          Rectangle extent in X direction (m) [default: 0.003]
  --size_y          Rectangle extent in Y direction (m) [default: 0.003]
  --scan_axis       Scan along "x" or "y" [default: x]
  --bidirectional   Alternating scan direction (default)
  --unidirectional  Same scan direction every track
  --hatch_spacing   Nominal hatch spacing (m), adjusted to fit rectangle exactly [default: 0.0001]
  --scan_speed      Laser scan speed (m/s) [default: 1.2]
  --turnaround_time Time with laser off between tracks (s) [default: 0.0005]
  --rotation_angle  CCW rotation of rectangle about its center (degrees) [default: 0]
  --output          Output .crs filename [default: toolpath.crs]
  --domain_x        Domain size in X for plot axis limits (m) [default: auto]
  --domain_y        Domain size in Y for plot axis limits (m) [default: auto]
  --test            Run all built-in test cases
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import argparse


def rotate_point(x, y, cx, cy, angle_deg):
    """Rotate (x, y) CCW by angle_deg about (cx, cy)."""
    angle_rad = np.radians(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    dx, dy = x - cx, y - cy
    return cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a


def rotate_points(xs, ys, cx, cy, angle_deg):
    """Rotate arrays of points CCW by angle_deg about (cx, cy)."""
    angle_rad = np.radians(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    dx, dy = xs - cx, ys - cy
    return cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a


def get_rotated_corners(start_x, start_y, size_x, size_y, cx, cy, angle_deg):
    """Get the 4 corners of the rotated rectangle, ordered as a polygon."""
    corners_x = np.array([start_x, start_x + size_x, start_x + size_x, start_x])
    corners_y = np.array([start_y, start_y, start_y + size_y, start_y + size_y])
    rx, ry = rotate_points(corners_x, corners_y, cx, cy, angle_deg)
    return rx, ry


def line_segment_intersect_horizontal(y_val, x1, y1, x2, y2):
    """Find x where horizontal line y=y_val intersects segment (x1,y1)-(x2,y2).
    Returns None if no intersection."""
    if abs(y2 - y1) < 1e-15:
        return None  # parallel
    t = (y_val - y1) / (y2 - y1)
    if -1e-10 <= t <= 1.0 + 1e-10:
        return x1 + t * (x2 - x1)
    return None


def line_segment_intersect_vertical(x_val, x1, y1, x2, y2):
    """Find y where vertical line x=x_val intersects segment (x1,y1)-(x2,y2).
    Returns None if no intersection."""
    if abs(x2 - x1) < 1e-15:
        return None  # parallel
    t = (x_val - x1) / (x2 - x1)
    if -1e-10 <= t <= 1.0 + 1e-10:
        return y1 + t * (y2 - y1)
    return None


def scan_line_intersections(corners_x, corners_y, scan_axis, scan_coord):
    """Find intersection points of an axis-aligned scan line with the rotated rectangle.

    For scan_axis='x': horizontal line at y=scan_coord, returns list of x values.
    For scan_axis='y': vertical line at x=scan_coord, returns list of y values.
    """
    n = len(corners_x)
    intersections = []
    for i in range(n):
        j = (i + 1) % n
        if scan_axis == "x":
            result = line_segment_intersect_horizontal(
                scan_coord, corners_x[i], corners_y[i], corners_x[j], corners_y[j]
            )
        else:
            result = line_segment_intersect_vertical(
                scan_coord, corners_x[i], corners_y[i], corners_x[j], corners_y[j]
            )
        if result is not None:
            intersections.append(result)
    # Remove near-duplicates (can happen at corners)
    if len(intersections) > 1:
        intersections.sort()
        unique = [intersections[0]]
        for v in intersections[1:]:
            if abs(v - unique[-1]) > 1e-12:
                unique.append(v)
        intersections = unique
    return intersections


def generate_toolpath(
    start_x, start_y, start_z,
    size_x, size_y,
    scan_axis="x",
    bidirectional=True,
    hatch_spacing=0.0001,
    scan_speed=1.0,
    turnaround_time=0.0005,
    output_filename="toolpath.crs",
    rotation_angle=0.0,
    domain_x=None,
    domain_y=None,
):
    """Generate a laser scan toolpath filling a rectangular region."""

    # Rectangle center (before rotation)
    cx = start_x + size_x / 2.0
    cy = start_y + size_y / 2.0

    # Get rotated corners
    corners_x, corners_y = get_rotated_corners(
        start_x, start_y, size_x, size_y, cx, cy, rotation_angle
    )

    # Determine perpendicular extent
    if scan_axis == "x":
        perp_min = np.min(corners_y)
        perp_max = np.max(corners_y)
    else:
        perp_min = np.min(corners_x)
        perp_max = np.max(corners_x)

    perp_extent = perp_max - perp_min

    # Compute adjusted hatch spacing
    n_tracks = max(round(perp_extent / hatch_spacing) + 1, 2)
    actual_hatch = perp_extent / (n_tracks - 1)

    # Starting point: bottommost vertex (x-scan) or leftmost vertex (y-scan)
    if scan_axis == "x":
        start_idx = np.argmin(corners_y)
        rot_start_x, rot_start_y = corners_x[start_idx], corners_y[start_idx]
    else:
        start_idx = np.argmin(corners_x)
        rot_start_x, rot_start_y = corners_x[start_idx], corners_y[start_idx]

    # Scan order: always start from perp_min (bottom for x-scan, left for y-scan)
    perp_positions = [perp_min + i * actual_hatch for i in range(n_tracks)]

    # Build waypoints
    waypoints = []  # (time, x, y, z, laser)
    time = 0.0

    # Line 1: all zeros (sentinel)
    waypoints.append((0.0, 0.0, 0.0, 0.0, 0))

    forward = True  # track scan direction for bidirectional
    first_track = True

    for i, perp_pos in enumerate(perp_positions):
        # Find intersections of this scan line with the rotated rectangle
        intersections = scan_line_intersections(corners_x, corners_y, scan_axis, perp_pos)

        if len(intersections) < 2:
            continue  # skip degenerate lines (at tips)

        scan_start = intersections[0]
        scan_end = intersections[-1]

        # Determine scan direction
        if bidirectional:
            if not forward:
                scan_start, scan_end = scan_end, scan_start
            forward = not forward
        # For unidirectional, always scan from scan_start to scan_end (low to high)

        # Compute track length and duration
        track_length = abs(scan_end - scan_start)
        if track_length < 1e-15:
            continue
        track_time = track_length / scan_speed

        # Track start/end coordinates
        if scan_axis == "x":
            track_start_x, track_start_y = scan_start, perp_pos
            track_end_x, track_end_y = scan_end, perp_pos
        else:
            track_start_x, track_start_y = perp_pos, scan_start
            track_end_x, track_end_y = perp_pos, scan_end

        if first_track:
            # Line 2: starting position at time=0, laser=0
            waypoints.append((0.0, track_start_x, track_start_y, start_z, 0))
            first_track = False
        else:
            # Turnaround: move to track start (laser off)
            time += turnaround_time
            waypoints.append((time, track_start_x, track_start_y, start_z, 0))

        # Track end (laser on)
        time += track_time
        waypoints.append((time, track_end_x, track_end_y, start_z, 1))

    # Final dwell (laser off)
    time += turnaround_time
    waypoints.append((time, waypoints[-1][1], waypoints[-1][2], start_z, 0))

    # Write .crs file
    output_dir = os.path.dirname(output_filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_filename, "w") as f:
        for wp in waypoints:
            f.write(
                f"  {wp[0]:14.8f}  {wp[1]:14.8f}  {wp[2]:14.8f}  {wp[3]:14.8f} {wp[4]} \n"
            )

    print(f"Written {output_filename}: {len(waypoints)} waypoints, {n_tracks} tracks, "
          f"actual hatch spacing = {actual_hatch*1e3:.4f} mm")

    # Starting point for visualization = line 2 (first track start)
    plot_start_x = waypoints[1][1]
    plot_start_y = waypoints[1][2]

    # Generate visualization
    plot_toolpath(
        waypoints, corners_x, corners_y,
        plot_start_x, plot_start_y,
        start_x, start_y, start_z,
        size_x, size_y, scan_axis, bidirectional,
        hatch_spacing, actual_hatch, scan_speed,
        turnaround_time, rotation_angle, n_tracks,
        output_filename, domain_x, domain_y,
    )

    return waypoints, corners_x, corners_y


def plot_toolpath(
    waypoints, corners_x, corners_y,
    rot_start_x, rot_start_y,
    start_x, start_y, start_z,
    size_x, size_y, scan_axis, bidirectional,
    hatch_spacing, actual_hatch, scan_speed,
    turnaround_time, rotation_angle, n_tracks,
    output_filename, domain_x=None, domain_y=None,
):
    """Plot the toolpath as arrows with rectangle outline."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    # Draw rotated rectangle outline
    poly_x = list(corners_x) + [corners_x[0]]
    poly_y = list(corners_y) + [corners_y[0]]
    ax.plot(
        [v * 1e3 for v in poly_x],
        [v * 1e3 for v in poly_y],
        "b-", linewidth=1.5,
    )

    # Draw scan tracks as arrows
    i = 0
    while i < len(waypoints) - 1:
        if waypoints[i + 1][4] == 1:
            x0, y0 = waypoints[i][1] * 1e3, waypoints[i][2] * 1e3
            x1, y1 = waypoints[i + 1][1] * 1e3, waypoints[i + 1][2] * 1e3
            ax.annotate(
                "", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color="red", lw=1.2),
            )
        i += 1

    # Mark starting point (green dot only)
    ax.plot(rot_start_x * 1e3, rot_start_y * 1e3, "go", markersize=8, zorder=5)

    # Axis labels
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # Set axis limits to domain size if provided
    if domain_x is not None:
        ax.set_xlim(0, domain_x * 1e3)
    if domain_y is not None:
        ax.set_ylim(0, domain_y * 1e3)

    # Parameter text block at the top (3 lines)
    bidir_str = "bidirectional" if bidirectional else "unidirectional"
    param_text = (
        f"start=({start_x*1e3:.2f}, {start_y*1e3:.2f}, {start_z*1e3:.2f}) mm    "
        f"size_x={size_x*1e3:.2f} mm    size_y={size_y*1e3:.2f} mm\n"
        f"scan_axis={scan_axis}    {bidir_str}    "
        f"rotation={rotation_angle:.1f} deg    n_tracks={n_tracks}\n"
        f"hatch={hatch_spacing*1e3:.4f} mm (actual={actual_hatch*1e3:.4f} mm)    "
        f"speed={scan_speed:.2f} m/s    turnaround={turnaround_time*1e3:.2f} ms"
    )
    # Compute layout: equal spacing between text bottom and plot top as between
    # text top and figure top. 3 lines of text at fontsize 10 ~ 0.06 of fig height.
    text_height = 0.06
    margin = 0.02  # gap above text and between text and plot
    top = 1.0 - margin - text_height - margin  # axes top position
    fig.subplots_adjust(top=top)
    fig.suptitle(param_text, fontsize=10, family="monospace",
                 y=1.0 - margin, va="top")

    # Save
    png_filename = os.path.splitext(output_filename)[0] + ".png"
    fig.savefig(png_filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Written {png_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rectangle toolpath generator")
    parser.add_argument("--start_x", type=float, default=0.0005)
    parser.add_argument("--start_y", type=float, default=0.0005)
    parser.add_argument("--start_z", type=float, default=0.0007)
    parser.add_argument("--center_x", type=float, default=None,
                        help="X center of rectangle (overrides start_x)")
    parser.add_argument("--center_y", type=float, default=None,
                        help="Y center of rectangle (overrides start_y)")
    parser.add_argument("--center_z", type=float, default=None,
                        help="Z coordinate (overrides start_z)")
    parser.add_argument("--size_x", type=float, default=0.003)
    parser.add_argument("--size_y", type=float, default=0.003)
    parser.add_argument("--scan_axis", type=str, default="x", choices=["x", "y"])
    parser.add_argument("--bidirectional", action="store_true", default=True)
    parser.add_argument("--unidirectional", action="store_true")
    parser.add_argument("--hatch_spacing", type=float, default=0.0001)
    parser.add_argument("--scan_speed", type=float, default=1.2)
    parser.add_argument("--turnaround_time", type=float, default=0.0005)
    parser.add_argument("--rotation_angle", type=float, default=0.0)
    parser.add_argument("--output", type=str, default="toolpath.crs")
    parser.add_argument("--domain_x", type=float, default=None,
                        help="Domain size in X for plot limits (m)")
    parser.add_argument("--domain_y", type=float, default=None,
                        help="Domain size in Y for plot limits (m)")
    args = parser.parse_args()

    # Convert center coordinates to start coordinates if provided
    sx = args.center_x - args.size_x / 2.0 if args.center_x is not None else args.start_x
    sy = args.center_y - args.size_y / 2.0 if args.center_y is not None else args.start_y
    sz = args.center_z if args.center_z is not None else args.start_z

    bidir = not args.unidirectional
    generate_toolpath(
        start_x=sx, start_y=sy, start_z=sz,
        size_x=args.size_x, size_y=args.size_y,
        scan_axis=args.scan_axis, bidirectional=bidir,
        hatch_spacing=args.hatch_spacing, scan_speed=args.scan_speed,
        turnaround_time=args.turnaround_time,
        output_filename=args.output, rotation_angle=args.rotation_angle,
        domain_x=args.domain_x, domain_y=args.domain_y,
    )
