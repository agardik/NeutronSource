#!/usr/bin/env python3
"""
Generate the electronics-strip variants (Top_Strip_1 ... Top_Strip_N) from one base GDML.

For every strip it
  * moves  <physvol name="PV_LV_Electronics">  to the strip's y centre, and
  * sets the y-size of the box used by LV_Electronics to the strip's width,
and writes one GDML file per strip. Everything else (tungsten Lid001, lid, detector...)
is copied unchanged from the base file.

Strips are stacked contiguously in -y, starting at Y_TOP_EDGE:
    centre_i = Y_TOP_EDGE - sum(widths of previous strips) - width_i / 2
With the defaults this gives y = -17.5, -47.5, -77.5, -104.5 (the last strip is 24 mm wide).

Usage:
    python3 gdml_strip_maker.py                       # uses the settings below
    python3 gdml_strip_maker.py base.gdml -o out_dir  # override input / output dir
"""

import argparse
import copy
import os
import xml.etree.ElementTree as ET

# =========================
# USER SETTINGS
# =========================

INPUT_FILE = "/home/kali/Software/Newtron source git/NeutronSource/Geometry/gdml_tungsten_strips/tungsten_shield_geometry_1cm_tungsten.gdml"   # base GDML (tungsten 1 cm version)
OUTPUT_DIR = "/home/kali/Software/Newtron source git/NeutronSource/Geometry/gdml_tungsten_strips/gdml_out"
OUTPUT_PREFIX = "Electronics_Top_Strip"                      # -> Electronics_Top_Strip_1.gdml, ...

# y position (mm) of the upper edge of strip 1 (strip 1 centre -17.5 mm, width 30 mm -> -2.5 mm)
Y_TOP_EDGE = -2.5

# y-size (mm) of each strip, in order; strips are stacked towards -y
STRIP_WIDTHS = [30.0, 30.0, 30.0, 24.0]

# Names used in the GDML
ELECTRONICS_PV = "PV_LV_Electronics"   # physvol that positions the strip
ELECTRONICS_LV = "LV_Electronics"      # logical volume whose solid is the strip box

# Keep the "xsi:" prefix on the schema attribute instead of ElementTree's "ns0:"
ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")

# =========================
# FUNCTIONS
# =========================


def strip_centres(widths, y_top):
    """Return the y centre of each strip when stacked contiguously from y_top towards -y."""
    centres = []
    edge = y_top
    for w in widths:
        centres.append(edge - w / 2.0)
        edge -= w
    return centres


def find_electronics_position(root):
    """
    Find the <position> node inside:
    <physvol name="PV_LV_Electronics">
    """
    for physvol in root.iter("physvol"):
        if physvol.get("name") == ELECTRONICS_PV:
            for child in physvol:
                if child.tag == "position":
                    return child
    return None


def find_electronics_box(root):
    """
    Find the <box> solid referenced by <volume name="LV_Electronics">
    (its name changes between files: Cube1, Cube015, Cube016, ...).
    """
    for volume in root.iter("volume"):
        if volume.get("name") == ELECTRONICS_LV:
            solidref = volume.find("solidref")
            if solidref is None:
                return None
            solid_name = solidref.get("ref")
            for box in root.iter("box"):
                if box.get("name") == solid_name:
                    return box
    return None


def indent(elem, level=0):
    """
    Pretty-print XML (optional but nice)
    """
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for child in elem:
            indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = i
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = i


# =========================
# MAIN SCRIPT
# =========================

def main():
    parser = argparse.ArgumentParser(description="Generate electronics-strip GDML variants.")
    parser.add_argument("input", nargs="?", default=INPUT_FILE, help="base GDML file")
    parser.add_argument("-o", "--output-dir", default=OUTPUT_DIR, help="output directory")
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load base GDML
    try:
        tree = ET.parse(args.input)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    root = tree.getroot()

    # Check that the electronics position and box exist
    if find_electronics_position(root) is None:
        raise RuntimeError(f"Could not find {ELECTRONICS_PV} position in GDML!")
    if find_electronics_box(root) is None:
        raise RuntimeError(f"Could not find the box solid of {ELECTRONICS_LV} in GDML!")

    centres = strip_centres(STRIP_WIDTHS, Y_TOP_EDGE)
    total_files = 0

    for n, (width, y_centre) in enumerate(zip(STRIP_WIDTHS, centres), start=1):

        # Deep copy full GDML structure
        new_root = copy.deepcopy(root)

        # Find the nodes in the copied tree
        pos = find_electronics_position(new_root)
        box = find_electronics_box(new_root)

        # Apply: move the strip and set its width (x and z are left unchanged)
        pos.set("y", f"{y_centre:.6f}")
        box.set("y", f"{width:.6f}")

        # Pretty formatting
        indent(new_root)

        # Output filename
        filename = os.path.join(args.output_dir, f"{OUTPUT_PREFIX}_{n}.gdml")

        # Write file
        ET.ElementTree(new_root).write(
            filename,
            encoding="utf-8",
            xml_declaration=True
        )

        total_files += 1

        print(f"Created: {filename} "
              f"(y={y_centre:.3f}, width={width:.1f}, "
              f"covers y={y_centre - width / 2:.2f} .. {y_centre + width / 2:.2f})")

    print(f"\n✅ Done! Generated {total_files} GDML files.")


# =========================

if __name__ == "__main__":
    main()