#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import copy
import os

# =========================
# USER SETTINGS
# =========================

INPUT_FILE = "/home/kali/Software/Newtron source git/NeutronSource/Geometry/gdml files/ElectronicsX0Y0.gdml"   # your original GDML file
OUTPUT_DIR = "gdml_out"

# Initial position from your GDML
X_START = 9.75
Y_START = -109.5

# Step sizes (mm)
DX = 9.0
DY = 10.0

# Number of steps
NX = 20+1   # x steps
NY = 10+1   # y steps

# =========================
# FUNCTIONS
# =========================

def find_electronics_position(root):
    """
    Find the <position> node inside:
    <physvol name="PV_LV_Electronics">
    """
    for physvol in root.iter("physvol"):
        if physvol.get("name") == "PV_LV_Electronics":
            for child in physvol:
                if child.tag == "position":
                    return child
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

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load original GDML
    try:
        tree = ET.parse(INPUT_FILE)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    root = tree.getroot()

    # Check that Electronics exists
    test_pos = find_electronics_position(root)
    if test_pos is None:
        raise RuntimeError("Could not find PV_LV_Electronics position in GDML!")

    total_files = 0

    # Loop over shifts
    for ix in range(NX):
        for iy in range(NY):

            # Deep copy full GDML structure
            new_root = copy.deepcopy(root)

            # Find position node in copied tree
            pos = find_electronics_position(new_root)

            # Compute new coordinates
            new_x = X_START + ix * DX
            new_y = Y_START + iy * DY

            # Apply shifts
            pos.set("x", f"{new_x:.6f}")
            pos.set("y", f"{new_y:.6f}")

            # Optional: keep z unchanged
            # pos.set("z", pos.get("z"))

            # Pretty formatting
            indent(new_root)

            # Output filename
            filename = os.path.join(
                OUTPUT_DIR,
                f"Electronics_x{ix:02d}_y{iy:02d}.gdml"
            )

            # Write file
            ET.ElementTree(new_root).write(
                filename,
                encoding="utf-8",
                xml_declaration=True
            )

            total_files += 1

            print(f"Created: {filename} (x={new_x}, y={new_y})")

    print(f"\n✅ Done! Generated {total_files} GDML files.")


# =========================

if __name__ == "__main__":
    main()