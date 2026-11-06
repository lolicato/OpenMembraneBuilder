"""The schematic system view of the protein insertion page.

It draws the box, the membrane and every protein as a simple picture (SVG), seen from the
side (x–z or y–z) or from the top (x–y). Each protein is drawn as the outline of its own
beads, placed and rotated exactly the way COBY will place it:

1. the protein is centred on the middle of its beads (centre of geometry),
2. it is rotated by the angles you entered (first about x, then y, then z, in the protein's
   own frame – COBY's convention),
3. it is moved to the position you entered, measured from the centre of the box.

The membrane is drawn with a typical Martini thickness. "Height above Membrane" is therefore
only approximate in the picture; the build measures the real membrane surface.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.spatial.transform import Rotation as R

VIEWS = {"x–z": (0, 2), "y–z": (1, 2), "x–y": (0, 1)}

# Typical Martini bilayer (nm from the midplane): outer edge of the head groups and phosphate
# planes. Used only for drawing and for the approximate "Height above Membrane" position.
MEMBRANE_HALF = 2.2
PHOSPHATE = 2.0


def coby_rotation(rx: float, ry: float, rz: float) -> np.ndarray:
    """The rotation COBY applies for the angles rx, ry, rz.

    Args:
        rx: Rotation about x, in degrees.
        ry: Rotation about y, in degrees.
        rz: Rotation about z, in degrees.

    Returns:
        A 3 × 3 rotation matrix (Rx·Ry·Rz).
    """
    return R.from_euler("XYZ", [rx, ry, rz], degrees=True).as_matrix()


def box_extent(box_type: str, bx: float, by: float, bz: float) -> tuple[float, float, float]:
    """The size of the box along x, y and z, in nm.

    Hexagonal boxes are not rectangular; the picture uses the rectangle around them.

    Args:
        box_type: `"rectangular"`, `"hexagonal"` or `"skewed_hexagonal"`.
        bx: Box length x, in nm.
        by: Box length y, in nm (rectangular boxes only).
        bz: Box height z, in nm (not used by skewed hexagonal boxes).

    Returns:
        `(x, y, z)` in nm.
    """
    if box_type == "hexagonal":
        return bx, bx * math.sqrt(3) / 2, bz
    if box_type == "skewed_hexagonal":
        return bx, bx * math.sqrt(3) / 2, bx * math.sqrt(2 / 3)
    return bx, by, bz


def place(coords_A: np.ndarray, cx: float, cy: float, cz: float, rx: float, ry: float, rz: float) -> np.ndarray:
    """Bead positions of a protein after COBY places it.

    Args:
        coords_A: Bead positions as read from the file, in Å.
        cx: Position along x, in nm from the box centre.
        cy: Position along y, in nm from the box centre.
        cz: Position along z, in nm from the box centre.
        rx: Rotation about x, in degrees.
        ry: Rotation about y, in degrees.
        rz: Rotation about z, in degrees.

    Returns:
        Bead positions in nm, measured from the box centre.
    """
    c = (coords_A - coords_A.mean(0)) @ coby_rotation(rx, ry, rz).T / 10.0
    return c + np.array([cx, cy, cz])


def height_above_membrane_cz(coords_A, rx, ry, rz, distance_to_mem: float) -> float:
    """The z position that puts the lowest bead a given distance above the membrane.

    Uses the typical membrane thickness of the picture (`MEMBRANE_HALF`).

    Args:
        coords_A: Bead positions as read from the file, in Å.
        rx: Rotation about x, in degrees.
        ry: Rotation about y, in degrees.
        rz: Rotation about z, in degrees.
        distance_to_mem: Gap between membrane surface and protein, in nm.

    Returns:
        The z position of the protein centre, in nm from the box centre.
    """
    c = (coords_A - coords_A.mean(0)) @ coby_rotation(rx, ry, rz).T / 10.0
    return MEMBRANE_HALF + distance_to_mem - c[:, 2].min()


def half_extent_x(coords_A, rx, ry, rz) -> float:
    """Half of the width of the rotated protein along x, in nm.

    Used when two proteins are placed by the gap between them.

    Args:
        coords_A: Bead positions as read from the file, in Å.
        rx: Rotation about x, in degrees.
        ry: Rotation about y, in degrees.
        rz: Rotation about z, in degrees.

    Returns:
        Half the x width, in nm.
    """
    c = (coords_A - coords_A.mean(0)) @ coby_rotation(rx, ry, rz).T / 10.0
    return (c[:, 0].max() - c[:, 0].min()) / 2.0


def make_grid(boxx: float, boxy: float, n_copies: int) -> list[tuple[float, float]]:
    """The grid on which copies "moving together" are placed.

    The box is divided into n × n equal cells (n = smallest whole number with
    n² ≥ number of copies); the copies go into the centres of randomly chosen cells.

    Args:
        boxx: Box length x, in nm.
        boxy: Box length y, in nm.
        n_copies: Number of copies.

    Returns:
        The cell centres `(x, y)` in nm from the box centre.
    """
    n = math.ceil(math.sqrt(n_copies))
    dx, dy = boxx / n, boxy / n
    return [((i + 0.5) * dx - boxx / 2, (j + 0.5) * dy - boxy / 2) for i in range(n) for j in range(n)]


def _hull2d(pts: np.ndarray) -> np.ndarray:
    """The outline (convex hull) of a set of points in the plane, in drawing order.

    Args:
        pts: Points `(h, v)`.

    Returns:
        The corner points of the outline.
    """
    pts = np.unique(np.round(np.asarray(pts, float), 3), axis=0)
    if len(pts) <= 2:
        return pts

    def cross(o, a, b):
        """Which side of the line o→a the point b lies on (positive: left)."""
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for q in map(tuple, pts):
        while len(lower) >= 2 and cross(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    for q in map(tuple, pts[::-1]):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    return np.array(lower[:-1] + upper[:-1])


def system_view_svg(extent: tuple[float, float, float], proteins: list[dict], membrane: bool, view: str = "x–z") -> str:
    """The picture of the system, seen along one axis, as SVG text.

    Args:
        extent: Box size `(x, y, z)` in nm.
        proteins: One entry per protein: `name`, `copies` (bead positions of every copy,
            in nm from the box centre) and `dashed` (`True` when the position or rotation
            is random, so the picture shows only one possibility).
        membrane: Draw the membrane in the middle of the box.
        view: `"x–z"`, `"y–z"` or `"x–y"`.

    Returns:
        SVG text that can be put straight into the page.
    """
    ha, va = VIEWS.get(view, (0, 2))
    bh, bv = extent[ha], extent[va]
    zview = va == 2
    W, H, pad = 420, 320, 14
    sc = min((W - 2 * pad) / bh, (H - 2 * pad) / bv)
    oh, ov = (W - bh * sc) / 2, (H - bv * sc) / 2
    X = lambda h: oh + (h + bh / 2) * sc  # noqa: E731
    Z = lambda v: ov + (bv / 2 - v) * sc  # noqa: E731
    lip = "#5a7a1a"
    g = [f'<rect x="{oh}" y="{ov}" width="{bh * sc}" height="{bv * sc}" fill="#2f6fb4" fill-opacity=".12" stroke="#9aa6b0"/>']
    if membrane:
        if zview:
            g.append(f'<rect x="{oh}" y="{Z(MEMBRANE_HALF)}" width="{bh * sc}" height="{2 * MEMBRANE_HALF * sc}" fill="{lip}" fill-opacity=".35"/>')
            for pz in (PHOSPHATE, -PHOSPHATE):
                g.append(f'<line x1="{oh}" x2="{oh + bh * sc}" y1="{Z(pz)}" y2="{Z(pz)}" stroke="{lip}" stroke-dasharray="4 3"/>')
            g.append(f'<text x="{oh + bh * sc - 3}" y="{Z(MEMBRANE_HALF / 2) + 3}" font-size="9" text-anchor="end" fill="#2f4a08">upper</text>')
            g.append(f'<text x="{oh + bh * sc - 3}" y="{Z(-MEMBRANE_HALF / 2) + 3}" font-size="9" text-anchor="end" fill="#2f4a08">lower</text>')
        else:
            g.append(f'<rect x="{oh}" y="{ov}" width="{bh * sc}" height="{bv * sc}" fill="{lip}" fill-opacity=".18"/>')
            g.append(f'<text x="{oh + 4}" y="{ov + 12}" font-size="10" fill="#5d6b77">membrane (z = 0, whole plane)</text>')
    for p in proteins:
        dashed = ' stroke-dasharray="3 2"' if p.get("dashed") else ""
        for i, c in enumerate(p["copies"]):
            hull = _hull2d(c[:, [ha, va]])
            pts = " ".join(f"{X(h):.1f},{Z(v):.1f}" for h, v in hull)
            g.append(f'<polygon points="{pts}" fill="none" stroke="#ffffff" stroke-width="3.5" stroke-opacity=".9"/>')
            g.append(f'<polygon points="{pts}" fill="#b5541a" fill-opacity="{".4" if dashed else ".6"}" stroke="#7a3410" stroke-width="1.5"{dashed}/>')
            if i == 0:
                n = len(p["copies"])
                g.append(f'<text x="{X(hull[:, 0].max()) + 3}" y="{Z(c[:, va].mean()) + 4}" font-size="10" fill="#b5541a">'
                         f'{p["name"]}{" ×" + str(n) if n > 1 else ""}</text>')
    g.append(f'<text x="{oh + bh * sc - 4}" y="{ov + bv * sc - 4}" font-size="10" text-anchor="end" fill="#5d6b77">{bh:g} × {bv:g} nm ({view})</text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="100%" '
            f'style="background:#f5f7f9;border:1px solid #d3d9df;border-radius:4px;max-width:720px">' + "".join(g) + "</svg>")
