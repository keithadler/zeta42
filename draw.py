"""Render unit-distance graphs to SVG, using the exact coordinates.

    ./.venv/bin/python draw.py

Writes spindle.svg and degrey-L.svg next to this file.  The geometry is the
real thing -- same points the SAT solver was given, just projected to float
at the last moment for drawing.
"""

import os

from field import get_field
from udg import moser_spindle
from colour import k_colourable
from degrey import build_L

HERE = os.path.dirname(os.path.abspath(__file__))

PALETTE = ["#d94a4a", "#3b7dd8", "#e8a33d", "#3aa66d", "#8b5cf6", "#14b8a6"]
BG = "#fbfaf8"
INK = "#2a2724"
EDGE = "#b9b2a8"


def render(pts, edges, colouring, width=520, height=520, pad=46,
           r=7, title="", subtitle=""):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    spanx = max(xs) - min(xs) or 1
    spany = max(ys) - min(ys) or 1
    top = 0 if not title else 54
    scale = min((width - 2 * pad) / spanx, (height - top - 2 * pad) / spany)
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2

    def place(p):
        x = width / 2 + (p[0] - cx) * scale
        y = (height + top) / 2 - (p[1] - cy) * scale     # flip: SVG y grows down
        return x, y

    xy = [place(p) for p in pts]
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
           f'width="{width}" height="{height}" font-family="ui-sans-serif, system-ui, sans-serif">',
           f'<rect width="{width}" height="{height}" fill="{BG}"/>']
    if title:
        out.append(f'<text x="{pad}" y="30" font-size="17" font-weight="600" '
                   f'fill="{INK}">{title}</text>')
    if subtitle:
        out.append(f'<text x="{pad}" y="48" font-size="12.5" fill="#7c746a">{subtitle}</text>')

    out.append(f'<g stroke="{EDGE}" stroke-width="1.1" stroke-linecap="round">')
    for i, j in edges:
        x1, y1 = xy[i]
        x2, y2 = xy[j]
        out.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}"/>')
    out.append('</g>')

    for k, (x, y) in enumerate(xy):
        fill = PALETTE[colouring[k] % len(PALETTE)] if colouring else "#8a827a"
        out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" fill="{fill}" '
                   f'stroke="{BG}" stroke-width="1.6"/>')
    out.append('</svg>')
    return "\n".join(out)


def main():
    # ---- Moser spindle: every line is exactly one unit long ----
    g = moser_spindle()
    _, col = k_colourable(g, 4)
    pts = [g.fld.pfloat(p) for p in g.points]
    svg = render(pts, g.edges, col, width=460, height=440, r=11,
                 title="Moser spindle — 7 points, 11 unit-length edges",
                 subtitle="4 colours suffice; 3 provably do not. This is why the plane needs ≥ 4.")
    open(os.path.join(HERE, "spindle.svg"), "w").write(svg)

    # ---- the same spindle, best possible 3-colouring ----
    ok3, col3 = k_colourable(g, 3)
    if not ok3:
        # no proper 3-colouring exists; show a greedy one and mark the clash
        greedy = [None] * len(g)
        adj = g.adjacency()
        for v in range(len(g)):
            used = {greedy[u] for u in adj[v] if greedy[u] is not None}
            greedy[v] = next(c for c in range(3) if c not in used) if len(used) < 3 else 0
        bad = [(i, j) for i, j in g.edges if greedy[i] == greedy[j]]
        svg3 = render(pts, g.edges, greedy, width=460, height=440, r=11,
                      title="The same 7 points with only 3 colours",
                      subtitle=f"{len(bad)} edge(s) must join matching colours — no arrangement avoids it.")
        open(os.path.join(HERE, "spindle-3.svg"), "w").write(svg3)

    # ---- de Grey's L ----
    f = get_field((3, 5, 7))
    L = build_L(f)
    lg = L.graph()
    _, lcol = k_colourable(lg, 4)
    lpts = [f.pfloat(p) for p in lg.points]
    svg = render(lpts, lg.edges, lcol, width=640, height=640, r=4.6,
                 title=f"de Grey's L — {len(lg)} points, {len(lg.edges)} unit-length edges",
                 subtitle="52 overlapping hexagons. Every edge here is exactly 1 unit long.")
    open(os.path.join(HERE, "degrey-L.svg"), "w").write(svg)

    # ---- de Grey's G, the 1581-vertex non-4-colourable graph ----
    from degrey_g import build_G
    G = build_G()
    _, gcol = k_colourable(G, 5)
    gpts = [G.fld.pfloat(p) for p in G.points]
    svg = render(gpts, G.edges, gcol, width=760, height=760, r=2.2,
                 title=f"de Grey's G — {len(G)} points, {len(G.edges)} unit-length edges",
                 subtitle="Not 4-colourable (verified). 5 colours suffice — hence the plane needs ≥ 5.")
    open(os.path.join(HERE, "degrey-G.svg"), "w").write(svg)

    for n in ("spindle.svg", "spindle-3.svg", "degrey-L.svg", "degrey-G.svg"):
        p = os.path.join(HERE, n)
        if os.path.exists(p):
            print(f"wrote {p} ({os.path.getsize(p)} bytes)")


if __name__ == "__main__":
    main()
