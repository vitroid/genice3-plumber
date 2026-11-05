# coding: utf-8
"""
GenIce3 exporter plugin: ice hydrogen-bond network as OpenSCAD pipework.

Usage:
    genice3 1h -e plumber > ice.scad
    genice3 1h --rep 2 2 2 -e "plumber :scale 40 :R 0.04"
"""

from __future__ import annotations

import sys
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import numpy as np

from genice3_plumber.pipework import bondfunc, pipework, torusfunc

desc = {
    "ref": {"OpenSCAD": "https://openscad.org"},
    "brief": "OpenSCAD pipework of the ice hydrogen-bond network.",
    "usage": """
Usage: genice3 UNITCELL -e plumber [:option value ...]

options:
    scale=50   Scale factor (nm to OpenSCAD units), applied to coordinates and R.
    R=0.05     Straight-pipe radius (nm). Joints are larger (tee 1.5R, collar 1.8R).
    fn=20      OpenSCAD $fn for pipes, elbows, and tees.
""",
}

format_desc = {
    "aliases": ["plumber"],
    "application": "[OpenSCAD](https://openscad.org)",
    "extension": ".scad",
    "water": "Center of mass",
    "solute": "none",
    "hb": "o",
    "remarks": (
        "Pipes with tee, elbow, and plug joints. "
        "Edges whose minimum image crosses a cell face are omitted. "
        "Nodes on a face are kept and are not duplicated onto the opposite face. "
        "A node left with no bonds is omitted."
    ),
    "suboptions": "scale, R, fn.",
}


def _scalar(v: Any) -> Any:
    """ネストしたリスト1要素ならスカラーに。"""
    while isinstance(v, (list, tuple)) and len(v) == 1:
        v = v[0]
    return v


def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """scale / R / fn を処理し、残りは unprocessed として返す。"""
    processed: Dict[str, Any] = {}
    unprocessed = dict(options)
    for key in ("scale", "R", "fn"):
        if key not in unprocessed:
            continue
        raw = _scalar(unprocessed.pop(key))
        if key == "fn":
            processed[key] = int(float(raw))
        else:
            processed[key] = float(raw)
    return processed, unprocessed


def edges_inside_cell(
    positions: np.ndarray,
    edges: Iterable[Sequence[int]],
    tol: float = 1e-6,
) -> List[Tuple[int, int]]:
    """最小像が格納座標と一致する辺だけ残す。

    一致しない辺は周期境界をまたいでいる。ノードはここでは落とさない。
    面の上（成分が 0）のサイトも呼び出し側の座標配列に残る。
    """
    pos = np.asarray(positions, dtype=float)
    kept: List[Tuple[int, int]] = []
    for i, j in edges:
        i = int(i)
        j = int(j)
        d = pos[j] - pos[i]
        d = d - np.floor(d + 0.5)
        if np.linalg.norm(pos[i] + d - pos[j]) > tol:
            continue
        kept.append((i, j))
    return kept


def dumps(genice, scale: float = 50.0, R: float = 0.05, fn: int = 20, command_line: str = "", **kwargs) -> str:
    pos = np.asarray(genice.lattice_sites, dtype=float)
    cell = np.asarray(genice.cell, dtype=float)
    edges = edges_inside_cell(pos, genice.graph.edges())
    coord = pos @ cell * float(scale)
    radius = float(R) * float(scale)
    fn = int(fn)
    lines: List[str] = []
    if command_line:
        for line in str(command_line).splitlines():
            lines.append(f"// {line}")
    lines.append(bondfunc.strip("\n"))
    lines.append(torusfunc.strip("\n"))
    body = pipework(coord, edges, radius, fn=fn).strip("\n")
    if body:
        lines.append(body)
    return "\n".join(lines) + "\n"


def dump(genice, file=sys.stdout, **options):
    file.write(dumps(genice, **options))
