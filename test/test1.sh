#!/bin/sh
set -e
out=$(python3 -m genice3.cli.genice 1h -e plumber)
echo "$out" | grep -q 'module bond'
echo "$out" | grep -q 'bond('
python3 << 'PY'
import io
from contextlib import redirect_stdout

buf = io.StringIO()
with redirect_stdout(buf):
    import genice3_plumber.pipework as pipework
assert buf.getvalue() == ""
assert "echo(" not in pipework.bondfunc

import numpy as np
from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter
from genice3_plumber.pipework import END_FOOT, END_PLUG, plug
from genice3_plumber.exporter.plumber import dumps, edges_inside_cell, parse_options

proc, rest = parse_options({"scale": ["40"], "R": ["0.04"], "fn": ["12"], "nope": 1})
assert proc == {"scale": 40.0, "R": 0.04, "fn": 12} and rest == {"nope": 1}
proc, rest = parse_options({})
assert proc == {} and rest == {}

pos = np.array([
    [0.0, 0.25, 0.25],
    [0.25, 0.25, 0.25],
    [0.75, 0.25, 0.25],
])
kept = edges_inside_cell(pos, [(0, 1), (0, 2), (1, 2)])
assert kept == [(0, 1)], kept

class Graph:
    def edges(self):
        return [(0, 1), (0, 2), (1, 2)]

class Ice:
    lattice_sites = pos
    cell = np.eye(3)
    graph = Graph()

scad = dumps(Ice(), scale=10, R=0.05, fn=12, command_line="genice3 test")
assert scad.startswith("// genice3 test\n")
assert "module bond" in scad and "bond(" in scad
assert "bond([0.0,2.5,2.5]" in scad
assert "7.5,2.5,2.5" not in scad
assert "bond([0.0,2.5,2.5],[7.5,2.5,2.5]" not in scad

coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
lines, progress = plug(coord, 0, 1, 0.1, nodetype={1: END_PLUG})
assert lines and abs(float(progress[0]) - 0.1) < 1e-12
lines, progress = plug(coord, 0, 1, 0.1, nodetype={1: END_FOOT})
assert lines and abs(float(progress[0]) - 0.15) < 1e-12

g = GenIce3(replication_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
g.unitcell = UnitCell("1h")
all_edges = list(g.graph.edges())
kept = edges_inside_cell(g.lattice_sites, all_edges)
assert 0 < len(kept) < len(all_edges)
ice = Exporter("plumber").dumps(g, command_line="genice3 1h -e plumber")
assert "module bond" in ice and "bond(" in ice
assert "module partial_bellow_torus" in ice
assert "np.float" not in ice
print("ok", len(kept), "of", len(all_edges), "edges;", len(ice), "chars")
PY
