# genice3-plumber

A GenIce3 exporter plugin that writes an ice hydrogen-bond network as an [OpenSCAD](https://openscad.org) pipe model (`.scad`).

Each oxygen site is a joint. Hydrogen bonds that lie inside the cell are straight pipes. The joint depends on how many of those bonds remain:

- degree 1: end plug
- degree 2: bellows elbow
- degree 3 or more: tee (sphere, arms, and collars)

Edges whose minimum-image segment leaves the cell are omitted. Sites that sit on a cell face are kept when they still have a bond, and are not copied onto the opposite face. A site left with no remaining bonds is omitted.

## Installation

```bash
pip install genice3-plumber
```

Development install:

```bash
pip install -e .
```

## Usage

### Command line

```bash
genice3 1h -e plumber > ice.scad
genice3 1h --rep 2 2 2 -e "plumber :scale 40 :R 0.04" > ice.scad
```

### Python API

```python
from genice3.plugin import UnitCell, Exporter
from genice3.genice import GenIce3

genice = GenIce3(replication_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
genice.unitcell = UnitCell("1h")

scad = Exporter("plumber").dumps(genice)
with open("output.scad", "w") as f:
    Exporter("plumber").dump(genice, file=f)
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `scale` | 50 | Scale factor applied to coordinates and to `R` (nm to OpenSCAD units). |
| `R` | 0.05 | Pipe radius (nm). Joints are larger than this: tee sphere and arms are `1.5 R`, collars are `1.8 R`. |
| `fn` | 20 | OpenSCAD `$fn` (facet count) for pipes, elbows, and tees. Plugs keep their own facet counts. |

## Requirements

- Python >= 3.11, < 3.14
- genice3 >= 3.0.a0
- numpy >= 2.0

## License

MIT
