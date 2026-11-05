# Adapted from Pipework (pipework.py). Geometry is unchanged.
# Statements are returned as text instead of printed.

from collections import defaultdict
import numpy as np
from itertools import combinations
from math import tan, acos, pi

END_PLUG = 0
END_FOOT = 1
ELBOW_PIPE = 0
ELBOW_BELLOW = 1


bondfunc = """module bond(p1=[0,0,0],p2=[1,1,1],r=1, fn=40)
{
    d = p2 - p1;
    H = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    R = sqrt(d[0]*d[0] + d[1]*d[1]);
    //atan2(y,x)
    theta = 90-atan2(d[2],R);
    phi   = atan2(d[1],d[0]);
    translate(p1)
    rotate([0,theta,phi])
    cylinder(r=r, h=H,$fn=fn);
}
"""

torusfunc = """
module torus(ri,ro, fn=40)
{
    rotate_extrude($fn=fn)
        translate([ro,0,0])
            circle(ri, $fn=fn);
}


module wedge0(r, h, a, fn=40)
{
    difference(){
        cylinder(r=r,h=h*2,center=true, $fn=fn);
        translate([-r,-r,-r])
            cube([r*2,r,r*2]);
        rotate([0,0,a])
            translate([-r,0,-r])
                cube([r*2,r,r*2]);
    }
}


module bellow_torus(ri,ro,N=100,fn=20)
{
    L = (ro - ri)*2*PI / N * 0.2;
    L2 = (ro + ri)*2*PI / N / 2;
    rotate([90,0,0])
    for(a=[1:N]){
      rotate([0,360/N*a,0])
        hull(){
            translate([ro,0,0]){
              cylinder(r=ri, h=L, $fn=fn);
            }
            rotate([0,360/N/2,0]){
                translate([ro,0,0]){
                  cylinder(r=ri-L2, h=L, $fn=fn);
                }
            }
            rotate([0,-360/N/2,0]){
                translate([ro,0,0]){
                  cylinder(r=ri-L2, h=L, $fn=fn);
                }
            }
        }
    }
}


module partial_torus(ri,ro,a,fn=40)
{
    // a: angle in 360 degree
    intersection(){
        torus(ri,ro,fn);
        wedge0(ri+ro,ri,a,fn=fn*2.5);
    }
}


module partial_bellow_torus(ri,ro,a,N=35,fn=40)
{
    // a: angle in 360 degree
    intersection(){
        bellow_torus(ri,ro,N=N,fn=fn);
        wedge0(ri+ro,ri,a,fn=fn*2.5);
    }
}
"""


def _clip_cosine(value):
    return float(np.clip(value, -1.0, 1.0))


def cylinder2(v1, v2, R, fn=40):
    a = [float(x) for x in v1]
    b = [float(x) for x in v2]
    return "bond([{0},{1},{2}],[{3},{4},{5}],{6},{7});".format(
        *a, *b, float(R), int(fn)
    )


def plug(coord, i, j, R, fn=20, nodetype={}):
    lines = []
    D = coord[j] - coord[i]
    L = np.linalg.norm(D)
    if j not in nodetype or nodetype[j] == END_PLUG:
        Lplug = R
        Rplug = R * 1.2
        Lplug2 = R / 2
        Rplug2 = R * 0.8
        lines.append(cylinder2(coord[i], coord[i] + D / L * Lplug, Rplug, fn=6))
        lines.append(cylinder2(coord[i], coord[i] - D / L * Lplug2, Rplug2, fn=4))
        return lines, [Lplug / L]
    elif nodetype[j] == END_FOOT:
        Lplug = R * 1.5
        Rplug = R * 1.3
        Lplug2 = R / 2
        Rplug2 = R * 3
        lines.append(cylinder2(coord[i], coord[i] + D / L * Lplug, Rplug, fn=fn))
        lines.append(cylinder2(coord[i], coord[i] - D / L * Lplug2, Rplug2, fn=fn))
        return lines, [Lplug / L]
    return lines, []


class sphere:
    def __init__(self, R, fn=40):
        self.s = "sphere({0}, $fn={1});".format(float(R), int(fn))

    def translate(self, V):
        xyz = [float(x) for x in V]
        self.s = "translate([{0},{1},{2}]){{ {3} }}".format(*xyz, self.s)
        return self

    def __str__(self):
        return self.s


def tee(coord, c, vertices, R, fn=20):
    lines = []
    Rarm = R * 1.4
    Rcollar = R * 1.0
    Lmax = R * 1.7
    for i, j in combinations(vertices, 2):
        Di = coord[i] - coord[c]
        Dj = coord[j] - coord[c]
        Di /= np.linalg.norm(Di)
        Dj /= np.linalg.norm(Dj)
        cosine = _clip_cosine(np.dot(Di, Dj))
        half = acos(cosine) / 2
        L = Rcollar / tan(half)
        if Lmax < L:
            Lmax = L
    Larm = Lmax + R / 8
    Lcollar = R / 2
    offset = []
    for v in vertices:
        D = coord[v] - coord[c]
        L = np.linalg.norm(D)
        lines.append(cylinder2(coord[c], coord[c] + D / L * Larm, Rarm, fn=fn))
        lines.append(
            cylinder2(
                coord[c] + D / L * (Larm + Lcollar),
                coord[c] + D / L * Larm,
                Rcollar,
                fn=fn,
            )
        )
        offset.append((Larm + Lcollar) / L)
    lines.append(sphere(Rarm, fn=fn).translate(coord[c]).s)
    return lines, offset


def elbow(coord, c, vertices, R, Relbow=0, fn=20, nodetype={}):
    lines = []
    if Relbow == 0:
        Relbow = R * 2
    Rcollar = R * 1.0
    Lcollar = R / 2
    Di = coord[vertices[0]] - coord[c]
    Dj = coord[vertices[1]] - coord[c]
    Li = np.linalg.norm(Di)
    Lj = np.linalg.norm(Dj)
    Ei = Di / Li
    Ej = Dj / Lj
    cosine = _clip_cosine(np.dot(Ei, Ej))
    angle = acos(cosine)
    half = angle / 2
    LL = Relbow / tan(half)
    Eh = Ei + Ej
    Eh /= np.linalg.norm(Eh)
    center = coord[c] + Eh * (LL**2 + Relbow**2) ** 0.5
    Ez = np.cross(Ei, Ej)
    Ez /= np.linalg.norm(Ez)
    Ex = coord[c] + LL * Ej - center
    Ex /= np.linalg.norm(Ex)
    Ey = np.cross(Ez, Ex)
    if c not in nodetype or nodetype[c] == ELBOW_BELLOW:
        torus = "partial_bellow_torus"
        Rbody = R
    elif nodetype[c] == ELBOW_PIPE:
        torus = "partial_torus"
        Rbody = R * 1.5
    else:
        torus = "partial_bellow_torus"
        Rbody = R
    lines.append(
        "multmatrix(m=[[{0},{1},{2},{3}],[{4},{5},{6},{7}],[{8},{9},{10},{11}],[{12},{13},{14},{15}]]) {20}({16},{17},{18},fn={19});".format(
            float(Ex[0]),
            float(Ey[0]),
            float(Ez[0]),
            float(center[0]),
            float(Ex[1]),
            float(Ey[1]),
            float(Ez[1]),
            float(center[1]),
            float(Ex[2]),
            float(Ey[2]),
            float(Ez[2]),
            float(center[2]),
            0,
            0,
            0,
            1,
            float(Rbody),
            float(Relbow),
            float(180 - angle * 180 / pi),
            int(fn),
            torus,
        )
    )
    offset = []
    for v in vertices:
        D = coord[v] - coord[c]
        L = np.linalg.norm(D)
        lines.append(
            cylinder2(
                coord[c] + D / L * (LL + Lcollar), coord[c] + D / L * LL, Rcollar, fn=fn
            )
        )
        offset.append((LL + Lcollar) / L)
    return lines, offset


def pipework(coord, edges, R, fn=20, nodetype={}):
    """
    coord: coordinates of the nodes
    edges: node pairs (undirected)
    R: radius of the straight pipe

    Returns OpenSCAD statements. Module definitions are not included.
    """
    lines = []
    vinfo = defaultdict(list)
    for e in edges:
        vinfo[e[0]].append(e[1])
        vinfo[e[1]].append(e[0])
    pipeends = dict()
    for e in edges:
        pipeends[e] = [0.0, 1.0]
    for v in vinfo:
        nei = vinfo[v]
        if len(nei) == 1:
            part, progress = plug(coord, v, nei[0], R, fn=fn, nodetype=nodetype)
        elif len(nei) == 2:
            part, progress = elbow(coord, v, nei, R, fn=fn)
        else:
            part, progress = tee(coord, v, nei, R, fn=fn)
        lines.extend(part)
        for i, p in enumerate(progress):
            e = (v, nei[i])
            if e in pipeends:
                pipeends[e][0] = p
            else:
                e = (nei[i], v)
                pipeends[e][1] = 1.0 - p
    for e in pipeends:
        p = pipeends[e]
        c = coord[e[0]]
        D = coord[e[1]] - c
        lines.append(cylinder2(c + D * p[0], c + D * p[1], R, fn=fn))
    return "\n".join(lines) + ("\n" if lines else "")
