from __future__ import annotations

def p1_gradient(points, values):
    (x0,y0),(x1,y1),(x2,y2) = ((float(x),float(y)) for x,y in points)
    f0,f1,f2 = (float(v) for v in values)
    twice_area = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
    if twice_area <= 0.0: raise ValueError("triangle must be strictly CCW")
    return ((f0*(y1-y2)+f1*(y2-y0)+f2*(y0-y1))/twice_area, (f0*(x2-x1)+f1*(x0-x2)+f2*(x1-x0))/twice_area)