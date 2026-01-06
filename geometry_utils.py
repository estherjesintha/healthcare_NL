# geometry_utils.py
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

class InvalidGeometryError(Exception):
    """Custom exception for invalid geometries."""
    pass

def make_square(x, y, size):
    """Return a square polygon centered on (x, y)."""
    # STEP 1: implement in Part 2
    if size <= 0:
            raise ValueError("Size must be positive.")
    half = size / 2
    coords = [
        (x - half, y - half),
        (x + half, y - half),
        (x + half, y + half),
        (x - half, y + half),
        (x - half, y - half)
    ]
    return Polygon(coords)

def is_valid_geometry(geom):
    """Return True if geom is not None and geom.is_valid."""
    # STEP 2: implement in Part 2
    if geom is None:
        return False
    return geom.is_valid

def union_geometries(geoms):
    """Return unary union of a list of geometries.
    Raise InvalidGeometryError if any is invalid.
    """
    # STEP 3: implement in Part 2
    for g in geoms:
        if not is_valid_geometry(g):
            raise InvalidGeometryError("Invalid geometry detected.")
    return unary_union(geoms)

def buffer_point(x, y, dist):
    """Create a point and buffer it by dist (dist >= 0)."""
    # STEP 4: implement in Part 2
    if dist < 0:
        raise ValueError("Distance must be non-negative.")
    return Point(x, y).buffer(dist)