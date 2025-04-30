#!/usr/bin/env python3
import geopandas
import numpy
import pandas
from libpysal.graph import Graph
from shapely import (
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import linemerge, polygonize, split

from .gap import gaps
from .hole import missing_interiors
from .overlap import is_overlapping, overlaps

__all__ = [
    "non_planar_edges",
    "planar_enforce",
    "is_planar_enforced",
    "fix_npe_edges",
    "insert_intersections",
    "self_intersecting_rings",
    "check_validity",
]


def non_planar_edges(gdf):
    """Find coincident nonplanar edges

    Parameters
    ----------

    gdf :  GeoDataFrame with polygon (multipolygon) GeoSeries


    Returns
    -------

    missing : libpysal.graph.Graph
        graph encoding nonplanar relationships between polygons

    Examples
    --------
    >>> c1 = [[0,0], [0, 10], [10, 10], [10, 0], [0, 0]]
    >>> p1 = Polygon(c1)
    >>> c2 = [[10, 2], [10, 8], [20, 8], [20, 2], [10, 2]]
    >>> p2 = Polygon(c2)
    >>> gdf = geopandas.GeoDataFrame(geometry=[p1, p2])
    >>> geoplanar.non_planar_edges(gdf).adjacency
    focal  neighbor
    0      0           0
    1      1           0
    Name: weight, dtype: int64

    """
    vertex_queen = Graph.build_contiguity(gdf, rook=False, strict=False)
    strict_queen = Graph.build_fuzzy_contiguity(gdf)
    return strict_queen.difference(vertex_queen)


def planar_enforce(gdf):
    uu = gdf.unary_union
    geoms = [uu.intersection(row.geometry) for index, row in gdf.iterrows()]
    return geopandas.GeoDataFrame(geometry=geoms)


def is_planar_enforced(gdf, allow_gaps=False):
    """Test if a geodataframe has any planar enforcement violations

    Parameters
    ----------
    gdf: GeoDataFrame with polygon geoseries for geometry
    allow_gaps: boolean
        If True, allow gaps in the polygonal coverage

    Returns
    -------
    boolean
    """

    if is_overlapping(gdf):
        return False
    if non_planar_edges(gdf):
        return False
    if not allow_gaps:
        _gaps = gaps(gdf)
        if _gaps.shape[0] > 0:
            return False
    return True


def fix_npe_edges(gdf, inplace=False):
    """Fix all npe intersecting edges in geoseries.

    Arguments
    ---------

    gdf: GeoDataFrame with polygon geoseries for geometry


    Returns
    -------
    gdf: GeoDataFrame with geometries respected planar edges.

    Examples
    --------
    >>> c1 = [[0,0], [0, 10], [10, 10], [10, 0], [0, 0]]
    >>> p1 = Polygon(c1)
    >>> c2 = [[10, 2], [10, 8], [20, 8], [20, 2], [10, 2]]
    >>> p2 = Polygon(c2)
    >>> gdf = geopandas.GeoDataFrame(geometry=[p1, p2])
    >>> geoplanar.non_planar_edges(gdf)
    defaultdict(set, {0: {1}})

    >>> gdf1 = geoplanar.fix_npe_edges(gdf)
    >>> geoplanar.non_planar_edges(gdf1)
    defaultdict(set, {})


    """
    if not inplace:
        gdf = gdf.copy()

    edges = non_planar_edges(gdf)
    unique_edges = edges.adjacency[
        ~pandas.DataFrame(numpy.sort(edges.adjacency.index.to_frame(), axis=1))
        .duplicated()
        .values
    ].index
    for i, j in unique_edges:
        poly_a = gdf.geometry.loc[i]
        poly_b = gdf.geometry.loc[j]
        new_a, new_b = insert_intersections(poly_a, poly_b)
        poly_a = new_a
        gdf.loc[i, gdf.geometry.name] = new_a
        gdf.loc[j, gdf.geometry.name] = new_b
    return gdf


def insert_intersections(poly_a, poly_b):
    """Correct two npe intersecting polygons by inserting intersection points
    on intersecting edges
    """
    overlapping_msg = (
        "Polygons are overlapping. Fix overlaps before fixing nonplanar edges."
    )
    pint = poly_a.intersection(poly_b)
    if isinstance(pint, LineString | MultiLineString | GeometryCollection):
        if isinstance(pint, GeometryCollection):
            for geom in pint.geoms:
                if isinstance(geom, Polygon | MultiPolygon):
                    raise ValueError(overlapping_msg)
        return poly_a.union(pint), poly_b.union(pint)
    elif isinstance(pint, Point | MultiPoint):
        new_polys = []
        for poly in [poly_a, poly_b]:
            if isinstance(poly, MultiPolygon):
                new_parts = []
                for part in poly.geoms:
                    exterior = LineString(list(part.exterior.coords))
                    splits = split(exterior, pint).geoms
                    if len(splits) > 1:
                        left, right = splits
                        exterior = linemerge([left, right])
                        part = Polygon(exterior)
                    new_parts.append(part)
                new_poly = MultiPolygon(new_parts)
                new_polys.append(new_poly)
            else:
                exterior = LineString(list(poly.exterior.coords))
                splits = split(exterior, pint).geoms
                if len(splits) > 1:
                    left, right = splits
                    exterior = linemerge([left, right])
                    new_poly = Polygon(exterior)
                    new_polys.append(Polygon(new_poly))
                else:
                    new_polys.append(poly)
        return new_polys
    else:  # intersection is Polygon
        raise ValueError(overlapping_msg)


def self_intersecting_rings(gdf):
    sirs = []
    for i, geom in enumerate(gdf.geometry):
        if not geom.is_valid:
            sirs.append(i)
    return sirs


def fix_self_intersecting_ring(polygon):
    p0 = polygon.exterior
    mls = p0.intersection(p0)
    polys = polygonize(mls)
    return MultiPolygon(polys)


def check_validity(gdf):
    gdfv = gdf.copy()
    sirs = self_intersecting_rings(gdf)
    if sirs:
        for i in sirs:
            fixed_i = fix_self_intersecting_ring(gdfv.geometry.iloc[i])
            gdfv.geometry.iloc[i] = fixed_i

    _gaps = gaps(gdfv)
    _overlaps = overlaps(gdfv)
    violations = {}
    violations["selfintersectingrings"] = sirs
    violations["gaps"] = _gaps
    violations["overlaps"] = _overlaps
    violations["nonplanaredges"] = non_planar_edges(gdfv)
    violations["missinginteriors"] = missing_interiors(gdfv)
    return violations
