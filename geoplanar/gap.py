#!/usr/bin/env python3
#
from collections import defaultdict

import geopandas
import numpy as np
import pandas as pd
import shapely
from esda.shape import isoperimetric_quotient
from packaging.version import Version

__all__ = ["gaps", "fill_gaps", "snap"]

GPD_GE_014 = Version(geopandas.__version__) >= Version("0.14.0")
GPD_GE_100 = Version(geopandas.__version__) >= Version("1.0.0dev")


def gaps(gdf):
    """Find gaps in a geodataframe.

    A gap (emply sliver polygon) is a set of points that:

    - are not contained by any of the geometries in the geoseries
    - are not contained by the external polygon

    Parameters
    ----------

    gdf :  GeoDataFrame with polygon (multipolygon) GeoSeries


    Returns
    -------

    _gaps : GeoDataFrame with gap polygons

    Examples
    --------
    >>> p1 = box(0,0,10,10)
    >>> p2 = Polygon([(10,10), (12,8), (10,6), (12,4), (10,2), (20,5)])
    >>> gdf = geopandas.GeoDataFrame(geometry=[p1,p2])
    >>> h = geoplanar.gaps(gdf)
    >>> h.area
    array([4., 4.])
    """

    polygons = geopandas.GeoSeries(
        shapely.get_parts(
            shapely.polygonize(
                [gdf.boundary.union_all() if GPD_GE_100 else gdf.boundary.unary_union]
            )
        ),
        crs=gdf.crs,
    )
    if GPD_GE_014:
        poly_idx, _ = gdf.representative_point().sindex.query(
            polygons, predicate="contains"
        )
    else:
        poly_idx, _ = gdf.representative_point().sindex.query_bulk(
            polygons, predicate="cotains"
        )

    return polygons.drop(poly_idx).reset_index(drop=True)


def fill_gaps(gdf, gap_df=None, strategy="largest", inplace=False):
    """Fill gaps in a GeoDataFrame by merging them with neighboring polygons.

    Note that as part of the gap filling heuristics, trim_overlaps is called internally.

    Parameters
    ----------
    gdf : GeoDataFrame
        A GeoDataFrame containing polygon or multipolygon geometries.

    gap_df : GeoDataFrame, optional
        A GeoDataFrame containing the gaps to be filled. If None, gaps will be
        automatically detected within `gdf`.

    strategy : {'smallest', 'largest', 'compact', None}, default 'largest'
        Strategy to determine how gaps are merged with neighboring polygons:
          - 'smallest': Merge each gap with the smallest neighboring polygon.
          - 'largest' : Merge each gap with the largest neighboring polygon.
          - 'compact' : Merge each gap with the neighboring polygon that results in
                        the new polygon having the highest compactness
                        (isoperimetric quotient).
          - None      : Merge each gap with the first available neighboring polygon
                        (order is indeterminate).

    inplace : bool, default False
        If True, modify the input GeoDataFrame in place. Otherwise, return a new
        GeoDataFrame with the gaps filled.

    Returns
    -------
    GeoDataFrame or None
        A new GeoDataFrame with gaps filled if `inplace` is False. Otherwise,
        modifies `gdf` in place and returns None.
    """
    if gap_df is None:
        gap_df = gaps(gdf)

    if not inplace:
        gdf = gdf.copy()

    if not GPD_GE_014:
        gap_idx, gdf_idx = gdf.boundary.sindex.query_bulk(
            gap_df.boundary, predicate="overlaps"
        )
    else:
        gap_idx, gdf_idx = gdf.boundary.sindex.query(
            gap_df.boundary, predicate="overlaps"
        )

    to_merge = defaultdict(set)

    for g_ix in range(len(gap_df)):
        neighbors = gdf_idx[gap_idx == g_ix]

        if strategy == "compact":
            # Find the neighbor that results in the highest IQ
            gap_geom = shapely.make_valid(gap_df.geometry.iloc[g_ix])
            best_iq = -1
            best_neighbor = None
            neighbor_geometries = gdf.geometry.iloc[neighbors].apply(shapely.make_valid)
            for neighbor, neighbor_geom in zip(neighbors, neighbor_geometries):
                combined_geom = shapely.union_all([neighbor_geom, gap_geom])
                iq = isoperimetric_quotient(combined_geom)
                if iq > best_iq:
                    best_iq = iq
                    best_neighbor = neighbor
            to_merge[best_neighbor].add(g_ix)
        elif strategy is None:  # don't care which polygon we attach cap to
            to_merge[gdf.index[neighbors[0]]].add(g_ix)
        elif strategy == "largest":
            # Attach to the largest neighbor
            to_merge[gdf.iloc[neighbors].area.idxmax()].add(g_ix)
        else:
            # Attach to the smallest neighbor
            to_merge[gdf.iloc[neighbors].area.idxmin()].add(g_ix)

    new_geom = []
    for k, v in to_merge.items():
        new_geom.append(
            shapely.union_all(
                [gdf.geometry.loc[k]] + [gap_df.geometry.iloc[i] for i in v]
            )
        )
    if isinstance(gdf, geopandas.GeoDataFrame):
        gdf.loc[list(to_merge.keys()), gdf.geometry.name] = new_geom
    else:
        gdf.loc[list(to_merge.keys())] = new_geom

    return gdf


def _get_parts(geom):
    """Get parts recursively to explode multi-part geoms in collections

    Parameters
    ----------
    geom : shapely.Geometry

    Returns
    -------
    np.array[shapely.Geometry]
    """
    parts = shapely.get_parts(geom)
    if (shapely.get_type_id(parts) > 3).any():
        return _get_parts(parts)
    return parts


def _snap(geometry, reference, threshold, segment_length):
    """Snap g1 to g2 within threshold

    Parameters
    ----------
    geometry : shapely.Polygon
        geometry to snap
    reference : shapely.Polygon
        geometry to snap to
    threshold : float
        max distance between vertices to snap
    segment_length : float
        max segment length parameter in segmentize

    Returns
    -------
    shapely.Polygon
        snapped geometry
    """
    # special case - corner touch. If the area of intersection of buffers is below
    # 110% of what a clean corner would cover, assume it is a corner that shall
    # not be snapped.
    if (
        geometry.buffer(threshold).intersection(reference.buffer(threshold)).area
        < threshold * 4 * 1.1
    ):
        return geometry

    # extract the shell and holes from the first geometry
    shell, holes = geometry.exterior, geometry.interiors
    # segmentize the shell and extract coordinates
    coords = shapely.get_coordinates(shapely.segmentize(shell, segment_length))
    # create a point geometry from the coordinates
    points = shapely.points(coords)
    # find the shortest line between the points and the second geometry
    lines = shapely.shortest_line(points, reference)
    # mask the coordinates where the distance is less than the threshold
    distance_mask = shapely.length(lines) < threshold

    # return the original geometry if no coordinates are within the threshold
    if not any(distance_mask):
        return geometry

    # update the coordinates with the snapped coordinates
    coords[distance_mask] = shapely.get_coordinates(lines)[1::2][distance_mask]
    # re-create the polygon with new coordinates and original holes and simplify
    # to remove any extra vertices.
    polygon = shapely.Polygon(coords, holes=holes)
    simplified = shapely.make_valid(polygon)
    # the function may return invalid and make_valid may return non-polygons
    # the largest polygon is the most likely the one we want
    if simplified.geom_type != "Polygon":
        parts = _get_parts(simplified)
        simplified = parts[np.argmax(shapely.area(parts))]


    # TODO: we need to add the same points to the exterior of reference, otherwise
    # floating point precision may end in overlaps.

    return simplified


def snap(geometry, threshold):
    """Snap geometries that are within threshold to each other

    Only one of the pair of geometries identified as nearby will be snapped,
    the one with the lower index.

    If the snapping heuristics leads to an invalid geometry, the function attempts
    to fix it using :func:`shapely.make_valid`, which may lead to multi-part geometries.
    If that happens, only the largest component is returned. Occasionally, this may
    lead to improper snapping.

    Parameters
    ----------
    geometry : GeoDataFrame | GeoSeries
        geometries to snap. Geometry type needs to be Polygon for all of them.
    threshold : float
        max distance between geometries to snap
        threshold should be ~10% larger than the distance between polygon edges to
        ensure snapping

    Returns
    -------
    GeoSeries
        GeoSeries with snapped geometries
    """
    # importing here to avoid circular imports
    from .planar import fix_npe_edges

    geometry = geometry.copy()
    if not GPD_GE_100:
        raise ImportError("geopandas 1.0.0 or higher is required.")

    nearby_a, nearby_b = geometry.sindex.query(
        geometry.geometry, predicate="dwithin", distance=threshold
    )
    overlap_a, overlap_b = geometry.boundary.sindex.query(
        geometry.boundary, predicate="overlaps"
    )

    self_mask = nearby_a != nearby_b
    nearby_a = nearby_a[self_mask]
    nearby_b = nearby_b[self_mask]

    self_mask = overlap_a != overlap_b
    overlap_a = overlap_a[self_mask]
    overlap_b = overlap_b[self_mask]

    nearby = pd.MultiIndex.from_arrays([nearby_a, nearby_b], names=("source", "target"))
    overlap = pd.MultiIndex.from_arrays(
        [overlap_a, overlap_b], names=("source", "target")
    )
    nearby_not_overlap = nearby.difference(overlap)
    if not nearby_not_overlap.empty:
        duplicated = pd.DataFrame(
            np.sort(np.array(nearby_not_overlap.to_list()), axis=1)
        ).duplicated()
        pairs_to_snap = nearby_not_overlap[~duplicated]

        new_geoms = []
        previous_geom = None
        snapped_geom = None
        for geom, ref in zip(
            geometry.geometry.iloc[pairs_to_snap.get_level_values("source")],
            geometry.geometry.iloc[pairs_to_snap.get_level_values("target")],
            strict=True,
        ):
            if previous_geom == geom:
                new_geoms.append(
                    _snap(
                        snapped_geom,
                        ref,
                        threshold=threshold,
                        segment_length=threshold,
                    )
                )
            else:
                snapped_geom = _snap(
                    geom, ref, threshold=threshold, segment_length=threshold
                )
                new_geoms.append(snapped_geom)
                previous_geom = geom

        if isinstance(geometry, geopandas.GeoDataFrame):
            gs = geometry.geometry
            gs.iloc[pairs_to_snap.get_level_values("source")] = new_geoms
            geometry.geometry = gs
        else:
            geometry.iloc[pairs_to_snap.get_level_values("source")] = new_geoms

    fixed_topology = fix_npe_edges(geometry)

    fixed_topology.geometry = fixed_topology.simplify_coverage(threshold / 100)
    return fixed_topology
