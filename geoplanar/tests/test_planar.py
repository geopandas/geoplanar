#!/usr/bin/env python3
import geopandas as gpd
import numpy as np
import pytest
from libpysal.graph import Graph
from numpy.testing import assert_equal
from shapely.geometry import MultiPolygon, Polygon

import geoplanar


@pytest.fixture
def test_polygons():
    # Polygon A: square
    coords1 = [[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]]
    # Polygon B: touches Polygon A at an edge but not properly snapped (non-planar)
    coords2 = [[10, 2], [10, 8], [20, 8], [20, 2], [10, 2]]
    poly1 = Polygon(coords1)
    poly2 = Polygon(coords2)
    gdf = gpd.GeoDataFrame(geometry=[poly1, poly2])
    return gdf


class TestPlanar:
    def setup_method(self):
        # Setup multipolygon with potential non-planar edge
        c1 = [[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]]
        p1 = Polygon(c1)

        c2 = [[10, 2], [20, 8], [20, 2], [10, 2]]
        p2a = Polygon(c2)
        p2b = Polygon([[21, 2], [21, 4], [23, 3]])
        p2 = MultiPolygon([p2a, p2b])

        self.gdf = gpd.GeoDataFrame(geometry=[p1, p2])
        self.gdf_str = self.gdf.set_index(np.array(["foo", "bar"]))

    def test_non_planar_edges(self):
        res = geoplanar.non_planar_edges(self.gdf)
        assert res.equals(Graph.from_dicts({0: [1], 1: [0]}))

        gdf1 = geoplanar.fix_npe_edges(self.gdf)
        assert_equal(
            gdf1.geometry[0].wkt,
            "POLYGON ((0 0, 0 10, 10 10, 10 2, 10 0, 0 0))"
        )

        res_str = geoplanar.non_planar_edges(self.gdf_str)
        assert res_str.equals(Graph.from_dicts({"foo": ["bar"], "bar": ["foo"]}))

        gdf1_str = geoplanar.fix_npe_edges(self.gdf_str)
        assert_equal(
            gdf1_str.geometry.iloc[0].wkt,
            "POLYGON ((0 0, 0 10, 10 10, 10 2, 10 0, 0 0))"
        )

    def test_check_validity(self, test_polygons):
        violations = geoplanar.check_validity(test_polygons)
        assert "nonplanaredges" in violations
        assert len(violations["nonplanaredges"].adjacency) > 0
        assert isinstance(violations["gaps"], gpd.GeoSeries)
        assert isinstance(violations["overlaps"], np.ndarray)
        assert violations["overlaps"].shape[1] == 0 or violations["overlaps"].shape[1] == 2
        

