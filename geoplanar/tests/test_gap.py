#!/usr/bin/env python3
import os.path

import geopandas
import numpy
import pytest
from geopandas.testing import assert_geodataframe_equal
from numpy.testing import assert_equal, assert_array_almost_equal
from packaging.version import Version
from shapely.geometry import Polygon, box

from geoplanar import fill_gaps, gaps, non_planar_edges, snap

HERE = os.path.abspath(os.path.dirname(__file__))
PACKAGE_DIR = os.path.dirname(os.path.dirname(HERE))
_TEST_DATA_DIR = os.path.join(PACKAGE_DIR, "geoplanar", "tests", "test_data")


class TestGap:
    def setup_method(self):
        self.p1 = box(0, 0, 10, 10)
        self.p2 = Polygon([(10, 10), (12, 8), (10, 6), (12, 4), (10, 2), (20, 5)])
        self.gdf = geopandas.GeoDataFrame(geometry=[self.p1, self.p2])
        self.gdf_crs = self.gdf.set_crs(3857)
        self.gdf_str = self.gdf_crs.set_index(numpy.array(["foo", "bar"]))

    def test_gaps(self):
        h = gaps(self.gdf_crs)
        assert_equal(h.area.values, numpy.array([4.0, 4.0]))
        assert self.gdf_crs.crs.equals(h.crs)

    def test_fill_gaps(self):
        gdf1 = fill_gaps(self.gdf)
        assert_equal(gdf1.area.values, numpy.array([108.0, 32.0]))

        gdf1 = fill_gaps(self.gdf_str)
        assert_equal(gdf1.area.values, numpy.array([108.0, 32.0]))

    def test_fill_gaps_smallest(self):
        gdf1 = fill_gaps(self.gdf, strategy="smallest")
        assert_equal(gdf1.area.values, numpy.array([100.0, 40.0]))

    def test_fill_gaps_none(self):
        gdf1 = fill_gaps(self.gdf, strategy=None)
        assert_equal(gdf1.area.values, numpy.array([108.0, 32.0]))

    def test_fill_gaps_gaps_df(self):
        gaps_df = gaps(self.gdf).loc[[0]]
        filled = fill_gaps(self.gdf, gaps_df)
        assert_equal(filled.area, numpy.array([104, 32]))

        filled = fill_gaps(self.gdf_str, gaps_df)
        assert_equal(filled.area, numpy.array([104, 32]))


@pytest.mark.skipif(
    Version(geopandas.__version__) < Version("1.0.0dev"),
    reason="requires geopandas 1.0",
)
class TestSnap:
    def setup_method(self):
        self.p1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10]])
        self.p2 = Polygon([(11, 0), (21, 0), (21, 20), (11, 20)])
        self.gdf = geopandas.GeoDataFrame(geometry=[self.p1, self.p2])
        self.gdf_str = self.gdf.set_index(numpy.array(["foo", "bar"]))

        self.p3 = Polygon([[0, 0], [10, 0], [13, 13], [3, 10]])
        self.p4 = Polygon([(10.7, 2), (23, 10), (15, 20)])
        self.p5 = Polygon([(10.7, 1.5), (23, 9), (10.3, 0)])

    def test_snap_below_threshold(self):
        gdf1 = snap(self.gdf, 0.5)
        assert_equal(gdf1.area.values, numpy.array([100.0, 200.0]))

        gdf1 = snap(self.gdf_str, 0.5)
        assert_equal(gdf1.area.values, numpy.array([100.0, 200.0]))

    def test_snap_above_threshold(self):
        gdf1 = snap(self.gdf, 1.1)
        assert_equal(gdf1.area.values, numpy.array([110.0, 200.0]))

        gdf1 = snap(self.gdf_str, 1.1)
        assert_equal(gdf1.area.values, numpy.array([110.0, 200.0]))

    def test_snap_reverse_order(self):
        gdf = geopandas.GeoDataFrame(geometry=[self.p2, self.p1])
        gdf1 = snap(gdf, 1.1)
        assert_equal(gdf1.area.values, numpy.array([210.0, 100.0]))

    def test_snap_complex_shapes(self):
        gdf = geopandas.GeoDataFrame(geometry=[self.p3, self.p4])
        gdf1 = snap(gdf, 0.5)
        assert_equal(
            numpy.round(gdf1.area.values, decimals=1), numpy.array([113.6, 93.5])
        )

    def test_snap_3shapes(self):
        gdf = geopandas.GeoDataFrame(geometry=[self.p3, self.p4, self.p5])
        gdf1 = snap(gdf, 1)
        assert_equal(
            numpy.round(gdf1.area.values, decimals=1), numpy.array([114.1, 102.3, 7.7])
        )

    def test_validity(self):
        df = geopandas.read_file(
            os.path.join(_TEST_DATA_DIR, "possibly_invalid_snap.gpkg")
        )
        snapped = snap(df, 0.5)
        assert snapped.is_valid.all()

    def test_corner_touch(self):
        p1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10]])
        p2 = Polygon([(10, 10), (20, 10), (20, 20), (10, 20)])
        gdf = geopandas.GeoDataFrame(geometry=[p2, p1])

        gdf2 = snap(gdf, threshold=8)

        # we detect corners as 'maybe snap' but nothing should happen
        assert_geodataframe_equal(gdf, gdf2)

    def test_skewed(self):
        p1 = Polygon([[0, 0], [9.9, 0], [10.1, 10], [0, 10]])
        p2 = Polygon([(11, 2), (21, 2), (21, 8), (11, 8)])

        # one direction
        gdf = geopandas.GeoDataFrame(geometry=[p1, p2])
        gdf1 = snap(gdf, 1.5)
        assert len(non_planar_edges(gdf1)) == 0
        assert_array_almost_equal(gdf1.area, [108, 60])

        # other direction
        gdf = geopandas.GeoDataFrame(geometry=[p2, p1])
        gdf1 = snap(gdf, 1.5)
        assert len(non_planar_edges(gdf1)) == 0
        assert_array_almost_equal(gdf1.area, [66.076427, 100])
