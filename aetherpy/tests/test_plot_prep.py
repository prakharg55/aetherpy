#!/usr/bin/env python
# Copyright 2021, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Unit tests for input utilities."""

import logging
import numpy as np
import pytest

from aetherpy.plot import data_prep


class TestSlice(object):
    """Unit tests for slicing functions."""

    def setup(self):
        """Initialize clean test environment."""

        self.in_coords = {"lon": np.arange(-180, 180, 1.0),
                          "lat": np.arange(-89, 90, 0.5),
                          "alt": np.arange(100, 10000, 100.0)}

        return

    def teardown(self):
        """Clean up the test environment."""
        del self.in_coords

    @pytest.mark.parametrize("cut_val, isgrid, cut_coord", [
        (0, True, "lon"), (0, True, "lat"), (0, True, "alt"),
        (50.0, False, "lon"), (50.0, False, "lat"), (600, False, "alt")])
    def test_get_cut_index(self, cut_val, isgrid, cut_coord):
        """Test sucessful construction of slicing indices.

        Parameters
        ----------
        cut_val : int or float
            Data value or grid number along which icut will be set
        isgrid : bool
            Flag that indicates `cut_val` is a grid index if True or that it is
            a data value if False
        cut_coord : str
            Expects one of 'lat', 'lon', or 'alt' and will return an index for
            that data, allowing a 2D slice to be created along other two
            coordinates

        """

        # Get the desired slice
        out = data_prep.get_cut_index(
            self.in_coords["lon"], self.in_coords["lat"],
            self.in_coords["alt"], cut_val, isgrid, cut_coord)

        # Test the output
        assert len(out) == 5, "unexpected number of returned variables"
        assert len(out[1]) == 3, "unexpected number of coordinates"

        if isgrid:
            assert out[0] == cut_val, "unexpected grid slice intersection"
        else:
            assert self.in_coords[cut_coord][out[0]] == cut_val, \
                "unexpected value slice intersection"

        iout_coords = 2
        for i, coord in enumerate(["lon", "lat", "alt"]):
            if isinstance(out[1][i], slice):
                assert np.all(out[iout_coords]
                              == self.in_coords[coord][out[1][i]]), \
                    "unexpected {:s} slice".format(coord)
                iout_coords += 1
            else:
                assert out[-1] == self.in_coords[coord][out[1][i]], \
                    "unexpected {:s} value".format(coord)
        
        return

    @pytest.mark.parametrize("cut_val", [-1, 400])
    @pytest.mark.parametrize("cut_coord", ["lon", "lat", "alt"])
    def test_bad_index_get_cut_index(self, cut_val, cut_coord):
        """Test raises ValueError when requested index is out of range.

        Parameters
        ----------
        cut_val : int
           Grid number along which icut will be set
        cut_coord : str
            Expects one of 'lat', 'lon', or 'alt' and will return an index for
            that data, allowing a 2D slice to be created along other two
            coordinates

        """

        with pytest.raises(ValueError) as verr:
            data_prep.get_cut_index(
                self.in_coords["lon"], self.in_coords["lat"],
                self.in_coords["alt"], cut_val, True, cut_coord)

        assert str(verr).find("Requested cut is outside the index range") >= 0
        return

    @pytest.mark.parametrize("cut_val", [-400, 400000])
    @pytest.mark.parametrize("cut_coord", ["lon", "lat", "alt"])
    def test_bad_value_get_cut_index(self, cut_val, cut_coord):
        """Test raises ValueError when requested index is out of range.

        Parameters
        ----------
        cut_val : float
            Data value along which icut will be set
        cut_coord : str
            Expects one of 'lat', 'lon', or 'alt' and will return an index for
            that data, allowing a 2D slice to be created along other two
            coordinates

        """

        with pytest.raises(ValueError) as verr:
            data_prep.get_cut_index(
                self.in_coords["lon"], self.in_coords["lat"],
                self.in_coords["alt"], cut_val, False, cut_coord)

        assert str(verr).find("Requested cut is outside the coordiante") >= 0
        return

    @pytest.mark.parametrize("lowlim", [True, False])
    @pytest.mark.parametrize("isgrid", [True, False])
    @pytest.mark.parametrize("cut_coord", ["lon", "lat", "alt"])
    def test_suspect_get_cut_index(self, caplog, lowlim, isgrid, cut_coord):
        """Test raises log warning when requested cut is not recommended.

        Parameters
        ----------
        lowlim : bool
            Data will be cut at the lower not-recommended limit if True, or
            along the higher not-recommended limit if False
        isgrid : bool
            Flag that indicates `cut_val` is a grid index if True or that it is
            a data value if False
        cut_coord : str
            Expects one of 'lat', 'lon', or 'alt' and will return an index for
            that data, allowing a 2D slice to be created along other two
            coordinates

        """

        # Get the desired cut value
        if lowlim:
            if cut_coord == "alt":
                # There are no problems at the lower altitude limit
                return
            else:
                if isgrid:
                    cut_val = 0
                else:
                    cut_val = self.in_coords[cut_coord][0]
        else:
            if isgrid:
                cut_val = len(self.in_coords[cut_coord]) - 1
            else:
                cut_val = self.in_coords[cut_coord][-1]

        # Raise the expected warning
        with caplog.at_level(logging.WARNING, logger="aetherpy"):
            data_prep.get_cut_index(
                self.in_coords["lon"], self.in_coords["lat"],
                self.in_coords["alt"], cut_val, isgrid, cut_coord)

        # Test the logger warning message
        captured = caplog.text
        if cut_coord == "alt":
             assert captured.find("Requested altitude slice is above the") >= 0
        else:
            assert captured.find("beyond the recommended limits") >= 0
        return
