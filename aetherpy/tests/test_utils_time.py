#!/usr/bin/env python
# Copyright 2021, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Unit tests for time conversion utilities."""

import datetime as dt
import pytest

from aetherpy.utils import time_conversion as tc


class TestEpochTime(object):
    """Unit tests for epoch time conversions."""

    @pytest.mark.parametrize("esec, dtime", [
        (-10, dt.datetime(1964, 12, 31, 23, 59, 50)),
        (0, dt.datetime(1965, 1, 1)),
        (1000.5, dt.datetime(1965, 1, 1, 0, 16, 40, 500000))])
    def test_epoch_to_datetime(self, esec, dtime):
        """Test conversion from epoch seconds to datetime.

        Parameters
        ----------
        esec : int or float
            Input epoch seconds
        dtime : dt.datetime
            Output datetime object

        """

        out_time = tc.epoch_to_datetime(esec)

        assert isinstance(out_time, dt.datetime)
        assert out_time == dtime
        return

    @pytest.mark.parametrize("dtime, esec", [
        (dt.datetime(1950, 1, 1), -473385600.0),
        (dt.datetime(1965, 1, 1), 0.0),
        (dt.datetime(2010, 2, 3, 4, 5, 6, 7), 1422936306.000007)])
    def test_datetime_to_epoch(self, dtime, esec):
        """Test conversion from datetime to epoch seconds.

        Parameters
        ----------
        dtime : dt.datetime
            Input datetime object
        esec : float
            Output epoch seconds

        """

        out_time = tc.datetime_to_epoch(dtime)

        assert isinstance(out_time, float)
        assert out_time == esec
        return
