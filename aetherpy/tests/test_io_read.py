#!/usr/bin/env python
# Copyright 2021, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Unit tests for I/O reading utilities."""

import logging
import os
import pytest
import sys
import tempfile

from aetherpy.io import read_routines


class TestIORead(object):
    """Unit tests for file reading routines."""

    def setup(self):
        """Initialize clean test environment."""
        return

    def teardown(self):
        """Clean up the test environment."""
        return

    @pytest.mark.parametrize("in_line, out_int",
                             [("1 2 3 4", 1),
                              ("42 NEU farm SpRiNg 45.0", 42)])
    @pytest.mark.parametrize("parse", [True, False])
    def test_parse_line_into_int_and_string(self, in_line, out_int, parse):
        """Test successful int/string line parsing.

        Parameters
        ----------
        in_line : str
            Input line with correct formatting
        out_int : int
            Expected integer output
        parse : bool
            `parse_string` kwarg input

        """

        # Parse the line
        lnum, lstr = read_routines.parse_line_into_int_and_string(in_line,
                                                                  parse)

        # Evaluate the retrieved values
        assert lnum == out_int, "unexpected integer value for first column"

        if parse:
            assert in_line.find(lstr) > 0, "unexpected line values"
        else:
            assert in_line == lstr, "unexpected line values"
        return

    @pytest.mark.parametrize("in_line",
                             ["1.0 2 3 4", "NEU farm SpRiNg 45.0"])
    def test_parse_line_into_int_and_string_bad_format(self, in_line):
        """Test successful int/string line parsing.

        Parameters
        ----------
        in_line : str
            Input line with incorrect formatting
        """

        # Parse the line and check error output
        with pytest.raises(ValueError) as verr:
            read_routines.parse_line_into_int_and_string(in_line)

        # Evaluate the retrieved values
        assert str(verr).find("invalid literal for int") >= 0
        return
