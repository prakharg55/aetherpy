#!/usr/bin/env python
# Copyright 2021, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Unit tests for plot-to-movie utilities."""

import logging
import numpy as np
import os
import pytest
import sys
import tempfile

from aetherpy.plot import movie_routines as mr


class TestMovie(object):
    """Unit tests for plot-to-movie functions."""

    def setup(self):
        """Initialize clean test environment."""

        # Create a temporary directory
        tkwargs = {}
        if sys.version_info.major >= 3 and sys.version_info.minor >= 10:
            tkwargs = {"ignore_cleanup_errors": True}
        self.tempdir = tempfile.TemporaryDirectory(**tkwargs)
        self.movie_dir = os.path.join(self.tempdir.name, "movie_dir")
        self.fileext = ".png"
        self.filebase = "test_"
        self.tempfiles = []

        return

    def teardown(self):
        """Clean up the test environment."""
        # Remove the created files and directories
        for filename in self.tempfiles:
            if os.path.isfile(filename):
                os.remove(filename)

        if os.path.isdir(self.movie_dir):
            os.rmdir(self.movie_dir)

        # Remove the temporary directory
        self.tempdir.cleanup()

        # Clear the test environment attributes
        del self.movie_dir, self.tempdir, self.tempfiles
        del self.fileext, self.filebase
        return

    def make_files(self):
        """Create a directory with temporary files."""
        os.makedirs(self.movie_dir)

        for i in range(4):
            out = tempfile.mkstemp(suffix=self.fileext, prefix=self.filebase,
                                   dir=self.movie_dir)
            self.tempfiles.append(out[1])

        return

    def test_setup_movie_dir_newdir(self):
        """Test sucessful creation of a new directory for movie files."""
        assert not os.path.isdir(self.movie_dir)
        mr.setup_movie_dir(self.movie_dir)
        assert os.path.isdir(self.movie_dir)

        return

    def test_setup_movie_dir_olddir(self):
        """Test sucessful creation of a new directory for movie files."""
        self.make_files()
        file_glob = "*".join([self.filebase, self.fileext])

        img_names = mr.setup_movie_dir(self.movie_dir, file_glob=file_glob)
        assert os.path.isdir(self.movie_dir)
        for filename in self.tempfiles:
            assert not os.path.isfile(filename), "old file not removed"

        assert img_names.find(self.filebase) >= 0, "unexpected file prefix"
        assert img_names.find(self.fileext) >= 0, "unexpected file extension"
        return

    def test_setup_movie_dir_no_overwrite(self):
        """Test raises IOError when conflicting files are present."""
        self.make_files()
        file_glob = "*".join([self.filebase, self.fileext])

        with pytest.raises(IOError) as ierr:
            mr.setup_movie_dir(self.movie_dir, file_glob=file_glob,
                               overwrite=False)

        assert str(ierr).find("files present in movie directory") >= 0
        return
