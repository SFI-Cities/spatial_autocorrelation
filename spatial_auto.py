# -*- coding: utf-8 -*-
"""
    Spatial AutoCorrelation
    ~~~~~~~~~

    :copyright: (c) 2015 by Joe Hand, Santa Fe Institute.
    :license: MIT
"""
import pickle

import numpy as np
import pandas as pd
import pysal


class Morans(object):
    """Morans

    Usage:
    from spatial_auto import Morans
    moran = Morans('shapefile_name') # Call with shapefile name, no extension
    moran.calculate_morans('column') # This can take a long time
    moran.print_results('densitypop')
    """

    def __init__(self, filename):
        super(Morans, self).__init__()
        self.filename = filename
        self.shapefile = filename + '.shp'
        self.dbf = filename + '.dbf'

        self.results = {}

        # Calculate the faster properties on init
        self._threshold = pysal.min_threshold_dist_from_shapefile(self.shapefile)
        self._data = pysal.open(self.dbf)
        self._columns = self._data.by_col

    @property
    def threshold(self):
        return self._threshold

    @property
    def data(self):
        return self._data

    @property
    def columns(self):
        return self._columns

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = value
        return self._weights

    def calculate_weights(self, threshold=None, *args, **kwargs):
        if threshold is None:
            if hasattr(self, 'threshold'):
                threshold=self.threshold
            else:
                raise ValueError("Must set threshold first")
        self.weights = pysal.threshold_binaryW_from_shapefile(
            self.shapefile, threshold, *args, **kwargs)
        return self.weights

    def calculate_morans(self, column, overwrite=False, *args, **kwargs):
        if not overwrite and column in self.results:
            return self.results[column]
        if not hasattr(self, 'weights'):
            self.calculate_weights(threshold=self.threshold)
        y = np.array(self.data.by_col(column))
        mi = pysal.Moran(y, self.weights, *args, **kwargs)
        self.results[column] = mi
        return mi

    def print_results(self, column):
        """ Quick way to nicely print results with Pandas Series
        """
        mi = self.results[column]
        results = {
            'COLUMN': column,
            "Moran's Index": mi.I,
            'Expected Index':mi.EI,
            'Variance':mi.VI_norm,
            'z-score':mi.z_norm,
            'p-value':mi.p_norm,
            'threshold' : self.threshold
        }
        pd.options.display.float_format = '{:10.7f}'.format
        return pd.Series(results)

    def write_pickle_weights(self, filename):
        if not hasattr(self, 'weights'):
            return
        with open(filename, 'wb') as f:
            pickle.dump( self.weights, f, protocol=pickle.HIGHEST_PROTOCOL)
        return

    def read_weights_pickle(self,filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)