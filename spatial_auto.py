# -*- coding: utf-8 -*-
"""
    Spatial AutoCorrelation
    ~~~~~~~~~

    :copyright: (c) 2015 by Joe Hand, Santa Fe Institute.
    :license: MIT
"""
import logging
import ntpath
import os
import pickle

import numpy as np
import pandas as pd
import pysal

from osgeo import ogr


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
        self._threshold = pysal.min_threshold_dist_from_shapefile(
            self.shapefile)
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
                threshold = self.threshold
            else:
                raise ValueError("Must set threshold first")
        logging.warning('Treshold = {}'.format(threshold))
        logging.info('Starting weight calculation (slow)')
        self.weights = pysal.threshold_binaryW_from_shapefile(
            self.shapefile, threshold, *args, **kwargs)
        return self.weights

    def calculate_morans(self, columns, overwrite=False, *args, **kwargs):
        if not hasattr(self, 'weights'):
            self.calculate_weights(threshold=self.threshold)
        for col in columns:
            if not overwrite and col in self.results:
                continue
            y = np.array(self.data.by_col(col))
            mi = pysal.Moran(y, self.weights, *args, **kwargs)
            self.results[col] = mi
        return self.results

    def print_results(self, column):
        """ Quick way to nicely print results with Pandas Series
        """
        mi = self.results[column]
        results = {
            'COLUMN': column,
            "Moran's Index": mi.I,
            'Expected Index': mi.EI,
            'Variance': mi.VI_norm,
            'z-score': mi.z_norm,
            'p-value': mi.p_norm,
            'threshold': self.threshold
        }
        pd.options.display.float_format = '{:10.7f}'.format
        return pd.Series(results)

    def pickle_results(self, column):
        # TODO
        pass


class ShapeFilter(object):

    """ ShapeFilter

    Filter single shapefile by a field,
    creating new shapefiles for each value.

    """

    def __init__(self, shapefile, filter_field, out_dir='tmp'):
        super(ShapeFilter, self).__init__()
        self.shapefile = shapefile
        self.field = filter_field

        self.input_ds = ogr.Open('{}'.format(shapefile))
        self.filename = self._get_filename()
        self.out_dir = self._get_create_out_dir(out_dir)

    def _get_create_out_dir(self, out_dir):
        path = os.path.dirname(self.shapefile)
        out_dir = os.path.join(path, out_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return out_dir

    def _get_filename(self):
        return os.path.splitext(ntpath.basename(self.shapefile))[0]

    def _create_filtered_shapefile(self, value):
        input_layer = self.input_ds.GetLayer()
        query_str = "'{}' = {}".format(self.field, value)
        # Filter by our query
        input_layer.SetAttributeFilter(query_str)
        driver = ogr.GetDriverByName('ESRI Shapefile')
        out_shapefile = "{}/{}.shp".format(self.out_dir, value)
        # Remove output shapefile if it already exists
        if os.path.exists(out_shapefile):
            driver.DeleteDataSource(out_shapefile)
        out_ds = driver.CreateDataSource(out_shapefile)
        out_layer = out_ds.CopyLayer(input_layer, str(value))
        del input_layer, out_layer, out_ds
        return out_shapefile

    def _get_unique_values(self):
        sql = "SELECT DISTINCT '{}' FROM {}".format(
            self.field, self.filename)
        layer = self.input_ds.ExecuteSQL(sql)
        values = []
        for feature in layer:
            values.append(feature.GetField(0))
        return values

    def create_all_shapefiles(self):
        shapefiles = []
        for val in self._get_unique_values():
            out_file = self._create_filtered_shapefile(val)
            shapefiles.append(out_file)
        return shapefiles


def run_single_morans(file, analysis_columns):
    named_path = os.path.splitext(file)[0]
    moran = Morans(named_path)
    moran.calculate_morans(analysis_columns)
    return moran


def run_moran_analysis(source_shapefile, analysis_columns, filter_column=None):
    """
    1. Filter Shapefile by filter_column
    2. Run Moran's analysis for each shapefile, each analysis column
    3. Return all results
    """
    results = {}
    if filter_column:
        logging.info('Running Shape Filter using: {}'.format(filter_column))
        shapefilter = ShapeFilter(source_shapefile, filter_column)
        files = shapefilter.create_all_shapefiles()
    else:
        logging.info('No Shapefilter, Analyzing Source Shapefile')
        files = [source_shapefile]
    for file in files:
        filename = os.path.splitext(os.path.basename(file))[0]
        logging.info('Starting Analysis of {}'.format(filename))
        results[filename] = run_single_morans(file, analysis_columns)
    return results
