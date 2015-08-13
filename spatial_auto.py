# -*- coding: utf-8 -*-
"""
    Spatial AutoCorrelation
    ~~~~~~~~~

    :copyright: (c) 2015 by Joe Hand, Santa Fe Institute.
    :license: MIT
"""
import logging
import multiprocessing as mp
import ntpath
import os
import pickle
import re
import time
import unicodedata

from datetime import timedelta
from functools import partial

from collections import OrderedDict

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
    moran.get_results('densitypop')
    """

    def __init__(self, filename, name=None):
        super(Morans, self).__init__()
        self.filename = filename
        self.shapefile = filename + '.shp'
        self.dbf = filename + '.dbf'

        if name:
            self.name = name
        else:
            self.name = os.path.splitext(ntpath.basename(self.filename))[0]

        self.results = {}

        # Calculate the faster properties on init
        self._threshold = pysal.min_threshold_dist_from_shapefile(
            self.shapefile)
        self._points_array = pysal.weights.util.get_points_array_from_shapefile(
            self.shapefile)
        self._data = pysal.open(self.dbf)
        self._columns = self._data.by_col

    @property
    def threshold(self):
        return self._threshold

    @property
    def points_array(self):
        return self._points_array

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

    def calculate_weights(self, threshold=None, p=2, *args, **kwargs):
        """
        Parameters
        ----------
        threshold  : float
                     distance band
        p          : float
                     Minkowski p-norm distance metric parameter:
                     1<=p<=infinity
                     2: Euclidean distance
                     1: Manhattan distance
        """
        if threshold is None:
            if hasattr(self, 'threshold'):
                threshold = self.threshold
            else:
                raise ValueError("Must set threshold first")
        logging.warning('{}: Treshold = {}'.format(self.name, threshold))
        logging.info('{}: Starting weight calculation'.format(self.name))
        t = time.process_time()

        self.weights = pysal.DistanceBand(
            self.points_array, threshold=threshold, p=p, *args, **kwargs)

        logging.debug('{}: Weight calculation elapsed time {}'.format(
            self.name, str(timedelta(seconds=time.process_time() - t))))
        return self.weights

    def calculate_morans(self, columns, overwrite=False, *args, **kwargs):
        if not hasattr(self, 'weights'):
            # TODO: add id variable here idVariable='ID'
            self.calculate_weights(threshold=self.threshold)
        for col in columns:
            if not overwrite and col in self.results:
                continue
            y = np.array(self.data.by_col(col))
            # TODO: is float always what we want? (morans breaks w/ string)
            y = y.astype(float)
            mi = pysal.Moran(y, self.weights, *args, **kwargs)
            self.results[col] = mi
        logging.info('{}: Finished Moran Calculation'.format(self.name))
        return self.results

    def get_results(self, column, print_results=True):
        """ Quick way to nicely print results with Pandas Series
        """
        mi = self.results[column]
        results = OrderedDict([
            ('COLUMN', column),
            ("Moran's Index", mi.I),
            ('Expected Index', mi.EI),
            ('Variance', mi.VI_norm),
            ('z-score', mi.z_norm),
            ('p-value', mi.p_norm),
            ('threshold', self.threshold)
        ])
        if print_results:
            results_string = '\n'
            for key, val in results.items():
                if isinstance(val, float):
                    results_string += '\t {}: {:>10.7f}'.format(key, val)
                else:
                    results_string += '\t {}: {}'.format(key, val)
                results_string += '\n'
            return results_string
        return results

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
        """ Return path for out_dir
            Creates directory if it doesn't exist
        """
        path = os.path.dirname(self.shapefile)
        out_dir = os.path.join(path, out_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return out_dir

    def _get_filename(self):
        """ Return filename for source shapefile
        """
        return os.path.splitext(ntpath.basename(self.shapefile))[0]

    def _slugify(self, value):
        """
        From Django source.
        Converts to lowercase, removes non-word characters (alphanumerics and
        underscores) and converts spaces to hyphens. Also strips leading and
        trailing whitespace.
        """
        value = unicodedata.normalize('NFKD', value).encode(
            'ascii', 'ignore').decode('ascii')
        value = re.sub('[^\w\s-]', '', value).strip().lower()
        return re.sub('[-\s]+', '-', value)

    def _create_filtered_shapefile(self, value):
        """ Return new shapefile path/name.shp
            Creates a shapefile from source, based on filtered value
        """
        input_layer = self.input_ds.GetLayer()
        query_str = '"{}" = "{}"'.format(self.field, value)
        # Filter by our query
        input_layer.SetAttributeFilter(query_str)
        driver = ogr.GetDriverByName('ESRI Shapefile')
        out_shapefile = self._value_to_fname_path(value)
        # Remove output shapefile if it already exists
        if os.path.exists(out_shapefile):
            driver.DeleteDataSource(out_shapefile)
        out_ds = driver.CreateDataSource(out_shapefile)
        out_layer = out_ds.CopyLayer(input_layer, str(value))
        del input_layer, out_layer, out_ds
        return out_shapefile

    def _get_unique_values(self):
        """ Return unique values of filter from source shapefile.
        """
        sql = 'SELECT DISTINCT "{}" FROM {}'.format(
            self.field, self.filename)
        layer = self.input_ds.ExecuteSQL(sql)
        values = []
        for feature in layer:
            values.append(feature.GetField(0))
        return values

    def _value_to_fname_path(self, value):
        """ Return full filename path for shapefile from query value
        """
        value = value.split('-')[0]  # Hack to make US City names prettier
        value = self._slugify(value)
        fname = "{}.shp".format(value)
        return os.path.join(self.out_dir, fname)

    def _shapefile_exists(self, value):
        """ Return boolean
            Does shapefile exist (uses query value, not fname).
        """
        return os.path.isfile(self._value_to_fname_path(value))

    def create_all_shapefiles(self, overwrite=False):
        """ Returns list of new shapefiles
            Creates shapefiles for filtered data from source shapefile.
        """
        shapefiles = []
        values = self._get_unique_values()
        logging.info('Creating {} Shapefiles'.format(len(values)))
        for val in values:
            # TODO: make this multiprocess also, too slow for big filters
            if overwrite or not self._shapefile_exists(val):
                out_file = self._create_filtered_shapefile(val)
                logging.debug('Shapefile created: {}'.format(val))
                shapefiles.append(out_file)
            else:
                logging.debug('Shapefile exists, skipped: {}'.format(val))
                shapefiles.append(self._value_to_fname_path(val))
        return shapefiles


def run_single_morans(file, analysis_columns):
    named_path = os.path.splitext(file)[0]
    filename = os.path.splitext(os.path.basename(file))[0]
    logging.info('{}: Starting Analysis'.format(filename.upper()))
    moran = Morans(named_path, name=filename.upper())
    moran_results = moran.calculate_morans(analysis_columns)

    results = {}
    for col in analysis_columns:
        results[col] = moran.get_results(col, print_results=False)
    return (filename, results)


def run_moran_analysis(source_shapefile, analysis_columns,
                       filter_column=None, mp=True):
    """
    1. Filter Shapefile by filter_column
    2. Run Moran's analysis for each shapefile, each analysis column
    3. Return all results
    """
    if filter_column:
        logging.info('Running Shape Filter using: {}'.format(filter_column))
        shapefilter = ShapeFilter(source_shapefile, filter_column)
        files = shapefilter.create_all_shapefiles()
        logging.info('Created {} new shapefiles: {}'.format(len(files), files))
    else:
        logging.info('No Shapefilter, Analyzing Source Shapefile')
        files = [source_shapefile]
    if mp:
        results = _moran_mp(files, analysis_columns)
    else:
        results = []
        for i, file in enumerate(files):
            filename = os.path.splitext(os.path.basename(file))[0]
            results.append(run_single_morans(file, analysis_columns))
            logging.debug('{} of {} done'.format(i, len(files)))
    return results


class Worker(mp.Process):

    def __init__(self, task_queue, done_q, counter, total):
        super(Worker, self).__init__()
        self.task_queue = task_queue
        self.done_q = done_q
        self.counter = counter
        self.total = total

    def run(self):
        while True:
            task = self.task_queue.get()
            if task is None:
                self.task_queue.task_done()
                break
            results = run_single_morans(*task)
            self.task_queue.task_done()
            self.done_q.put(results)
            self.counter.value += 1
            logging.debug(
                '\n\t{} of {} Done\n'.format(self.counter.value, self.total))
        return


def _moran_mp(files, cols):
    """ Runs Morans code with multiprocessing module
        processes number could be tuned to computer

        Returns ALL the results at the end.
    """
    num_threads = 20

    tasks = mp.JoinableQueue()
    results_queue = mp.Queue()
    count = mp.Value('i', 0)

    # Start consumers
    workers = [Worker(tasks, results_queue, count, len(files))
               for i in range(num_threads)]
    for w in workers:
        w.start()

    # Enqueue jobs
    for i, f in enumerate(files):
        tasks.put((f, cols))
    # Add a stop marker for each worker
    for i in range(num_threads):
        tasks.put(None)
    # Wait for all of the tasks to finish
    tasks.join()

    results = [results_queue.get() for f in files]
    return results
