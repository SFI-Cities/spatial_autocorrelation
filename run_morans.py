# -*- coding: utf-8 -*-
"""
    Spatial AutoCorrelation - Run Script
    ~~~~~~~~~

    Run with:
        python run_morans.py shapefiles/my_shapes.shp \
            'Density_Col' 'Population_Col' 'etc_col' -f 'Filter_Col'

    :copyright: (c) 2015 by Joe Hand, Santa Fe Institute.
    :license: MIT
"""
import argparse
import logging
import pickle
import time

from datetime import timedelta

import pandas as pd

from spatial_auto import run_moran_analysis

parser = argparse.ArgumentParser(
    description="Run Moran's on Shapefile (optional filter)")
parser.add_argument('shapefile', type=str,
                    help='source shapefile for analysis')
parser.add_argument('analysis_vars', nargs='+',
                    help='columns to run Morans I analysis on')
parser.add_argument('-f', '--filter', type=str,
                    help='Filter Shapefile by Column (Optional)')
parser.add_argument('--show-logs', dest='log', action='store_true')
parser.add_argument('--no-logs', dest='log', action='store_false')
parser.set_defaults(log=True)
parser.add_argument('--logs-file', dest='log_file', action='store_true')
parser.set_defaults(log_file=False)
args = parser.parse_args()

if __name__ == '__main__':
    t = time.process_time()
    if args.log_file:
        logging.basicConfig(format='%(asctime)s \n \t %(message)s',
                            filename='morans.log', level=logging.DEBUG,
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    elif args.log:
        logging.basicConfig(format='%(asctime)s \n \t %(message)s',
                            level=logging.DEBUG,
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(format='%(asctime)s \n \t %(message)s',
                            level=logging.WARNING,
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    if args.filter:
        filter_col = args.filter
    else:
        filter_col = None

    logging.info('Starting Analysis')
    try:
        results = run_moran_analysis(
            args.shapefile, args.analysis_vars, filter_column=filter_col)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        logging.exception('\n\nError: ')

    logging.info('Finished Calculations \n\n')

    try:
        results_df = []
        keys = []
        for shapefile, values in results:
            df = pd.DataFrame(values).transpose()
            if shapefile in keys:
                val = 1
                while True:
                    shapefile = '{}_{}'.format(shapefile,val)
                    if shapefile not in keys:
                        break
                    val += 1
            keys.append(shapefile)
            del(df['COLUMN']) # add this as a key later
            results_df.append(df)

            results_log = '{} RESULTS \n'.format(shapefile.upper())
            results_log += df.to_string()
            results_log += '\n'
            logging.debug(results_log)
        results_df = pd.concat(results_df, keys=keys, names=['CITY', 'COLUMN'], axis=0)
        results_df.to_csv('results.csv')
        pickle.dump(results, open( "results.p", "wb" ))
    except:
        logging.exception('Some error exporting results: ')
    logging.debug('Total elapsed time {}'.format(
        str(timedelta(seconds=time.process_time() - t))))
