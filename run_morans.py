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
args = parser.parse_args()

if __name__ == '__main__':
    if args.log:
        logging.basicConfig(format='%(asctime)s \n \t %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(asctime)s \n \t %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
    if args.filter:
        filter_col = args.filter
    else:
        filter_col = None

    logging.info('Starting Analysis')
    results = run_moran_analysis(
        args.shapefile, args.analysis_vars, filter_column=filter_col)
    logging.info('Finished Calculations \n\n')
    for shapefile in results.keys():
        for var in args.analysis_vars:
            results_log = 'Results: {} - {}'.format(shapefile, var)
            results_log += results[shapefile].print_results(var)
            results_log += '\n'
            logging.info(results_log)
