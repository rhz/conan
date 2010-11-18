#!/usr/bin/env python

from sys import argv
from math import sqrt, log
import csv

from matplotlib import rcParams

golden_mean = (sqrt(5) - 1.0) / 2.0
rcParams.update( \
  {'backend': 'Agg',
   'axes.labelsize': 18,
   'text.fontsize': 10,
   'xtick.labelsize': 18,
   'ytick.labelsize': 18,
   'font.family': 'serif',
   'figure.figsize': [ 12.0, golden_mean * 12.0 ],
   'text.usetex': False })

#from numpy import arange, array, concatenate, sqrt, add
from scipy import stats
# pylab = matplotlib
from pylab import *

column_number = {}
column_number['T'], column_number['num_vertices'], column_number['num_edges'], \
    column_number['asp'], column_number['entropy'], column_number['dd_entropy'] = range(6)
# dd_entropy = degree distribution entropy

for filename in argv[1:]:
  data = {}
  reader = csv.reader(open(filename, 'rb'), delimiter='\t') # read inputfile as csv

  # collect data from csv reader
  for row in reader:
    T = row[column_number['T']]

    if T not in data:
      data[T] = {}
      data[T]['asp'] = {}
      data[T]['entropy'] = {}
      data[T]['dd_entropy'] = {}

    num_vertices = int(row[column_number['num_vertices']])

    data[T]['asp'][num_vertices] = float(row[column_number['asp']])
    data[T]['entropy'][num_vertices] = float(row[column_number['entropy']])
    data[T]['dd_entropy'][num_vertices] = float(row[column_number['dd_entropy']])

  for T in sorted(data.keys()):
    lines = []
    names = []
    for prop in data[T].keys():
      x_array = []
      y_array = []

      line_type = ''
      if prop == 'asp':
        line_type = 'r-'
      elif prop == 'entropy':
        line_type = 'b-'
      elif prop == 'dd_entropy':
        line_type = 'g-'

      for num_vertices in sorted(data[T][prop].keys()):
        x_array.append(log(num_vertices))
        y_array.append(log(data[T][prop][num_vertices]))

      lines.append(plot(x_array, y_array, line_type))
      names.append(prop)

    legend(lines, names)
    savefig('log' + T.replace('-', '_') + '.png', dpi=150)
    clf()
