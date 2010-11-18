#!/usr/bin/env python

from sys import argv, exit
import csv

# Default values
num_bootstrap_steps = 1000
threshold = 3 # * sd

if len(argv) < 2 or len(argv) > 4:
  print 'Usage:', argv[0], '<input_file> [num_bootstrap_steps] [threshold]'
  exit()
elif len(argv) == 3:
  num_bootstrap_steps = int(argv[2])
elif len(argv) == 4:
  num_bootstrap_steps = int(argv[2])
  threshold = float(argv[3])

# Read input data and store it in data_matrix
inputfile = open(argv[1], 'r')
reader = csv.reader(inputfile, delimiter=' ')

data_matrix = []
num_pajaros = 0

first_row = True
for row in reader:
  if num_pajaros != len(row) and not first_row:
    print 'Warning: number of birds are not the same for each year'

  num_pajaros = len(row)

  if first_row:
    data_matrix = [[] for pajaro in range(num_pajaros)]
    first_row = False

  for pajaro_index in range(num_pajaros):
    if row[pajaro_index] == 'NA': # handle special case
      data_matrix[pajaro_index].append( None )
    else:
      data_matrix[pajaro_index].append( float(row[pajaro_index]) )

inputfile.close()
print 'Data matrix read'

# First filter
def remove_first_bird_without_enough_observations(matrix, min_num_obs):
  for bird_index in range(len(matrix)):
    num_obs = 0
    for obs in matrix[bird_index]:
      if obs != None:
        num_obs += 1

    if num_obs < min_num_obs:
      del matrix[bird_index]
      return True

  return False

min_num_obs = 10
while remove_first_bird_without_enough_observations(data_matrix, min_num_obs):
  pass

print 'All birds with at least', min_num_obs, 'will be analized'

# Process data matrix
import conan

g = conan.inference.maximum_entropy_with_bootstrapping( \
      data_matrix, \
      num_bootstrap_steps, \
      threshold, \
      coexpr_measure = "covariance_pairwise_complete_obs", \
      min_num_complete_obs = 5, \
      normalize = False, \
      shuffle_by_column = True \
      )

# Write output
#input_base_filename = argv[1].split('.')[0]
#output_base_filename = input_base_filename + '_ns' + str(num_bootstrap_steps) + '_threshold' + str(threshold)

#g.write_dotfile(output_base_filename + '.dot')
