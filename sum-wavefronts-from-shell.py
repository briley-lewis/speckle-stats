import statistics as stats
import sys

"""just an easy framework for calling sum_wavefronts from command line - python3 sum-wavefronts-from-shell.py filename"""

filename = sys.argv[1]
stats.sum_wavefronts(filename)