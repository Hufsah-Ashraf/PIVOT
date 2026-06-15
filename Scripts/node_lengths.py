import pyximport
import argparse
import re
from collections import defaultdict
from collections import Counter
import numpy as np




def parse_gfa(filename):

	"""
	Read nodes and their sequences from GFA
	and outputs their length.
	"""

	with open(filename, 'r') as gfa_file:
		for line in gfa_file:
			if not line[0] in ['S']:
				# we are only interested in the segments
				continue
			fields = line.split()
			print(fields[1], '\t', len(fields[2]))
		

	
if __name__== "__main__":
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-gfa', metavar='GFA', required=True, help='graph in GFA format.')
	args = parser.parse_args()
	parse_gfa(args.gfa)
	
