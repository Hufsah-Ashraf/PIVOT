import pyximport
import argparse
import re
from collections import defaultdict
import json
import itertools
from pathlib import Path
import glob



def store_bubbles(bubbles_input):
	# Opening JSON file with bubbles
	#print('storing bubbles')
	f = open(bubbles_input)
	bubbles_data = json.load(f)
	bubble_types = {}
	bubbles = {}
	ends = []
	#store just the bubbles i.e. ends and insides
	for i in bubbles_data:
		for b in bubbles_data[i]['bubbles']:
			#print(b['ends'], ' ', b['type'])
			#print(b['inside'])
			assert tuple(b['ends']) not in bubbles #there should be one bubble per node pair
			if b['type'] != 'super':
				bubbles[tuple(b['ends'])] = b['inside']
			#if b['type'] not in bubble_types:
			#	bubble_types[b['type']]= True
	#print('done storing bubbles')
	#for bubble in bubbles:
	#	print (bubble)
	#	print (bubbles[bubble])
	return bubbles	


def flatten(l):
    return [item for sublist in l for item in sublist]		


def clean_bubbles_old(filename,bubbles,excluded_haps):
	with open(filename, 'r') as gfa_file:
		#print ('finding clean bubbles')
		paths = defaultdict(list)
		for line_r in gfa_file:
			line = line_r.split('\t')
			if  line[0] in ['P']:
				sample_info = line[1].split('#')
				sample = sample_info[0]
				haplotype = sample_info[1]
				if f"{sample}#{haplotype}" not in excluded_haps:
					path = [x for x in line[2].split(',')]
					paths[(sample,haplotype)] += path
	
	sample = None 
	haplotype = None
	#count_record= {}
	bubbles_repeated= {}
	bubbles_not_repeated = {}
	for hap in paths:
		nodes_read = defaultdict(list)
		insides = []
		path = paths[hap]
		#print ('path has a length ', len(path))
		#counter = 0
		sample = hap[0]
		haplotype = hap[1]
		#print (sample, ' ' , haplotype)
		for node in path:
			#counter+=1
			if node != 'N':
				node_id = re.split('\+|-', node.strip())[0]
				direction = ['F' if i == '+' else 'R' for i in node if i in ['+', '-']][0]
				nodes_read[node_id].append(direction)
		for t in bubbles.keys():
			assert(len(t)==2)
			if t[0] in nodes_read and t[1] in nodes_read:
				key1= len(nodes_read[t[0]])
				key2= len(nodes_read[t[1]])
				#print(t, ' ' , key1, ' ', key2)
				if key1 > 1 and key2 >1:
					bubbles_repeated[t]= True
					if t in bubbles_not_repeated:
						del bubbles_not_repeated[t]
				elif key1 ==1 and key2 ==1 and t not in bubbles_repeated:
					bubbles_not_repeated[t] = bubbles[t]
					
	#for bubble in bubbles_not_repeated:
	#	print(bubble, ' not repeated')
	#	print(bubbles[bubble])
	#for bubble in bubbles_repeated:
	#	print(bubble)
	#	print(bubbles[bubble])
	
	print(bubbles_not_repeated)




		   


					
			
					
if __name__== "__main__":
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-bubbles', metavar='BUBBLES', required=True, help='A json file with bubbles')
	#parser.add_argument('-out', metavar='OUTPUT', required=True, help='Output file to write the Forward Reverse counts for the nodes')
	parser.add_argument('-gfa', metavar='GFA', required=True, help='graph in GFA format.')
	parser.add_argument('-exhaps', metavar='EXCLUDED', type=str, required=True, help='The haplotypes to be excluded for the analysis')

	args = parser.parse_args()
	excluded_haps = {x.strip() for x in args.exhaps.split(',')}
	bubbles= store_bubbles(args.bubbles)
	clean_bubbles_old(args.gfa,bubbles, excluded_haps)
