import argparse
import re
from collections import defaultdict
import json
import itertools
from pathlib import Path
import glob
import ast





def flatten(l):
    return [item for sublist in l for item in sublist]		


def counts_per_haplotype_cont_contigs(hap_records, bubbles_file):
	haplotype_records = open(hap_records, 'r')
	#writer4 = open(outputfile, 'w')
	with open(bubbles_file) as f:
		data = f.read()
	nonrep_bubbles = ast.literal_eval(data)
	sample = None 
	haplotype = None
	count_record= {}
	broken_contigs = 0
	unbroken_contigs = 0
	for hap in haplotype_records:
		haplotype_info = hap.split('\t')
		if 'chrom' not in haplotype_info[0] : #skip the header
			nodes_read = defaultdict(list)
			informative_nodes = {'F': {}, 'R':{}}
			insides = []
			path = haplotype_info[6].split(' ')
			if 'N' not in path:
				unbroken_contigs +=1
				#print ('path has a length ', len(path))
				#counter = 0
				sample = haplotype_info[3]
				haplotype = haplotype_info[4]

				for node in path:
					#counter+=1
					if node != 'N':
						node_id = re.split('\+|-', node.strip())[0]
						direction = ['F' if i == '+' else 'R' for i in node if i in ['+', '-']][0]
						nodes_read[node_id].append(direction)
				#for n_id in nodes_read:
				#	print (n_id, '\t', 'in_hap_chunks')
				for t in nonrep_bubbles.keys():

					assert(len(t)==2)
					
					if t[0] in nodes_read or t[1] in nodes_read:
						#print(sample, ' ', haplotype, 'new bubble', t)
						insides = flatten([nonrep_bubbles[t]])
						#if '1755790' in insides : 
						for n in insides:
							if n in nodes_read:
								#print(n, '\t', 'nonrep_bubble_inside')
								#print(n, ' is inside node and read in direction', nodes_read[n])
								for d in nodes_read[n]:
									for i in insides:
										if i not in informative_nodes[d]: #we only want to add the new bubble nodes that we see
											#print('initializing ', i, ' in ', d)
											informative_nodes[d][i] = 0
									#print('updating count of ',n ,' in ',d,  ' to ',  nodes_read[n].count(d))
									informative_nodes[d][n]= nodes_read[n].count(d)
									#print (n, ' ', d,  ' ', informative_nodes[d][n])
							
								
				#to make sure no node is present in one but also shown as absent in the other dict 
				#because of some bubble nodes in different direction as compared to the rest
				#f_keys = list(informative_nodes['F'].keys())
				#r_keys = list(informative_nodes['R'].keys())
				
				#for f in f_keys:
				#	if informative_nodes['F'][f] == 0 and f in informative_nodes['R'] and informative_nodes['R'][f] > 0:
				#		del informative_nodes['F'][f]
						
				#for f in r_keys:
				#	if informative_nodes['R'][f] == 0 and f in informative_nodes['F'] and informative_nodes['F'][f] > 0:
				#		del informative_nodes['R'][f]
				
					
				assert((sample,haplotype) not in count_record)
				count_record[(sample,haplotype)] = (informative_nodes['F'], informative_nodes['R'])
				#print(sample,'\t',haplotype, '\t', 'F', '\t', informative_nodes['F'])    
				#print(sample,'\t',haplotype, '\t', 'R', '\t',informative_nodes['R'])  
			else:
				broken_contigs +=1
			   	
	
	return count_record, broken_contigs, unbroken_contigs	
					
def build_tables(count_record, output_file, broken_contigs, unbroken_contigs):
	tables =  {}
	writer = open(output_file, 'w')
	writer.write('node' + '\t' + 'FA' + '\t' + 'FP' + '\t' + 'RA' + '\t' + 'RP' + '\t' + 'FA_samples' + '\t' + 'FP_samples'+ '\t' + 'RA_samples'+ '\t' + 'RP_samples'+ '\t'+ 'broken_contigs'+'\t'+ 'unbroken_contigs'+'\n')
	samples_dict = {}
	for line in count_record:
		sample = line[0]
		haplotype = line[1]
		F = count_record[line][0]
		R = count_record[line][1]
		nodes_to_tables = set(list(F.keys()) + list(R.keys()))
		
		#print ('node',' ','FA' ,' ' , 'FP' , ' ' , 'RA' , ' ', 'RP' )
		for node in nodes_to_tables:
			node_table = [None, None, None, None] #FA,FP,RA,RP
			f_count = F[node] if node in F else None  #store the count in forward
			r_count = R[node] if node in R else None #store the count in reverse
			if node not in samples_dict:
				samples_dict[node] = {'FA':[],'FP':[],'RA':[],'RP':[]}
			
			if f_count is not None and r_count is not None: #if a node exists in both directions in the same haplotype it's not informative for us, neitehr would be its accompannying inside nodes that would appear as absent in both forward and reverse
				pass
				
			else:
				if f_count is None:
				
					node_table[0] = 0 
					node_table[1] = 0
					
				else:
					if f_count > 0:
						node_table[0] = 0
						node_table[1] = 1
						samples_dict[node]['FP'].append(sample + '_' +haplotype)
						
					else:
						node_table[0] = 1
						node_table[1] = 0
						samples_dict[node]['FA'].append(sample + '_' +haplotype)
			
			
				if r_count is None:
				
					node_table[2] = 0
					node_table[3] = 0
				else:
				
					if r_count > 0:
						node_table[2] = 0
						node_table[3] = 1
						samples_dict[node]['RP'].append(sample + '_' +haplotype)
					else:
						node_table[2] = 1
						node_table[3] = 0
						samples_dict[node]['RA'].append(sample + '_' +haplotype)

			
				#print (sample, '\t', haplotype, '\t', node, '\t', node_table[0], ' ', node_table[1], ' ', node_table[2], ' ', node_table[3])
				if node in tables:
					tables[node] = [x + y for x, y in zip(tables[node], node_table)]
				else:
					tables[node] = node_table
		
	for node in tables:
		writer.write(str(node) + '\t' + str(tables[node][0]) + '\t' + str(tables[node][1]) + '\t' + str(tables[node][2]) + '\t' + str(tables[node][3])+ '\t'+ ','.join(samples_dict[node]['FA']) + '\t'+ ','.join(samples_dict[node]['FP']) +'\t'+ ','.join(samples_dict[node]['RA']) + '\t'+ ','.join(samples_dict[node]['RP'])+ '\t'+ str(broken_contigs) + '\t' + str(unbroken_contigs) +'\n')				
					
if __name__== "__main__":
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-bubbles', metavar='BUBBLES', required=True, help='A json file with bubbles')
	parser.add_argument('-haps', metavar='HAPLOTYPES', required=True, help='The haplotype paths found using grch38 anchors')
	parser.add_argument('-out', metavar='OUTPUT', required=True, help='Output file to write the Forward Reverse counts for the nodes')
	args = parser.parse_args()
	counts_record, broken_contigs, unbroken_contigs = counts_per_haplotype_cont_contigs(args.haps, args.bubbles)
	build_tables(counts_record, args.out, broken_contigs, unbroken_contigs )
