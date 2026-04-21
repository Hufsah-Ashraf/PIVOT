###6.1..2025...in this version I changed the criteria for anchor node selection.Previously it only added a new anchor node to the list if it was either on the left side of the leftmost already found or to the right side of the rightmost already found (which I think is unnecessary..everything inside the inversion+flank existing in both orientations should classify as anchor)
import argparse
import re
from collections import defaultdict
from collections import Counter
import numpy as np 




def parse_gfa(filename):

	"""
	Read nodes and their sequences from GFA
	and store them.
	"""
	node_lengths = {}
	node_unique_chars = {}
	with open(filename, 'r') as gfa_file:
		for line in gfa_file:
			if not line[0] in ['S']:
				# we are only interested in the segments
				continue
			fields = line.split()
			node_lengths[fields[1]] = len(fields[2]) #store the length of a node
			node_unique_chars[fields[1]] = len(set(fields[2])) # for later use during safe node finding

	return (node_lengths, node_unique_chars)
	

def find_safe_nodes_1(filename, node_lengths, node_unique_chars, safe_len_limit, excluded_haps):

	"""
	Read paths from the GFA and find nodes greater than a set threshold that exist only once across all assemblies
	
	"""
	
	with open(filename, 'r') as gfa_file:
		print ('storing paths')
		paths = defaultdict(list)
		for line_r in gfa_file:
			line = line_r.split('\t')
			if  line[0] in ['P']:
				sample_info = line[1].split('#')
				sample = sample_info[0]
				haplotype = sample_info[1]
				if f"{sample}#{haplotype}" not in excluded_haps:
					path = [x[:-1].strip() for x in line[2].split(',') if (node_lengths[x[:-1].strip()] >= safe_len_limit and node_unique_chars[x[:-1].strip()] >= 3)] #atleast 3 unique characters to exclude homopolymers
					paths[(sample,haplotype)] += path
		paths_unique = defaultdict(list)
		print('done storing paths')
		return paths

def find_safe_nodes_2(paths, node_lengths, safe_length):

	"""
	Read paths from the GFA and find nodes greater than a set threshold that exist only once across all assemblies
	
	"""
	print ('finding safe nodes')
	safe_nodes = []
	intersection_counter = 0
	for sample in paths:
		samp_node_list = []
		sample_dict = Counter(paths[sample])#get how many times each element appears in the path
		for it in sample_dict:
			if sample_dict[it] == 1 and node_lengths[it]>=safe_length:
				samp_node_list.append(it)
		if len(safe_nodes) == 0 and intersection_counter == 0:
	        #print ('safe_node is empty', safe_nodes)
			safe_nodes = samp_node_list
	        #print ('so adding', safe_nodes)
		elif len(safe_nodes) == 0 and intersection_counter > 0:
	        #print('since two samples showed no common node we can already finish')
			break
		elif len(safe_nodes) > 0:
	        #print ('intersecting' , set(safe_nodes), 'and' , set(samp_node_list))
			safe_nodes = list(set(safe_nodes) & set(samp_node_list))
			intersection_counter +=1
			#print ('leading to ', safe_nodes)

	#print(safe_nodes)
	safe_nodes_dict ={}
	for safe in safe_nodes:
		safe_nodes_dict[safe] = True
	return (safe_nodes_dict)					
			
	
def read_reference_find_nodes(reference, paths, filename, chromosome, inv_start, inv_end, node_lengths, safe_nodes, flank):
	"""
	Read the reference path for the inversion
	belonging chromosome and keep a track of the 
	nodes and indices to reach the inversion start
	"""
	gfa_file= open(filename, 'r')
	start = int(inv_start)-flank
	end = int(inv_end)+flank
	node_counts = defaultdict(int)
	global safe_length
	global safe_lim
	print ('safe_length is ', safe_length, ' and safe_lim is ', safe_lim)

	for line_r in gfa_file:
		line = line_r.split('\t')
		nodes_of_interest = defaultdict(list)
		nodes_traversal = defaultdict(list)
		ref_path = []
		if  line[0] in ['P'] and reference in line[1] and chromosome in line[1] :
			right_safe_node = None
			left_safe_node = None
			lock_left = False
			inside_inv = []
			path= line[2].split(',')
			ref_start = -1
			ref_end = -1
			node_num = 0 #node number in the path being traversed
			for node in path:
				
				node_num += 1 #we read the next node in the path
				direction = ['+' if i == '+' else '-' for i in node if i in ['+', '-']][0]
				node_id = re.split('\+|-', node.strip())[0]
				node_counts[node_id] += 1
				node_length = int(node_lengths[node_id])#check the length of this node
				ref_start = ref_end + 1 #move to the next node's first character
				ref_end += node_length #see where this new node ends
				nodes_traversal[node_num] = (node_id, direction, ref_start, ref_end)  # keep a record of when this node appeared in the path and in which direction)
				#now check where we are w.r.t to the point of start
				if node_id in safe_nodes and lock_left == False:
					left_safe_node = (node_id,direction)
				if (ref_end >= start and ref_start <= end): ##we are inside the inversion+flank region
					#ref_path.append(node_id)
					lock_left = True
					if (ref_end >= inv_start and ref_start <= inv_end): #we are inside the inversion breakpoints
						inside_inv.append(node_id)
					
					if direction not in nodes_of_interest[node_id]:
						nodes_of_interest[node_id].append(direction) # keep unique orientations for each node traversed
				elif (ref_start > end and lock_left == True): #now we need to look for the safe node on right side
					if node_id in safe_nodes:
						right_safe_node = (node_id,direction)
						break
			if (left_safe_node == None or right_safe_node == None):
			    if(safe_length > safe_lim): #only run this if here is margin to reduce safe length
			        safe_length -= 5
			        print ('since one side had no safe node, looking for safe nodes with lower length ', safe_length)
			        safe_nodes = []
			        safe_nodes = find_safe_nodes_2(paths, node_lengths, safe_length)
			        read_reference_find_nodes(reference, paths,filename, chromosome, inv_start, inv_end, node_lengths, safe_nodes, flank)
			   
			break	
		else:
			pass
	print (left_safe_node, right_safe_node)
		
	return (find_anchors_from_nodes(nodes_of_interest, nodes_traversal), left_safe_node, right_safe_node, inside_inv)
		
			
def find_anchors_from_nodes(nodes_of_interest, nodes_traversal):
	left, right = None, None
	#potential_anchors = [] 
	#this keeps a record of conflicting nodes in case the inverted flanks don't look the way they should i.e. not as 1+,2+,3+......3-,2-,1-.
	#So it would e.g. append a node in case it's first occurrence is to the left but it's last occurrence is before the anchor already found
	#e.g. 1+,2+, ......1-,2- so in this case 1 and 2 would be conflicting anchors
	#final_anchor = tuple() #this keeps only one anchor which can be used as such in case there are no conflicts
	anchor_path = {}
	for r in nodes_of_interest:
		if len(nodes_of_interest[r]) > 1: #these would be the nodes that occurred in both forward and reverse orientation
			anchor_path[r]= True
			# #occurrences= [k for k, v in nodes_traversal.items() if (v == (r,'+') or v == (r,'-'))]
			# #leftmost_occurrence = min([k for k, v in nodes_traversal.items() if v[0] == r]) #where was the node first found in the path
			# #leftmost_occurrence_direction = [v[1] for k, v in nodes_traversal.items() if k == leftmost_occurrence][0] #what was the direction of the node when it was first found
			# ##print('leftmost_occ_direction is', leftmost_occurrence_direction)
			# dir_to_find = ['-' if leftmost_occurrence_direction == '+' else '+'][0]
			# #print('direction to find is', dir_to_find)
			# rightmost_occurrence = max([k for k, v in nodes_traversal.items() if (v[0],v[1]) == (r,dir_to_find)]) #where was the node last found in a reverse orientation
			# #leftmost_occurrence_refstart= [v[2] for k, v in nodes_traversal.items() if k == leftmost_occurrence][0] # refstart and refend for the leftmost ocurrence of the respective node r
			# #rightmost_occurrence_refend= [v[3] for k, v in nodes_traversal.items() if k == rightmost_occurrence][0]

			# if (left is None and right is None) or (leftmost_occurrence < left and rightmost_occurrence > right): #it's either the first node found or its both directions are outside the ones already found
			# 	left = int(leftmost_occurrence)
			# 	right = int(rightmost_occurrence)
			# 	#final_anchor = (r, leftmost_occurrence, rightmost_occurrence, leftmost_occurrence_refstart, rightmost_occurrence_refend)
			# 	anchor_path[r] = True
			# elif (leftmost_occurrence > left  and rightmost_occurrence < right): # both directions of this node exist inside the structure we already have found (a different version of the following if statement to allow for variants in the flanks)
			# 	anchor_path[r] = True	
				
			#elif (leftmost_occurrence == left + len(anchor_path) and rightmost_occurrence == right - len(anchor_path)): # both directions of this node exist inside the structure we already have found
			#	anchor_path[r] = True
				
			#elif (leftmost_occurrence == left + len(anchor_path) and rightmost_occurrence == right + len(anchor_path)):#for cases like chr8 where flanks don't have the conventional structure
			#	anchor_path[r] = True
				
			#elif leftmost_occurrence < left or rightmost_occurrence > right:
				#potential_anchors.append((r, leftmost_occurrence, rightmost_occurrence, leftmost_occurrence_refstart, rightmost_occurrence_refend))
			#	left = int(leftmost_occurrence)
			#	right = int(rightmost_occurrence)
			# else:
			# 	pass
	#print(anchor_path)
	return(anchor_path)
			


def find_anchors_in_haplotypes(gfa_file, anchor_path, writer3, chrom, inv_start, inv_end, left_safe_node, right_safe_node, flank, excluded_haps):
	"""
	Read the paths for all other haplotypes and try to locate the region using anchor nodes
	"""
	#after the table building is decided and in place. we want to go back to the question of how to find the haplotypes
	#this implementation doesn't do the alignment because apparently it is not needed at least for the moment
	#the above implementations also have some issue because they weren't picking up some haplotype chunks even though they existed e.g. the chm13 one
	gfa= open(gfa_file, 'r')
	haplotypes_record = defaultdict(list)
	seen_samples = defaultdict(int)
	paths = defaultdict(list)
	#sample_cons={('HG02145','2')}
	for line_r in gfa:
		line = line_r.split('\t')
		if  line[0] in ['P'] :
			nodes_of_interest = defaultdict(list)
			nodes_traversal = defaultdict(list)
			sample_info = line[1].split('#')
			sample = sample_info[0]
			haplotype = sample_info[1]
			if f"{sample}#{haplotype}" not in excluded_haps:
				keep_direction_left = None
				keep_direction_right = None
				if (seen_samples[(sample,haplotype)] == 1):
					#print(sample, ' ', haplotype, ' occurred more than once but we already had found a chunk ')
					pass
				elif (seen_samples[(sample,haplotype)] == 0 or seen_samples[(sample,haplotype)] is None):
					if( seen_samples[(sample,haplotype)] == 0):
						#print(sample, ' ', haplotype, ' occurred more than once and we still havent found a chunk ')
						pass
					else:
					
						#print(sample, ' ', haplotype, ' seen for the first time ')
						pass
					print(sample, ' ', haplotype)
					path = line[2].split(',')
					paths[(sample,haplotype)].append(path)
					structure = []
					node_num = 0 #node number in the path being traversed
					for node in path:
						node_num += 1 #we read the next node in the path
						direction = ['+' if i == '+' else '-' for i in node if i in ['+', '-']][0]
						node_id = re.split('\+|-', node.strip())[0]
						#now check where we are w.r.t to the point of start
						if node_id in anchor_path:
							#print (node, ' ', node_num)
							nodes_traversal[node_num] = (node_id, direction)  # keep a record of when this node appeared in the path and in which direction) 
							if direction not in nodes_of_interest[node_id]:

								nodes_of_interest[node_id].append(direction) # keep unique orientations for each node traversed
						elif left_safe_node is not None and node_id == left_safe_node[0]:
							if (direction == left_safe_node[1]):
								keep_direction_left = True

							else :
								keep_direction_left = False

							print (node_id, ' is the left_safe with direction ', direction)
						elif right_safe_node is not None and node_id == right_safe_node[0]:
							if direction == right_safe_node[1] :
								keep_direction_right = True
							else:
								keep_direction_right = False
							print (node_id, ' is the right_safe with direction ', direction)
								
					left, right = None, None
					for r in nodes_of_interest:
						if len(nodes_of_interest[r]) > 1:
							#occurrences= [k for k, v in nodes_traversal.items() if (v == (r,'+') or v == (r,'-'))]
							leftmost_occurrence = min([k for k, v in nodes_traversal.items() if v[0] == r]) #where was the node first found in the path
							leftmost_occurrence_direction = [v[1] for k, v in nodes_traversal.items() if k == leftmost_occurrence][0] #what 														was the direction of the node when it was first found
							#print('leftmost_occ_direction is', leftmost_occurrence_direction)
							dir_to_find = ['-' if leftmost_occurrence_direction == '+' else '+'][0]
							#print('direction to find is', dir_to_find)
							rightmost_occurrence = max([k for k, v in nodes_traversal.items() if (v[0],v[1]) == (r,dir_to_find)]) #where was 												the node last found in a reverse orientation

							if (left is None and right is None) or (leftmost_occurrence < left and rightmost_occurrence > right):
								left = int(leftmost_occurrence)
								right = int(rightmost_occurrence)
								
							elif leftmost_occurrence < left and rightmost_occurrence < right:
								#print (sample, ' ', haplotype, 'left occurrence to the left but rightmost not to the right')
								#potential_anchors.append((r, leftmost_occurrence, rightmost_occurrence, leftmost_occurrence_refstart, 												rightmost_occurrence_refend))
								left = int(leftmost_occurrence)
								
							elif leftmost_occurrence > left and rightmost_occurrence > right:
								#print (sample, ' ', haplotype, 'right occurrence to the right but leftmost not to the left')
								right = int(rightmost_occurrence)
							else:
								pass
					if left is not None and right is not None:

						print ('structure found by just using anchor nodes')
						####we still need to figure out in which direction to read the structure
						structure = path[left-1: right] #left and right come from node numbers in the path that is 1 based
						if ((keep_direction_left == True and keep_direction_right == True) or (keep_direction_left == True and keep_direction_right == None) or (keep_direction_right == True and keep_direction_left == None)): #STRICTER CONDITION + a bit flexible condition for cases where we don't see a safe node on one side
							print('and would be kept in same orientation')
							
						elif ((keep_direction_left == False and keep_direction_right == False) or (keep_direction_left == False and keep_direction_right == None) or (keep_direction_right == False and keep_direction_left == None)): 
							print('and would be flipped')
							inverted_path = []
							for s in structure:

								direction = ['+' if i == '-' else '-' for i in s if i in ['+', '-']][0]
								inverted_path.append(s[:-1]+direction) #same node inverted orientation

							structure = inverted_path
							#structure.reverse()
							#chunk = ' '.join(structure)
							#writer3.write(str(chrom) + '\t' + str(inv_start) + '\t' + str(inv_end) +  '\t'  + str(sample) + '\t' + str(haplotype) +'\t' + str(flank) + '\t' + chunk +  '\n' )
							#seen_samples[(sample,haplotype)] = 1
							
						else:
							#raise Exception("Sorry, the safe nodes are not helping")
							print(sample, haplotype, 'but safe nodes didnt help so considering it with broken contig cases')
							seen_samples[(sample,haplotype)] = 0
							continue
							
						chunk = ' '.join(structure)
						writer3.write(str(chrom) + '\t' + str(inv_start) + '\t' + str(inv_end) +  '\t'  + str(sample) + '\t' + str(haplotype) +'\t' + str(flank) + '\t' + chunk +  '\n' )
						seen_samples[(sample,haplotype)] = 1
						
						
					else:
						print ('structure not found in this contig')
						seen_samples[(sample,haplotype)] = 0
					
	gfa.close()
	for samples in seen_samples:
		
		if seen_samples[samples] == 0 and right_safe_node is not None and left_safe_node is not None:
			print(samples)
			print('structure needs to be found from broken contigs usings safe nodes')
			structure = find_chunks_using_safe_nodes(paths[samples], left_safe_node, right_safe_node, anchor_path)
			if structure is not None:
				writer3.write(str(chrom) + '\t' + str(inv_start) + '\t' + str(inv_end) +  '\t'  + str(samples[0]) + '\t' + str(samples[1]) + '\t'
				+ str(flank) + '\t' + structure +  '\n' )
			else:
				writer3.write(str(chrom) + '\t' + str(inv_start) + '\t' + str(inv_end) +  '\t'  + str(samples[0]) + '\t' + str(samples[1]) + '\t'
				+ str(flank) + '\t' + 'more than one left or right structures' +  '\n' )				
			
def find_chunks_using_safe_nodes(paths, left_safe_node, right_safe_node, anchor_path):
	left_chunks = []
	right_chunks = []
	path_num=0
	print(len(paths), 'for this haplotyoe')
	for path in paths:
		path_num+=1
		print (path_num)
		left_safe = None
		right_safe = None
		anchor_found = None
		keep_direction_left = None
		keep_direction_right = None
		chunk = []
		node_num = 0
		inverted_path = []
		for node in path:
			node_num += 1
			direction = ['+' if i == '+' else '-' for i in node if i in ['+', '-']][0]
			node_id = re.split('\+|-', node.strip())[0]
			inv_direction = ['-' if i == '+' else '+' for i in node if i in ['+', '-']][0]
			inv_node = node_id + inv_direction #same node but with opposite direction to be kept for alignment later
			inverted_path.append(inv_node) # to be used in case we have to reverse the contig
			if node_id in anchor_path and anchor_found is None:
				print (node_id, ' is the anchor')
				anchor_found = node_num - 1
				
			elif node_id == left_safe_node[0]:
				if (direction == left_safe_node[1]):
					keep_direction_left = True

				else :
					keep_direction_left = False

				left_safe = node_num - 1
				print (node_id, ' is the left_safe with direction ', direction)
			elif node_id == right_safe_node[0]:
				if direction == right_safe_node[1] :
					keep_direction_right = True
				else:
					keep_direction_right = False
				right_safe = node_num - 1
				print (node_id, ' is the right_safe with direction ', direction)
		#assert((left_safe is None and right_safe is not None) or (right_safe is None and left_safe is not None) and anchor_found is not None) #I think both safe nodes shouldn't be in the same contig	

		if left_safe is not None:
			if keep_direction_left == True: #means the contig is in the correct orientation so we need to store the chunk from left safe node to the end of the contig as it is
				chunk = path[left_safe:]
			elif keep_direction_left == False: #means the contig is in reverse orientation so we need to store the chunk from the start of the contig to the left safe node and reverse its' orientation
				chunk = inverted_path[0 : left_safe+1]
				chunk.reverse()
			left_chunks.append(chunk)
				
		if right_safe is not None:
			if keep_direction_right == True: #means the contig is in the correct orientation so we need to store the chunk from left safe node to the end of the contig as it is
				chunk = path[0:right_safe+1]
			elif keep_direction_right ==False: #means the contig is in reverse orientation so we need to store the chunk from the start of the contig to the left safe node and reverse its' orientation
				chunk = inverted_path[right_safe:]
				chunk.reverse()
			right_chunks.append(chunk)
	
	if len(left_chunks) == 1 and len(right_chunks) == 1:#which ideally should be the case
		print('left_chunk starts with ', left_chunks[0][0], ' and ends at ', left_chunks[0][-1])
		print('right_chunk starts with ', right_chunks[0][0], ' and ends at ', right_chunks[0][-1])
		structure = left_chunks[0] + ['N'] + right_chunks[0]
		return (' '.join(structure))
	else:
		print ('we found ', len(left_chunks), ' left chunks and ', len(right_chunks), ' right chunks')
		return None	
		
		
		

	
if __name__== "__main__":
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-invs', metavar='INVS', required=True, help='A bed file with inversion co-ordinates')
	parser.add_argument('-gfa', metavar='GFA', required=True, help='graph in GFA format.')
	parser.add_argument('-anchor', metavar='FINAL', required=True, help='One anchor path per region')
	parser.add_argument('-hapchunks', metavar='HAP', required=True, help='Output file to write the haplotype paths flanked by anchor nodes')
	parser.add_argument('-flank', metavar='FLANK', type=int, required=True, help='The length of the flanking region to consider')
	parser.add_argument('-limit', metavar='FLANK_L', type=int, required=True, help='The atmost length of the flanking region to consider')
	parser.add_argument('-safe_len', metavar='SAFE', type=int, required=True, help='The length of the safe node to consider')
	parser.add_argument('-safe_len_limit', metavar='SAFE_LIM', type=int, required=True, help='The smallest length allowed for finding safe nodes')
	parser.add_argument('-ref', metavar='REFERENCE', type=str, required=True, help='The haplotype to consider as reference')
	parser.add_argument('-exhaps', metavar='EXCLUDED', type=str, required=True, help='The haplotypes to be excluded for the analysis')


	args = parser.parse_args()
	
	
	#parse the GFA
	#print('Reading sequence information from GFA file...')
	excluded_haps = {x.strip() for x in args.exhaps.split(',')}
	node_lengths, node_unique_chars = parse_gfa(args.gfa)
	paths= find_safe_nodes_1(args.gfa, node_lengths, node_unique_chars, args.safe_len_limit, excluded_haps)
	safe_nodes = []
	safe_length = args.safe_len
	safe_lim = args.safe_len_limit
	while len(safe_nodes)==0 and safe_length >=safe_lim: #to make sure the safe nodes we find are atleast a certain length
		safe_nodes = find_safe_nodes_2(paths, node_lengths, safe_length)
		print ('looking for safe nodes with length atleast ', safe_length)
		safe_length -=5
	safe_length +=5 #because this is the value we last checked
	#print ('looking for safe nodes with length atleast ', args.safe_len)
	#print('Done reading GFA.')
	#start looking for each inversion
	inversions= open(args.invs, 'r')
	writer1 = open(args.anchor, 'w')
	writer1.write('chrom' + '\t' + 'inv_start' + '\t' + 'inv_end' +  '\t' +'anchor_path' +  '\t' +'inside_inv_nodes'+  '\n' )
	for invs in inversions:
		inv = invs.strip().split('\t')
		#print('Working on ', invs)
		chrom, inv_start, inv_end = inv[0], int(inv[1]), int(inv[2])
		anchor_path = {}
		i=0
		while (len(anchor_path)==0 and int(args.flank)+i <= int(args.limit)):
		    print ('looking for anchors with flank ', int(args.flank)+i, ' and safe length ', safe_length)
		    anchor_path, left_safe_node, right_safe_node, inside_inv =  read_reference_find_nodes(args.ref, paths, args.gfa, chrom, inv_start, inv_end, node_lengths, safe_nodes, int(args.flank)+i) #because we update the safe length before after last use but here we need to give the one that was actually used)
		    i+=10000
			
		writer1.write(str(chrom) + '\t' + str(inv_start) + '\t' + str(inv_end) +  '\t' + ' '. join(anchor_path.keys()) + '\t'+ ' '. join(inside_inv) +  '\n' )
		#as we only deal with one inversion in one iteration of this script, we can delete the stored dictionaries to save space
		del node_lengths
		del node_unique_chars
		del safe_nodes
		writer3 = open(args.hapchunks,'w')
		writer3.write('chrom' + '\t' + 'inv_start' + '\t' + 'inv_end' + '\t' +'sample' + '\t' + 'haplotype' + '\t' + 'flank_value'+ '\t' +'hap_chunk' +  '\n' )	
		if (len(anchor_path) != 0):
			find_anchors_in_haplotypes(args.gfa, anchor_path, writer3, chrom, inv_start, inv_end, left_safe_node, right_safe_node, int(args.flank)+i, excluded_haps)
	inversions.close()
	writer1.close()
	writer3.close()
			
			
			
			
