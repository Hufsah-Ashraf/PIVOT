# Hufsah Ashraf (latest update 5.1.2026)
from collections import defaultdict
import glob
import re
configfile: "Snakefile.json"


w_dir = config["working_dir"]
path_to_gfas = config["path_to_graphs"]
path_to_bubblegun = config["bubble_detection"] 
path_to_invs = config["path_to_invs"]
chroms= config["chromosomes"]
names = glob_wildcards(path_to_invs + "/{n}.bed")
inv_per_chrom = defaultdict(list)
for n in names[0]:
	t = n.strip().split(':')
	inv_per_chrom[t[0]].append(t[1])

final_tables_simple = []
final_tables_all = []
final_tables_cont_contigs_simple = []
final_tables_cont_contigs_all = []
node_lengths = []


for chrom in inv_per_chrom :
    temp_list_3 = (expand(w_dir+"/final_tables/simple_bubbles/{c}:{i}.txt", c=chrom, i=inv_per_chrom[chrom]))
    final_tables_simple +=temp_list_3
    temp_list_3_1 = (expand(w_dir+"/final_tables/all_bubbles/{c}:{i}.txt", c=chrom, i=inv_per_chrom[chrom]))
    final_tables_all +=temp_list_3_1
    temp_list_8 = (expand(w_dir+"/final_tables_cont_contigs/simple_bubbles/{c}:{i}.txt", c=chrom, i=inv_per_chrom[chrom]))
    final_tables_cont_contigs_simple +=temp_list_8
    temp_list_8_1 = (expand(w_dir+"/final_tables_cont_contigs/all_bubbles/{c}:{i}.txt", c=chrom, i=inv_per_chrom[chrom]))
    final_tables_cont_contigs_all +=temp_list_8_1
    temp_list_5= expand(w_dir+"/node_lengths/{c}.txt", c=chrom, i=inv_per_chrom[chrom])
    node_lengths += temp_list_5

rule all:
    input:
    	w_dir+"/summary_table/Summary_cont_contigs_all_bubbles.tsv",
    	w_dir+"/summary_table/Summary_all_contigs_all_bubbles.tsv",
    	w_dir+"/summary_table/Summary_cont_contigs_simple_bubbles.tsv",
    	w_dir+"/summary_table/Summary_all_contigs_simple_bubbles.tsv"
        

#Step 1
def get_safe_lim(wildcards):
    if ({wildcards.c} in [{'chrX'},{'chrY'}]):
        return 5    
    else:
    	return 15
    	
def get_safe_len(wildcards):
    if ({wildcards.c} in [{'chrX'},{'chrY'}]):
        return 15    
    else:
    	return 30

rule find_anchors:
    input:
        inv = path_to_invs+"/{c}:{i}.bed",
	    gfa = path_to_gfas+"/{c}.gfa", #for example, pattern needs to be adjusted according to the input gfas being used
    output:
        node_sanity = w_dir+"/node_sanity/{c}:{i}-node_sanity.txt"
    params:
        hap_chunks = w_dir+"/hap_chunks/{c}:{i}-hap_chunks.txt",
        final_anchor = w_dir+"/anchors/{c}:{i}-anchor.txt",
        safe_lim = get_safe_lim,
        safe_len= get_safe_len
    resources:
    	runtime_hrs = 160,
    	mem_total_mb = 400000
    shell:
        """
        python3 Scripts/inv_locator.py -invs {input.inv} -ref CHM13 -gfa -gfa {input.gfa} -hapchunks {params.hap_chunks} -anchor {params.final_anchor} -flank 70000 -safe_len {params.safe_len} -limit 100000 -safe_len_limit {params.safe_lim} -exhaps HG03456#1,HG03456#2 > {output.node_sanity}
        """
        
#Step 2.1       
rule get_bubbles:
    input:
        gfa = path_to_gfas+"/{c}.gfa",
    output:
        out = path_to_gfas+"/bubbles/{c}-bubbles.json",
        log = path_to_gfas+"/bubbles/{c}.log"
    resources:
        runtime_hrs=10,
        mem_total_mb=30000
    shell:
        """
        {path_to_bubblegun} -g {input.gfa} --log_file {output.log} bchains --bubble_json {output.out} 
        """

#Step 2.2                 
rule get_non_rep_simple_bubbles:
    input:
        gfa = path_to_gfas+"/{c}.gfa",
        bubbles = path_to_gfas+'/bubbles/{c}-bubbles.json'
    output:
    	out = path_to_gfas+'/bubbles/{c}-nonrep-simple-bubbles.txt'
    resources:
        runtime_hrs=10,
        mem_total_mb=70000
    shell:
        """
        python3 {path_to_scripts}/nonrep_bubbles_simple.py -g {input.gfa}  -bubbles {input.bubbles} > {output.out}
        """  
#Step 2.3                 
rule get_non_rep_bubbles:
    input:
        gfa = path_to_gfas+"/{c}.gfa",
        bubbles = path_to_gfas+'/bubbles/{c}-bubbles.json'
    output:
    	out = path_to_gfas+'/bubbles/{c}-nonrep-all-bubbles.txt'
    resources:
        runtime_hrs=10,
        mem_total_mb=70000
    shell:
        """
        python3 {path_to_scripts}/nonrep_bubbles_all.py -g {input.gfa}  -bubbles {input.bubbles} > {output.out}   
        """      
#Step 3.1          
rule final_counts_simple:
    input:
    	bubbles = path_to_gfas+'/bubbles/{c}-nonrep-simple-bubbles.txt',
    	node_sanity = w_dir+"/node_sanity/{c}:{i}-node_sanity.txt"
    output:
        out_files = w_dir+"/final_tables/simple_bubbles/{c}:{i}.txt"
    params:
        hap_chunks = w_dir+"/hap_chunks/{c}:{i}-hap_chunks.txt",
    resources:
        runtime_hrs = 3,
        mem_total_mb = 50000      
    run:
        if os.path.isfile(params.hap_chunks):
            shell("python3 {path_to_scripts}/hapchunks_processing_allinone_bubbleinput.py -haps {params.hap_chunks} -bubbles {input.bubbles} -out {output.out_files}")
        else:
            shell("echo 'no haplotype chunks found' && touch {output.out_files}")
            
rule produce_table_all_contigs_simple_bubbles:
    input:
        final_tables_simple,
        node_lengths
    output:
        out = w_dir+"/summary_table/Summary_all_contigs_simple_bubbles.tsv"
    params:
        main = w_dir,
        specific = w_dir+ "final_tables/simple_bubbles/"
    resources:
        runtime_hrs=2,
        mem_total_mb=30000
    conda: "envs/r4.yaml"
    shell:
        """
        Rscript Scripts/Summary_table.R -m {params.main} -s {params.specific} -o {output.out}
        """
   
#Step 3.1      (alternative)     
rule final_counts_cont_contigs_simple:
    input:
    	bubbles = path_to_gfas+'/bubbles/{c}-nonrep-simple-bubbles.txt',
    	node_sanity = w_dir+"/node_sanity/{c}:{i}-node_sanity.txt"
    output:
        out_files = w_dir+"/final_tables_cont_contigs/simple_bubbles/{c}:{i}.txt"
    params:
        hap_chunks = w_dir+"/hap_chunks/{c}:{i}-hap_chunks.txt",
    resources:
        runtime_hrs = 3,
        mem_total_mb = 50000      
    run:
        if os.path.isfile(params.hap_chunks):
            shell("python3 {path_to_scripts}/hapchunks_processing_allinone_bubbleinput_cont_contigs.py -haps {params.hap_chunks} -bubbles {input.bubbles} -out {output.out_files}")
        else:
            shell("echo 'no haplotype chunks found' && touch {output.out_files}")

rule produce_table_cont_contigs_simple_bubbles:
    input:
        final_tables_cont_contigs_simple,
        node_lengths
    output:
        out = w_dir+"/summary_table/Summary_cont_contigs_simple_bubbles.tsv"
    params:
        main = w_dir,
        specific = w_dir+ "final_tables_cont_contigs/simple_bubbles/"
    resources:
        runtime_hrs=2,
        mem_total_mb=30000
    conda: "envs/r4.yaml"
    shell:
        """
        Rscript Scripts/Summary_table.R -m {params.main} -s {params.specific} -o {output.out}
        """

#Step 3.2           
rule final_counts_all:
    input:
    	bubbles = path_to_gfas+'/bubbles/{c}-nonrep-all-bubbles.txt',
    	node_sanity = w_dir+"/node_sanity/{c}:{i}-node_sanity.txt"
    output:
        out_files = w_dir+"/final_tables/all_bubbles/{c}:{i}.txt" 
    params:
        hap_chunks = w_dir+"/hap_chunks/{c}:{i}-hap_chunks.txt",
    resources:
        runtime_hrs = 3,
        mem_total_mb = 50000      
    run:
        if os.path.isfile(params.hap_chunks):
            shell("python3 {path_to_scripts}/hapchunks_processing_allinone_bubbleinput.py -haps {params.hap_chunks} -bubbles {input.bubbles} -out {output.out_files}")
        else:
            shell("echo 'no haplotype chunks found' && touch {output.out_files}")

rule produce_table_all::
    input:
        final_tables_all,
        node_lengths
    output:
        out = w_dir+"/summary_table/Summary_all_contigs_all_bubbles.tsv"
    params:
        main = w_dir,
        specific = w_dir+ "final_tables/all_bubbles/"
    resources:
        runtime_hrs=2,
        mem_total_mb=30000
    conda: "envs/r4.yaml"
    shell:
        """
        Rscript Scripts/Summary_table.R -m {params.main} -s {params.specific} -o {output.out}
        """

#Step 3.2      (alternative)     
rule final_counts_cont_contigs_all:
    input:
    	bubbles = path_to_gfas+'/bubbles/{c}-nonrep-all-bubbles.txt',
    	node_sanity = w_dir+"/node_sanity/{c}:{i}-node_sanity.txt"
    output:
        out_files = w_dir+"/final_tables_cont_contigs/all_bubbles/{c}:{i}.txt"
    params:
        hap_chunks = w_dir+"/hap_chunks/{c}:{i}-hap_chunks.txt",
    resources:
        runtime_hrs = 3,
        mem_total_mb = 50000      
    run:
        if os.path.isfile(params.hap_chunks):
            shell("python3 {path_to_scripts}/hapchunks_processing_allinone_bubbleinput_cont_contigs.py -haps {params.hap_chunks} -bubbles {input.bubbles} -out {output.out_files}")
        else:
            shell("echo 'no haplotype chunks found' && touch {output.out_files}")

rule produce_table_all::
    input:
        final_tables_cont_contigs_all,
        node_lengths
    output:
        out = w_dir+"/summary_table/Summary_cont_contigs_all_bubbles.tsv"
    params:
        main = w_dir,
        specific = w_dir+ "final_tables_cont_contigs/all_bubbles/"
    resources:
        runtime_hrs=2,
        mem_total_mb=30000
    conda: "envs/r4.yaml"
    shell:
        """
        Rscript Scripts/Summary_table.R -m {params.main} -s {params.specific} -o {output.out}
        """

#Independent step       
rule node_lengths:
    input:
        gfa = path_to_gfas+"/{c}.gfa",
    output:
        out =w_dir+"/node_lengths/{c}.txt",
    resources:
        runtime_hrs=10,
        mem_total_mb=30000
    shell:
        """
        python3 {path_to_scripts}/node_lengths.py -gfa {input.gfa} > {output.out}
        """ 
