# PIVOT: Pangenome graph based analysis of inversion toggling
PIVOT is an extended version of the [tiSNPs-based method for inversion recurrence detection](https://doi.org/10.1016/j.cell.2022.04.017) and operates within a pangenome-graph framework. Currently, PIVOT only works with [PGGB](https://www.nature.com/articles/s41592-024-02430-3) graphs. The paths (P lines) in the graph should follow the pattern sample#haplotype#contigid e.g. HG00096#1#chrY_random000000.

This workflow uses [Snakemake](https://bitbucket.org/snakemake/snakemake) to execute all steps in order. Before running the pipeline, the file paths should be added to Snakefile.json:

  1. **working_dir**: the folder where all the output files would be written 
  2. **path_to_invs**: the folder containing all the inversion regions to be tested. Each region should be specified in a tab-separate (.bed) file stating chromosome, inversion start and inversion end (no header) and the file should be named as {chrom}:{inversion start}-{inversion end}.bed e.g. chr1:13104252-13122521.bed.
  3. **path_to_graphs**: the folder with pggb graphs (.gfa) files for each chromosomes. The files should be named as {chromosomes}.gfa. Alternatively the naming pattern should be adjusted in the Snakefile.
  4. PIVOT uses [Bubblegun](https://github.com/fawaz-dabbaghieh/bubble_gun.git) for detecting bubbles in the graph. So, before running the pipeline it should be installed on the system and the path should be provided in **bubble_detection**
  
By default, PIVOT considers CHM13 as the reference. In case another haplotype is to be used as the reference, the **-ref** parameter in rule **find_anchors** needs to be adjusted.
If certain haplotypes in the graph need to be ignored, all of them need to be specified (comma-separated) in the parameter **-exhaps** in rule **find_anchors** as sample#haplotype e.g. HG03456#1,HG03456#2. 
The length of the flanking inverted repeat should be specified by **-flank** (default:70000bp). If an inverted repeat is not found in that range PIVOT keeps increasing the search window by 10kbp until the repeats are found or a length limit, specified under **-limit** (default:100000bp) is reached. 
The minimum length for safe nodes can be specified using **-safe_len**. By default it is 30bp for autosomes and 15bp for sex chromosomes. In certain cases, no node fulfilling the length criteria is "safe" so PIVOT keeps lowering the length by 5bp until safe nodes are found or the lowest allowed length, specified under **-safe_len_limit** (default: 15bp for autosomes and 5bp for sex chromosomes) is reached.

