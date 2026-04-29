# PIVOT: Pangenome graph based analysis of inversion toggling

This workflow uses [Snakemake](https://bitbucket.org/snakemake/snakemake) to execute all steps in order. Before running the pipeline, the file paths should be added to Snakefile.json:

  1. **path_to_graphs**: the folder with alignmnet (.bam) files for each sample. The folder structure should be "path_to_bams/{sample}/{one bam file per cell for this sample}"
  
Pivot works only with pggb graphs with paths named as sample#haplotype#contigid e.g. HG00096#1#chrY_random000000.
By default, PIVOT considers CHM13 as the reference. In case another haplotype is to be used as reference, the **-ref** parameter in rule **find_anchors** needs to be adjusted
If certain haplotypes in the graph need to be ignored, all of them need to be specified (comma-separated) in the parameter **-exhaps** in rule **find_anchors** as sample#haplotype e.g. HG03456#1,HG03456#2
https://github.com/fawaz-dabbaghieh/bubble_gun.git
