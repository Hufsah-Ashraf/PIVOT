oldw <- getOption("warn")
options(warn = -1)

suppressMessages(library("optparse"))
library(data.table)
library(dplyr)
library(stringr)

option_list = list( 
  make_option(c("-m", "--mainfolder"), type="character", default=NULL, 
              help="main folder with all the outputs", metavar="character"),
  make_option(c("-s", "--specificfolder"), type="character", default=NULL, 
              help="specific input w.r.t to bubbles or broken contigs", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="output file", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
main_folder<-opt$mainfolder
files_link<-opt$specificfolder

files <- list.files(files_link,pattern = "txt")
full_df<-data.frame()
warning(
  for (f in 1:length(files)){
    df<-fread(paste0(files_link,files[f]))
    chrom<- str_split_fixed(files[f], pattern = ':', n=Inf)[1]
    node_lengths<-fread(paste0(main_folder,'node_lengths/',chrom, '.txt'))
    colnames(node_lengths)<-c('node', 'node_length')
    df<-left_join(df, node_lengths)
    df$sum<-rowSums(df[,c(2:5)])
    tidf<-df[(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1),]
    
    sample_list <- list()
    if (nrow(tidf)>=1){
      for (t in 1:nrow(tidf)) {
        tidf[t,'min'] <-min(tidf[t,'FA'], tidf[t,'FP'], tidf[t,'RA'], tidf[t,'RP'])
        tidf[t,'min_haps']<-''
        if (tidf[t,'FA'] == tidf[t,'min'] ) {
          tidf[t,'min_haps']<-paste0(tidf[t,'min_haps'],',FA')
          #tidf[t,'FA_to_remove']<-tidf[t,'FA_samples']
          #FA_samples<-c(FA_samples, str_split_fixed(tidf[t,'FA_samples'],pattern = ',',n=Inf))
          
        }
        if (tidf[t,'FP'] == tidf[t,'min'] ) {
          tidf[t,'min_haps']<-paste0(tidf[t,'min_haps'],',FP')
          #tidf[t,'FP_to_remove']<-tidf[t,'FP_samples']
          #FP_samples<-c(FP_samples, str_split_fixed(tidf[t,'FP_samples'],pattern = ',',n=Inf))
        }
        if (tidf[t,'RA'] == tidf[t,'min'] ) {
          tidf[t,'min_haps']<-paste0(tidf[t,'min_haps'],',RA')
          #tidf[t,'RA_to_remove']<-tidf[t,'RA_samples']
          #RA_samples<-c(RA_samples, str_split_fixed(tidf[t,'RA_samples'], pattern = ',',n=Inf))
        }
        if (tidf[t,'RP'] == tidf[t,'min']) {
          tidf[t,'min_haps']<-paste0(tidf[t,'min_haps'],',RP')
          #tidf[t,'RP_to_remove']<-tidf[t,'RP_samples']
          #RP_samples<-c(RP_samples, str_split_fixed(tidf[t,'RP_samples'], pattern = ',',n=Inf))
        }
        
        haps_to_cons<-str_split_fixed(tidf[t,'min_haps'], pattern = ',', n=Inf)[2:length(str_split_fixed(tidf[t,'min_haps'], pattern = ',', n=Inf))]
        len<-NA
        for (h in haps_to_cons){
          if (is.na(len)){
            k<-paste0(as.character(h),'_samples')
            len <- length(unique(c(sample_list,str_split_fixed(tidf[t,..k], pattern=',', n=Inf))))
            orientation<-paste0(as.character(h),'_samples')
          }else{
            k<-paste0(as.character(h),'_samples')
            if (length(unique(c(sample_list,str_split_fixed(tidf[t,..k], pattern=',', n=Inf)))) < len) {
              len <- length(unique(c(sample_list,str_split_fixed(tidf[t,..k], pattern=',', n=Inf))))
              orientation<-paste0(as.character(h),'_samples')
              
            }
          }
        }
        sample_list<-unique(c(sample_list,str_split_fixed(tidf[t,..orientation], pattern=',', n=Inf)))
      }
    }

    new_df<- data.frame(Inv = str_split_fixed(files[f], pattern = '.txt', n=Inf)[1] )
    new_df$haplotypes<-max(df$sum)
    new_df$nodes_considered<-length(df$node)
    new_df$length_nodes_considered<-sum(df$node_length)
    new_df$num_tinodes<-length(which(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1))
    new_df$tinodes_len_sum<-sum(df[(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1),]$node_length)
    new_df$tiSNPs<- length(which(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1 & df$node_length==1))
    new_df$ti_smallvar<- length(which(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1 & df$node_length>1 & df$node_length<50))
    new_df$ti_smallvar_len<- sum(df[(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1 & df$node_length>1 & df$node_length<50),]$node_length)
    new_df$ti_structvar<- length(which(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1 & df$node_length>=50))
    new_df$ti_structvar_len<- sum(df[(df$FA>=1 & df$FP >=1 & df$RA >=1 & df$RP >=1 & df$node_length>=50),]$node_length)
    new_df$broken_contigs<-as.numeric(df[1,'broken_contigs'])
    new_df$unbroken_contigs<-as.numeric(df[1,'unbroken_contigs'])
    new_df$min_haps_to_remove<-length(sample_list)
    
    #new_df$min_haps_to_remove<- min(length(FA_samples),length(FP_samples), length(RA_samples), length(RP_samples))
    new_df$num_tinodes_2<-length(which(df$FA>=2 & df$FP >=2 & df$RA >=2 & df$RP >=2))
    new_df$tinodes_len_sum_2<-sum(df[(df$FA>=2 & df$FP >=2 & df$RA >=2 & df$RP >=2),]$node_length)
    
    new_df$num_tinodes_3<-length(which(df$FA>=3 & df$FP >=3 & df$RA >=3 & df$RP >=3))
    new_df$tinodes_len_sum_3<-sum(df[(df$FA>=3 & df$FP >=3 & df$RA >=3 & df$RP >=3),]$node_length)
    
    new_df$num_tinodes_4<-length(which(df$FA>3 & df$FP >3 & df$RA >3 & df$RP >3))
    new_df$tinodes_len_sum_4<-sum(df[(df$FA>3 & df$FP >3 & df$RA >3 & df$RP >3),]$node_length)
    new_df$max_min<-max(tidf$min)#the maximum value seen as min in the 4 cell table across all tinodes
    inv<-str_split_fixed(files[f], pattern = '.txt', n=Inf)[1]
    ###since I am not storing the exact anchor value anywhere except the hapchunks file, it is not possible to extract
    #it if the safe nodes didn't help, even if the anchor was found. So, for now I am implementing it in a way that
    # it just checks if the anchor was found or not but it should be updated in the anchor finding script to store the value in the anchor file
    anchor<-fread(paste0(main_folder,'anchor/',inv,'-anchor.txt'), nrows = 1)
    if(nrow(anchor)>0) {
      new_df$anchor_region<-as.character(!is.na(anchor[1,'anchor_path']))
    }
    
    full_df<-rbind(full_df,new_df)
  })
colnames(full_df)<-c(colnames(full_df)[1:ncol(full_df)-1],'anchor')
write.table(full_df, opt$output_file, col.names = T, row.names = F, quote = F, sep = '\t')
