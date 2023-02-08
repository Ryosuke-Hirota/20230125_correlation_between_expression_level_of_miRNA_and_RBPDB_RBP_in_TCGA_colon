# This script is to investigate correlation between expression level of miRNA and RBP and draw plot about it
# made 2023/01/25

library(stringr)
library(ggplot2)

# set function for calculating outlier
cal.outlier <-function(x){
  q <-as.numeric(quantile(x))
  iqr <-IQR(x)
  outlier1 <-q[2]-iqr*1.5
  outlier2 <-q[4]+iqr*1.5
  outliers <-append(outlier1,outlier2)
  return(outliers)
}

# set function for calculating minimum value of expression level
find.min <-function(x,y){
  for (i in 1:y) {
    m <-x[x[,i]!=0,i]
    mv <-min(m)
    if(i==1){
      mvs <-mv
    }else{
      mvs <-append(mvs,mv)
    }}
  return(min(mvs))
}

# make new directory
setwd("C:/Rdata")
dir.create("20230125_correlation_between_expression_level_of_miRNA_and_RBPDB_RBP_in_TCGA_colon")

# import table of gene counts
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230117_TCGA_colon_gene_counts"
setwd("C:/Rdata/20230105_TCGA_colon_miRNA_quantification")
miRNA.quant.table <-read.table("table_of_TCGA_colon_miRNA_quantifications.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# calculate minimum value of miRNA expression level
m.table <-miRNA.quant.table[,-1]
min.miRNA <-find.min(m.table,467)

# import table of gene counts
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230117_TCGA_colon_gene_counts"
setwd("C:/Rdata/20230117_TCGA_colon_gene_counts")
gene.counts.table <-read.table("table_of_TCGA_colon_gene_counts.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# calculate minimum value of gene expression level
g.table <-gene.counts.table[,c(-1,-2)]
min.gene <-find.min(g.table,293)

# import list of RBPDB RBP
# this list is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\refereces"
setwd("C:/Rdata/references")
RBPDB.RBP.list <-read.csv("RBPDB_v1.3.1_proteins_human_2012-11-21.csv",sep=",",header = F,stringsAsFactors = F)
RBPDB.RBP <-unique(RBPDB.RBP.list[RBPDB.RBP.list[,5]!="",5])

# import list of RBP2GO RBP
# this list is located at "https://www.dropbox.com/home/Okamura%20Lab%20share%20folder/Hirota/results_and_matterials/RBP2GO's_RBPs_contained_to_CCLE"
setwd("C:/Rdata/RBP2GO's_RBPs_contained_to_CCLE")
RBP2GO.RBP.list <-read.table("all_RBP2GO_RBPs_which_merged_with_CCLE.txt",sep="\t",header = T,stringsAsFactors = F)

# import CCLE RBP list
# this list is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
df <-read.table("RNAseq_of_RBP.txt",sep="\t",header = T,stringsAsFactors = F)
ccle.RBP <-colnames(df)
ccle.RBP <-ccle.RBP[c(-1,-412)]
m1 <-match(ccle.RBP,RBP2GO.RBP.list[,2])
RBP <-RBP2GO.RBP.list[m1,4]
m2 <-match(RBP,gene.counts.table[,2])
RBP <-gene.counts.table[m2,2]
RBP <-sort(as.character(unique(na.omit(RBP))))

# import correspondence table of TCGA file
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230117_TCGA_colon_gene_counts"
setwd("C:/Rdata/20230116_TCGA_colon_make_correspond_gene_counts_files")
cor.table <-read.table("TCGA_colon_gene_counts_correspondence_table.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# import list of CCLE pri-miRNA
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_transcript_that_intersect_with_miRNA"
setwd("C:/Rdata/20230110_TCGA_transcript_that_intersect_with_miRNA")
mir.list <-read.table("list_of_TCGA_pri-miRNA_that_exsit_in_CCLE.txt",sep="\t",header = T,stringsAsFactors = F)
miRNA <-unique(mir.list[,1])

# make empty summary
sm <-as.data.frame(matrix(nrow = 278*402,ncol = 5))
colnames(sm) <-c("miRNA","RBP","r","p.value","number_of_sample")

# calculate correlation coefficient and draw plot
for (i in 1:length(miRNA)) {
  for (k in 1:length(RBP)) {
    
    if(k==1){
      setwd("C:/Rdata/20230125_correlation_between_expression_level_of_miRNA_and_RBPDB_RBP_in_TCGA_colon")
      dir.create(miRNA[i])
    }
    
    # extract expression levelof certain miRNA from summarized table
    miRNA.df <-miRNA.quant.table[miRNA.quant.table[,1]==miRNA[i],]
    miRNA.df <-as.data.frame(t(miRNA.df),stringsAsFactors = F)
    miRNA.df[,2] <-rownames(miRNA.df)
    miRNA.df <-miRNA.df[-1,]
    rownames(miRNA.df) <-NULL
    colnames(miRNA.df) <-c(miRNA[i],"miRNA_file_name")
    
    # extract expression level of a certain gene from summarized table  
    gene.counts.df <-gene.counts.table[gene.counts.table[,2]==RBP[k],]
    gene.counts.df <-as.data.frame(t(gene.counts.df),stringsAsFactors = F)
    gene.counts.df[,2] <-rownames(gene.counts.df)
    gene.counts.df <-gene.counts.df[c(-1:-2),]
    rownames(gene.counts.df) <-NULL
    colnames(gene.counts.df) <-c(RBP[k],"filename")
    
    # add name
    f <-match(gene.counts.df[,2],cor.table[,2])
    gene.counts.df[,3] <-cor.table[f,6]
    colnames(gene.counts.df)[3] <-"miRNA_file_name"
    
    # merge residual dataframe and gene expression dataframe
    merged.df <-merge(miRNA.df,gene.counts.df,by="miRNA_file_name")
    merged.df <-merged.df[,c(2,3,1,4)]
    
    # if a sample hasn't miRNA expression level, substitute minmum value 
    merged.df[,1] <-ifelse(merged.df[,1]==0,min.miRNA,merged.df[,1])
    
    # if a sample hasn't RBP expression level, substitute minmum value 
    merged.df[,2] <-ifelse(merged.df[,2]==0,min.gene,merged.df[,2])
    
    # normalize RBP expression
    merged.df[,1] <-as.numeric(merged.df[,1])
    merged.df[,2] <-as.numeric(merged.df[,2])
    merged.df[,1] <-log2(merged.df[,1])
    merged.df[,2] <-log2(merged.df[,2])
    
    # calculate outliers
    m.outlier <-cal.outlier(merged.df[,1])
    p.outlier <-cal.outlier(merged.df[,2])
    
    # remove outliers
    merged.df <-merged.df[merged.df[,1]>m.outlier[1]&merged.df[,1]<m.outlier[2],]
    merged.df <-merged.df[merged.df[,2]>p.outlier[1]&merged.df[,2]<p.outlier[2],]
    
    # calculate correlation coefficient
    r <-try(cor.test(merged.df[,1],merged.df[,2],method = "spearman"),silent = T)
    
    if(class(r)!="try-error"){
      # draw plot
      p <-ggplot(data = merged.df,aes(x=merged.df[,1],y=merged.df[,2]))+
        geom_point(color="blue")+
        geom_smooth(data=merged.df,mapping = aes(x=merged.df[,1],y=merged.df[,2]),method="lm",formula='y~x',se=FALSE,colour="black",size=0.5)+
        labs(title=paste0("R =",signif(r$estimate,3),", p = ",signif(r$p.value,3),", n = ",nrow(merged.df)),x=miRNA[i],
             y=RBP[k])+ 
        theme_bw()+
        theme(legend.background = element_rect(fill = "white", colour = "black"))
      
      setwd(paste0("C:/Rdata/20230125_correlation_between_expression_level_of_miRNA_and_RBPBD_RBP_in_TCGA_colon/",miRNA[i]))
      ggsave(filename=paste0("plot_of_correlation_between_",miRNA[i],"_and_",RBP[k],".pdf"),plot = p)
      
      # write summary
      sm[402*(i-1)+k,1] <-miRNA[i]
      sm[402*(i-1)+k,2] <-RBP[k]
      sm[402*(i-1)+k,3] <-r$estimate
      sm[402*(i-1)+k,4] <-r$p.value
      sm[402*(i-1)+k,5] <-nrow(merged.df)
    }else{
      # write summary
      sm[402*(i-1)+k,1] <-miRNA[i]
      sm[402*(i-1)+k,2] <-RBP[k]
      sm[402*(i-1)+k,3] <-NA
      sm[402*(i-1)+k,4] <-NA
      sm[402*(i-1)+k,5] <-nrow(merged.df)
    }
  }}

# write summary table
setwd("C:/Rdata/20230125_correlation_between_expression_level_of_miRNA_and_RBPDB_RBP_in_TCGA_colon")
write.table(sm,"table_of_correlation_between_expression_level_of_miRNA_and_RBPDB_RBP_in_TCGA_colon.txt",sep="\t",row.names = F,quote = F)

