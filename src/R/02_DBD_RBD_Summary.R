############################################################################
############################### Script 2 ###################################
############################################################################

## Function : Summarize DNA-binding domains (DBDs) and RNA-binding domains (RBDs)
#             based on DBRs and RBRs from Script 1.

## Inputs : 
#   (1) WD/path : Path for basic repo DIRectory
#   (1) project : Project folder (where files generated in script 1 are located)
#   (1) LncRNA data file : File + path from DIRectory for lncRNA bed file
#   (1) Gene data file : File + path from DIRectory for gene bed file

## Output : Generates two main groups of files within your project DIRectory.
#           Group 1 - TSV and Excel sheets for all DBDs + bed file with DBD coordinates.
#           Group 2 - TSV and Excel sheets for all RBDs + bed file with RBD coordinates.
#           Note: Script generates a 'noLists' version for both DBD and RBD files which
#                 significantly saves on space by removing the DBR_RBR list column from the file.

## Libraries and helper functions
library(tidyverse)
library(ggplot2)
library(openxlsx)

args = commandArgs(trailingOnly=TRUE)
DIR <- paste0(args[1])
#DIR <- '/projectnb/wax-dk/alan/2022_TriplexAnalysis'

source(paste0(DIR,"/src/helpers/functions.R"))

#### USER SPECIFIC CODE ####

param <- read_tsv(paste0(DIR,'/data/parameters.txt'))[1,]
#param <- read_tsv(paste0(DIR,'/data/parameters.txt'))[args[2],]

# Path to bed files 
if(param$Region_Bed == "NULL"){
  region_bed_file <- list.files(path = paste0(DIR,'/data/',param$Input),pattern = 'regions.bed')
}else{
  region_bed_file <- param$Region_Bed
}
if(param$LncRNA_Bed == "NULL"){
  lncRNA_bed_file <- '00_REF/09_RefSeqLncRNA76k_longestIsoform_geneOnly_mm10.bed'
}else{
  lncRNA_bed_file <- param$LncRNA_Bed
}

## Parameters for filtering low scoring TTS prior to grouping them into regions
TFO_ADJUSTED_VALUE = ifelse(param$TFO_ADJUSTED_VALUE == "NULL",NA,param$TFO_ADJUSTED_VALUE) #20.01
TFO_COUNT = ifelse(param$TFO_COUNT == "NULL",NA,param$TFO_COUNT) #2 

# Number of DBR
REGION_MIN <-  ifelse(param$REGION_MIN == "NULL",NA,param$REGION_MIN) #5

# Generated diagnostic figures for DBR and RBR regions
DIAGNOSTIC_FIGURES = param$DIAGNOSTIC_FIGURES

#############################################################
#### Script shouldn't need modification below this point ####
#############################################################

# Bed file with lncRNA coordinates
lnc_data <-read_tsv(paste0(DIR,'/data/',lncRNA_bed_file),col_names = F)
lncRNA <- lnc_data[,c(4,1,2)]
colnames(lncRNA) <- c("lnc","chr","start")

dir.create(paste0(DIR,"/output/",param$Output,"/02"),recursive = T,showWarnings = F)

#### Calculate DBDs ####

# Read in DBRs
df <- read_tsv(paste0(DIR,"/output/",param$Output,"/01/alllncRNA_DBRs.txt"))

## Diagnostic figures
if(DIAGNOSTIC_FIGURES == T){
  # Distribution of average triplex formation scores (Adj1)
  ggplot(df,aes(x=DBR_TFO_abj2)) +
    geom_histogram(bins = 500) +
    geom_vline(xintercept = c(1,1.2,1.5),
               colour=c('red','blue','green')) +
    coord_trans(x="log2") +
    scale_y_continuous(trans = 'log2') +
    labs(x="TFO Adjusted Value (abj2)",
         y="Count") +
    theme_bw()
  
  # Distribution of average DBR scores
  ggplot(df,aes(x=DBR_score)) +
    geom_histogram(bins=500) +
    # geom_vline(xintercept = c(1,1.2,1.5),
    #            colour=c('red','blue','green')) +
    #coord_trans(x="log2") +
    scale_y_continuous(trans = 'log2') +
    labs(x="Average TFO Score",
         y="Count") +
    theme_bw()
  
  ggplot(df,aes(x=DBR_RBR_span)) +
    geom_histogram(bins=length(unique(df$DBR_RBR_span)),binwidth = 1) +
    # geom_vline(xintercept = c(1,1.2,1.5),
    #            colour=c('red','blue','green')) +
    #coord_trans(x="log2") +
    scale_y_continuous(trans = 'log2') +
    labs(x="DBR Span",
         y="Count") +
    theme_bw()
}

## Filter DBRs with only a few TTSs
if (!is.na(TFO_ADJUSTED_VALUE) | !is.na(TFO_COUNT)){
  if(is.na(TFO_ADJUSTED_VALUE)){
    df2 <- df[df$DBR_TFO_count <= TFO_COUNT,]
  }else{
    if(is.na(TFO_COUNT)){
      df2 <- df[df$DBR_TFO_abj2 >= TFO_ADJUSTED_VALUE,]
    }else{
      df2 <- df[df$DBR_TFO_abj2 >= TFO_ADJUSTED_VALUE & df$DBR_TFO_count >= TFO_COUNT,]
    }
  }
}else{df2 <- df}

## Calculate DBDs
f3 <- df2 %>% 
  arrange(lnc,DBR_end) %>% 
  mutate(lag_end = lag(DBR_end),
         lag_lnc = lag(lnc)) %>% 
  replace_na(list(lag_end = 10000000000)) %>%
  mutate(overlap = DBR_start-lag_end,
         lnc_overlap = as.numeric(lnc == lag_lnc)) %>%
  replace_na(list(lnc_overlap = 1)) %>%
  mutate(DBD = find_domain2(overlap,lnc_overlap))

# Generate alternate file with DBRs that contain DBD id
f3_alt <- f3 %>%
  relocate(DBD,.after=gname) %>%
  select(c(-lag_end,-lag_lnc,-overlap,-lnc_overlap))
write_tsv(f3_alt,paste0(DIR,"/output/",param$Output,"/02/alllncRNA_DBRs_withDBDs.txt"))

# Some of the reference bed files have duplicate entries. This removes those duplicates.
lncRNA <- lncRNA[!duplicated(lncRNA$lnc),]
# Merge lncRNA coordinates with DBR information
f3 <- left_join(f3,lncRNA,by=c('lnc'='lnc'))

# Summarize DBDs
f4 <- f3 %>%
  group_by(lnc,DBD) %>%
  summarise(DBD_chr=chr[1],
            DBD_start=min(DBR_start),
            DBD_end=max(DBR_end),
            DBD_start_mm10=min(DBR_start)+min(start),
            DBD_end_mm10=max(DBR_end)+min(start),
            DBD_span=DBD_end - DBD_start,
            DBD_geneRegion_Count=n_distinct(gname),
            DBD_DBR_Count=n_distinct(DBR_ID),
            DBD_DBR_RBR_Count=sum(DBR_RBR_count),
            DBD_TFO_Count=sum(DBR_TFO_count),
            DBD_TFO_Avg=mean(DBR_TFO_count),
            DBD_TFO_Max=max(DBR_TFO_count),
            DBD_abj1_Total=sum(DBR_TFO_abj1),
            DBD_abj1_Avg=mean(DBR_TFO_abj1),
            DBD_abj1_Max=max(DBR_TFO_abj1),
            DBD_abj1_Avg_byRegion=DBD_abj1_Avg/DBD_geneRegion_Count,
            DBD_abj2_Total=sum(DBR_TFO_abj2),
            DBD_abj2_Avg=mean(DBR_TFO_abj2),
            DBD_abj2_Avg_byRegion=DBD_abj2_Avg/DBD_geneRegion_Count,
            DBD_score = mean(DBR_score),
            DBD_error_rate = mean(DBR_error_rate),
            DBD_gc_error = mean(DBR_gc_error),
            DBD_geneRegion_list=paste0(unique(gname),collapse = ","),
            DBD_DBR_list=paste0(DBR_ID,collapse = ","),
            DBD_BR_list=paste0(BR_list,collapse = ","))

## Minimum DBR per DBD Threshold 
if(!is.na(REGION_MIN)){
  f4_filter <- f4[f4$DBD_DBR_Count >= REGION_MIN,]
}else{f4_filter <- f4}

if(DIAGNOSTIC_FIGURES == T){
  f4_filter$geneRegion_ls <- strsplit(f4_filter$DBD_geneRegion_list,",")
  diag_df <- data.frame(lncRNA=unique(f4_filter$lnc))
  diag_df$unique_genes <- rep(1,times=length(diag_df$lncRNA))
  diag_df$max_interaction <- rep(1,times=length(diag_df$lncRNA))
  for(i in 1:length(diag_df$lncRNA)){
    diag_df$unique_genes[i] <- length(unique(unlist(f4_filter$geneRegion_ls[grep(pattern = diag_df$lncRNA[i],f4_filter$lnc)])))
    diag_df$max_interaction[i] <- max(table(unlist(f4_filter$geneRegion_ls[grep(pattern = diag_df$lncRNA[i],f4_filter$lnc)])))
  }
  
  ## DBD diagnostic plots
  p1 <- data.frame(table(f4_filter$lnc)) %>%
    arrange(Freq) %>%
    mutate(name=factor(Var1, levels=Var1)) %>%
    ggplot(aes(x=name,y=Freq)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=1, color="orange") +
    coord_flip() +
    theme_bw() +
    labs(x="LncRNA",y="Number of DBDs")
  
  ggsave(plot=p1,
         width = 10,
         height = 25,
         dpi = 800, 
         filename=paste0(DIR,"/output/",param$Output,"/02/DBD_NumDBDbyLncRNA.pdf"))
  
  p1 <- diag_df %>%
    arrange(unique_genes) %>%
    mutate(name=factor(lncRNA, levels=lncRNA)) %>%
    ggplot(aes(x=name,y=unique_genes)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=1, color="orange") +
    coord_flip() +
    theme_bw() +
    labs(x="LncRNA",y="Number of Unique Regions")
  ggsave(plot=p1,
         width = 10,
         height = 25,
         dpi = 800, 
         filename=paste0(DIR,"/output/",param$Output,"/02/DBD_NumRegionsbyLncRNA.pdf"))
  
  p1 <- diag_df %>%
    arrange(max_interaction) %>%
    mutate(name=factor(lncRNA, levels=lncRNA)) %>%
    ggplot(aes(x=name,y=max_interaction)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=1, color="orange") +
    coord_flip() +
    theme_bw() +
    labs(x="LncRNA",y="Max number of DBRs with single region")
  p1
  ggsave(plot=p1,
         width = 10,
         height = 25,
         dpi = 800, 
         filename=paste0(DIR,"/output/",param$Output,"/02/DBD_MaxRegionDBRsbyLncRNA.pdf"))
}

## Generate files
list_remove <- c('DBD_geneRegion_list','DBD_DBR_list','DBD_BR_list')
write_tsv(f4_filter,file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_DBDs.txt"))
write_tsv(re_assign(f4_filter,list_remove),file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_DBDs_noLists.txt"))
write.xlsx(f4_filter,file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_DBDs.xlsx"),keepNA=T)
write.xlsx(re_assign(f4_filter,list_remove),file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_DBDs_noLists.xlsx"),keepNA=T)
dbd_bed <- data.frame(chr=f4_filter$DBD_chr,start=f4_filter$DBD_start_mm10,end=f4_filter$DBD_end_mm10,dbr=paste0(f4_filter$lnc,"_",f4_filter$DBD))
write_tsv(dbd_bed,file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_DBDs_mm10coordinates.bed"),col_names = F)

#### Calculate RBDs ####
df <- read_tsv(paste0(DIR,"/output/",param$Output,"/01/alllncRNA_RBRs.txt"),col_names = T)

## Filter DBRs with only a few TTSs
if (!is.na(TFO_ADJUSTED_VALUE) | !is.na(TFO_COUNT)){
  if(is.na(TFO_ADJUSTED_VALUE)){
    df2 <- df[df$RBR_TFO_count_total <= TFO_COUNT,]
  }else{
    if(is.na(TFO_COUNT)){
      df2 <- df[df$RBR_TFO_abj1 >= TFO_ADJUSTED_VALUE,]
    }else{
      df2 <- df[df$RBR_TFO_abj1 >= TFO_ADJUSTED_VALUE & df$RBR_TFO_count_total >= TFO_COUNT,]
    }
  }
}else{df2 <- df}

# Read in gene
genes <- read_tsv(paste0(DIR,'/data/',param$Input,'/',region_bed_file),col_names = F)
genes <- genes[,c(4,1,2)]
colnames(lncRNA) <- c("gene","chr","start")

f3 <- df2 %>% 
  arrange(gname,RBR_end) %>% 
  mutate(lag_end = lag(RBR_end),
         lag_gname = lag(gname)) %>% 
  replace_na(list(lag_end = 10000000000)) %>%
  mutate(overlap = RBR_start-lag_end,
         gname_overlap = as.numeric(gname == lag_gname)) %>%
  replace_na(list(gname_overlap = 1)) %>%
  mutate(RBD = find_domain2(overlap,gname_overlap))

f3_alt <- f3 %>%
  relocate(RBD,.after=gname) %>%
  select(c(-lag_end,-lag_gname,-overlap,-gname_overlap))
write_tsv(f3_alt,paste0(DIR,"/output/",param$Output,"/02/alllncRNA_RBRs_withRBDs.txt"))

f4 <- f3 %>%
  group_by(gname,RBD) %>%
  summarize(RBD_chr=RBR_chr[1],
            RBD_start=min(RBR_start),
            RBD_end=max(RBR_end),
            RBD_span=RBD_end - RBD_start,
            RBD_lncRNA_Count=n_distinct(lnc),
            RBD_RBR_Count=n_distinct(RBR_ID),
            RBD_DBR_RBR_Count=sum(RBR_DBR_count),
            RBD_TFO_Total=sum(RBR_TFO_count_total),
            RBD_TFO_Avg=mean(RBR_TFO_count_total),
            RBD_TFO_Max=max(RBR_TFO_count_total),
            RBD_abj1_Total=sum(RBR_TFO_abj1),
            RBD_abj1_Avg=mean(RBR_TFO_abj1),
            RBD_abj1_Max=max(RBR_TFO_abj1),
            RBD_abj2_Total=mean(RBR_TFO_abj2),
            RBD_abj2_Avg=mean(RBR_TFO_abj2),
            RBD_score = mean(RBR_score),
            RBD_error_rate = mean(RBR_error_rate),
            RBD_gc_error = mean(RBR_gc_error),
            RBD_lncRNA_list=paste0(unique(lnc),collapse = ","),
            RBD_RBR_list=paste0(RBR_ID,collapse = ","),
            RBD_BR_list=paste0(BR_list,collapse = ","))

## Minimum DBR per DBD Threshold 
if(!is.na(REGION_MIN)){
  f4_filter <- f4[f4$RBD_RBR_Count >= REGION_MIN,]
}else{f4_filter <- f4}

list_remove <- c('RBD_geneRegion_list','RBD_RBR_list','RBD_BR_list')
write_tsv(f4_filter,file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_RBDs.txt"))
write_tsv(re_assign(f4_filter,list_remove),file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_RBDs_noLists.txt"))
write.xlsx(f4_filter,file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_RBDs.xlsx"),keepNA=T)
write.xlsx(re_assign(f4_filter,list_remove),file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_RBDs_noLists.xlsx"),keepNA=T)
rbd_bed <- data.frame(chr=f4_filter$RBD_chr,start=f4_filter$RBD_start,end=f4_filter$RBD_end,dbr=paste0(f4_filter$gname,"_",f4_filter$RBD))
write_tsv(rbd_bed,file = paste0(DIR,"/output/",param$Output,"/02/alllncRNA_RBDs_mm10coordinates.bed"),col_names = F)

### Meta data text file ###
meta_list <- list("#### Input Data #####",
                  paste0("# Region Data File: ", region_bed_file),
                  paste0("# Lnc Data File: ", lncRNA_bed_file),
                  "#### Parameters #####",
                  paste0("# Adjusted TFO Threshold (adj1; TFO_ADJUSTED_VALUE): ",TFO_ADJUSTED_VALUE),
                  paste0("# Minimum TFO Count (TFO_COUNT): ",TFO_COUNT),
                  paste0("# Minimum number of DBR in a DBD (REGION_MIN): ",REGION_MIN ),
                  paste0("# Stored Diagnostic Plots (DIAGNOSTIC_FIGURES): ",DIAGNOSTIC_FIGURES))
writeLines(unlist(meta_list),paste0(DIR,"/output/",param$Output,'/02/metadata.txt'))
