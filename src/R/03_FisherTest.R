############################################################################
############################### Script 03 ##################################
############################################################################

## Function : Rank and summarize individual Gene-lncRNA interactions

## Inputs : 
#   (1) WD/path : Path for basic repo directory
#   (2) project : Project folder (where files generated in script 1 are located)
#   (3) LncRNA data file : File + path from directory for lncRNA bed file
#   (4) Gene data file : File + path from directory for gene bed file

## Output : Generates two main groups of files within your project directory.
#           Group 1 - TSV and Excel sheets for all DBDs + bed file with DBD coordinates.
#           Group 2 - TSV and Excel sheets for all RBDs + bed file with RBD coordinates.
#           Note: Script generates a 'noLists' version for both DBD and RBD files which
#                 significantly saves on space by removing the DBR_RBR list column from the file.

## Libraries and helper functions
library(tidyverse)
library(cowplot)
source("/projectnb/wax-dk/alan/2022_TDFAnalysislnRNA/src/helpers/functions.R")

# Set main directory path and project folder name
dir <- "/projectnb/wax-dk/alan/2022_TDFAnalysislnRNA/"
project <- "Investigation10_SB_ATAC/"

ATAC_regions <- 'ATAC_regions.bed'
ATAC_Categories <- 'ATAC_region_treatment_Final.tsv'  

regions <- read_tsv(paste0(dir,"data/",project,ATAC_regions),col_names = F)
categories <- read_tsv(paste0(dir,"data/",project,ATAC_Categories),col_names = T)

categories_labels <- c('ATAC_Region','Categories') # Select column names from 'categories' that you would like to keep. 
# Note the region name that matches the gname in the dbr file should be the first item in vector.

#############################################################
#### Script shouldn't need modification below this point ####
#############################################################

dir.create(paste0(dir,"output/",project,"04"),recursive = T,showWarnings = F)

dbr <- read_tsv(paste0(dir,"output/",project,"/02/alllncRNA_DBRs_withDBDs.txt"),col_names = T)
rbr <- read_tsv(paste0(dir,"output/",project,"/02/alllncRNA_RBRs_withRBDs.txt"),col_names = T)


dbr_cat <- left_join(dbr,categories[,categories_labels],by=c('gname' = categories_labels[1]))
regions_cat <- left_join(regions,categories[,categories_labels],by=c('X4' = categories_labels[1]))

dbr_cat_summary <- dbr_cat %>%
  distinct(gname,.keep_all = TRUE) %>%
  count(Categories) %>%
  rename(Triplex_ATAC_Regions=n)

dbr_cat_summary <- dbr_cat %>%
  distinct(gname,.keep_all = TRUE) %>%
  group_by(Categories) %>%
  summarise(Triplex_ATAC_Regions=n(),
            DBR_span_mean=mean(DBR_span),
            DBR_span_max=max(DBR_span),
            DBR_TFO_count_mean = mean(DBR_TFO_count),
            DBR_TFO_count_max = max(DBR_TFO_count),
            Adj1_TFO_mean=mean(DBR_TFO_abj1),
            Adj1_TFO_max=max(DBR_TFO_abj1))

ATAC_cat_summary <- regions_cat %>%
  distinct(X4,.keep_all = TRUE) %>%
  group_by(Categories) %>%
  summarise(Included_ATAC_Regions=n())

ATAC_region_cat_summary <- categories %>%
  select(-ATAC_Region) %>%
  group_by(Categories,DAR_SB,DAR_GH_response,Gene_Waxman_Category,Gene_Cat_Description,ATAC_Waxman_Category,ATAC_Cat_Description) %>%
  summarize(All_ATAC_Regions=n()) %>%
  left_join(ATAC_cat_summary) %>%
  left_join(dbr_cat_summary)

comparison_list_target <- list(list(1:7),
                               list(8:15),
                               list(c(1:7,8:15)),
                               list(1),
                               list(8))
comparison_list_background <- list(list(15:17),
                                   list(15:17),
                                   list(15:17),
                                   list(2),
                                   list(9))
comparison_list_name <- c("M_specific (1-7) vs Background (15-17) ",
                          "F_specific (8-15) vs Background (15-17) ",
                          "Sex-bias (1-7,8-15) vs Background (15-17) ",
                          "M_specific, GH_repressed, M_specific Genes (1) vs M_specific, GH_repressed, non SB genes (2)",
                          "F_specific, GH_induced, F_specific Genes (8) vs F_specific, GH_induced, non SB genes (9)")

DBD_summary <- data.frame(LncRNA=NULL,DBD=NULL,Comparison=NULL,
                          Target_DBR_span_mean=NULL,Target_DBR_span_max=NULL,
                          Target_DBR_TFO_count_mean=NULL,Target_DBR_TFO_count_max=NULL,
                          Target_DBR_Adj1_TFO_mean=NULL,Target_DBR_Adj1_TFO_max=NULL,
                          Target_DBR_Triplex_Count=NULL,Target_DBR_Total=NULL,Background_DBR_span_mean=NULL,
                          Background_DBR_span_max=NULL,Background_DBR_TFO_count_mean=NULL,Background_DBR_TFO_count_max=NULL,
                          Background_DBR_Adj1_TFO_mean=NULL,Background_DBR_Adj1_TFO_max=NULL,Background_DBR_Triplex_Count=NULL,
                          Background_DBR_Total=NULL,
                          Adj_TFO_Ratio = NULL,
                          Adj_TFO_Ratio_Scaled = NULL,
                          Target_Triple_Ratio = NULL,
                          Odd_Ratio=NULL,Fisher_Test_Pval=NULL)

for(i in 1:length(unique(dbr_cat$lnc))){
  print(paste0('LncRNA ',i,' of ',length(unique(dbr_cat$lnc))))
  unique_DBD <- unique(dbr_cat[dbr_cat$lnc == unique(dbr_cat$lnc)[i],]$DBD)
  for(k in unique_DBD){
    print(paste0('DBD ',k,' of ',length(unique_DBD)))
    for(j in 1:length(comparison_list_target)){
      target <- dbr_cat[dbr_cat$lnc == unique(dbr_cat$lnc)[i] &
                          dbr_cat$DBD == k &
                          dbr_cat$Categories %in% unlist(comparison_list_target[[j]]),]
      target_total <- sum(ATAC_region_cat_summary[ATAC_region_cat_summary$Categories %in% unlist(comparison_list_target[[j]]),'Included_ATAC_Regions'])
      background <- dbr_cat[dbr_cat$lnc == unique(dbr_cat$lnc)[i] &
                              dbr_cat$DBD == k &
                              dbr_cat$Categories %in% unlist(comparison_list_background[[j]]),]
      background_total <- sum(ATAC_region_cat_summary[ATAC_region_cat_summary$Categories %in% unlist(comparison_list_background[[j]]),'Included_ATAC_Regions'])
      
      if(nrow(target) != 0){
        t_sum <- target %>%
          group_by(lnc) %>%
          summarise(Target_DBR_span_mean=mean(DBR_span),
                    Target_DBR_span_max=max(DBR_span),
                    Target_DBR_TFO_count_mean = mean(DBR_TFO_count),
                    Target_DBR_TFO_count_max = max(DBR_TFO_count),
                    Target_DBR_Adj1_TFO_mean=mean(DBR_TFO_abj1),
                    Target_DBR_Adj1_TFO_max=max(DBR_TFO_abj1),
                    Target_DBR_Triplex_Count=n()) %>%
          mutate(Target_DBR_Total=target_total) %>%
          select(-lnc)
      }else{
        t_sum <- data.frame(Target_DBR_span_mean=0,
                            Target_DBR_span_max=0,
                            Target_DBR_TFO_count_mean = 0,
                            Target_DBR_TFO_count_max = 0,
                            Target_DBR_Adj1_TFO_mean=0,
                            Target_DBR_Adj1_TFO_max=0,
                            Target_DBR_Triplex_Count=0,
                            Target_DBR_Total=0)
      }
      
      if(nrow(background) != 0){
        b_sum <- background %>%
          group_by(lnc) %>%
          summarise(Background_DBR_span_mean=mean(DBR_span),
                    Background_DBR_span_max=max(DBR_span),
                    Background_DBR_TFO_count_mean = mean(DBR_TFO_count),
                    Background_DBR_TFO_count_max = max(DBR_TFO_count),
                    Background_DBR_Adj1_TFO_mean=mean(DBR_TFO_abj1),
                    Background_DBR_Adj1_TFO_max=max(DBR_TFO_abj1),
                    Background_DBR_Triplex_Count=n()) %>%
          mutate(Background_DBR_Total=background_total) %>%
          select(-lnc)
      }else{
        b_sum <- data.frame(Background_DBR_span_mean=0,
                            Background_DBR_span_max=0,
                            Background_DBR_TFO_count_mean = 0,
                            Background_DBR_TFO_count_max = 0,
                            Background_DBR_Adj1_TFO_mean=0,
                            Background_DBR_Adj1_TFO_max=0,
                            Background_DBR_Triplex_Count=0,
                            Background_DBR_Total=0)
      }
      
      if(nrow(target) != 0 & nrow(background) != 0){
        temp <- data.frame(LncRNA=unique(dbr_cat$lnc)[i],DBD=unique(dbr_cat$DBD)[k],Comparison=comparison_list_name[j],t_sum,b_sum) %>%
          mutate(Adj_TFO_Ratio = Target_DBR_Adj1_TFO_mean/Background_DBR_Adj1_TFO_mean,
                 Adj_TFO_Ratio_Scaled = Adj_TFO_Ratio/(Target_DBR_Total/Background_DBR_Total),
                 Target_Triple_Ratio = Target_DBR_Triplex_Count/Target_DBR_Total,
                 Odd_Ratio=(Target_DBR_Triplex_Count/Target_DBR_Total)/(Background_DBR_Triplex_Count/Background_DBR_Total),
                 Fisher_Test_Pval=fisher.test(matrix(c(Target_DBR_Triplex_Count,
                                                       Target_DBR_Total,
                                                       Background_DBR_Triplex_Count,
                                                       Background_DBR_Total),
                                                     nrow = 2),
                                              alternative = 'greater')$p.value)
      }else{
        if(nrow(target) != 0){
          temp <- data.frame(LncRNA=unique(dbr_cat$lnc)[i],
                             DBD=unique(dbr_cat$DBD)[k],
                             Comparison=comparison_list_name[j],t_sum,b_sum) %>%
            mutate(Adj_TFO_Ratio = Inf,
                   Adj_TFO_Ratio_Scaled = Inf,
                   Target_Triple_Ratio = Target_DBR_Triplex_Count/Target_DBR_Total,
                   Odd_Ratio=Inf,
                   Fisher_Test_Pval=NA)
        }else{
          temp <- data.frame(LncRNA=unique(dbr_cat$lnc)[i],
                             DBD=unique(dbr_cat$DBD)[k],
                             Comparison=comparison_list_name[j],t_sum,b_sum) %>%
            mutate(Adj_TFO_Ratio = 0,
                   Adj_TFO_Ratio_Scaled = 0,
                   Target_Triple_Ratio = 0,
                   Odd_Ratio=0,
                   Fisher_Test_Pval=NA)
        }
      }
      DBD_summary <- rbind(DBD_summary,temp)
    }
  }
}

write_tsv(DBD_summary,file = paste0(dir,"output/",project,"/04/LncRNA_ComparisonSummary.txt"))