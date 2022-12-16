##### Functions #####

## Reassign - Removed certain select column names
re_assign <- function(df, drop) {
  df <- df [, ! names(df) %in% drop, drop = FALSE]
  df
}

## Calculates non-overlapping domains based on a vector of spans
find_domain <- function(l,M){
  counter <- 1
  r <- l
  for(i in seq_along(l)){
    if (l[i] <= M) {
      r[i] <- counter
    } else {
      counter <- counter+1
      r[i] <- counter
    }
  }
  r
}

find_domain2 <- function(l,chr){
  counter <- 1
  r <- l
  for(i in seq_along(l)){
    if (l[i] < 0 & chr[i] == 1) {
      r[i] <- counter
    } else {
      counter <- counter+1
      r[i] <- counter
    }
  }
  r
}

## Calculates the number of unique nucleotides from a series of a vector of start and stop sites
span_cal <- function(x,y){
  z <- cbind(x,y)
  s <- NULL 
  for(i in 1:nrow(z)){
    a <- seq(z[i,1],z[i,2],by=1)
    s <- unique(c(s,a))
  }
  return(length(unique(s))-1)
}

grep_large <- function(pattern,values){
  if(length(pattern) <= 1000){
    return(grep(paste(pattern,collapse = "|"),values,value=T))
  }else{
    start <- seq(from=1,to=length(pattern),by=1000)
    regions <- NULL
    for(i in start){
      if(i < max(start)){
        regions <- c(regions,grep(paste(pattern[i:(i+999)],collapse = "|"),values,value=T))
      }else{
        regions <- c(regions,grep(paste(pattern[i:length(pattern)],collapse = "|"),values,value=T))
      }
    }
    return(regions)
  }
}

### Select specific value from character delimited string

return_value <- function(x,split="_",num_cols=1,col_return=1){
  y <- matrix(unlist(strsplit(x,split = split)),ncol = num_cols,byrow = T)
  return(y[,col_return])
}

#### Vizualization Functions #####

### Function generates simply visualization for DBR data using TFO adbj1 1kb and gene coverage as x and y axis

# Need ggplot2, cowplot, and ggpubr libraries

# dat = merge dataframe of target and background dbrs
# Num_cand = number of candidates
# min_rank = minimum rank value (base on DBR_TFO_abj1_sum_1kbstandardized_rank) for EITHER target or background regions
# title = title for plot
#temp <- dat_merge[which(dat_merge$sum_rank_geneProp_TFO < label_rank & dat_merge$Background == F),]

DBREnrichment_viz <- function(dat,Num_Cand,min_rank,label_rank,title="NULL"){
  
  #dat <- dat_merge
  dat$Background_factor <- as.factor(ifelse(dat$Background == T,"0","1"))
  
  p1 <- ggplot(dat[dat$sum_rank_geneProp_TFO <= min_rank,],
               aes(x=gene_prop*100,
                   y=DBR_TFO_abj1_sum_1kbstandardized,
                   colour=Background_factor,
                   alpha=Background_factor)) +
    geom_point() +
    geom_label_repel(data=dat[which(dat$sum_rank_geneProp_TFO <= label_rank & dat$Background_factor == 1),],aes(label = paste0(lnc,"\n",gname)),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     size=3,
                     segment.color = 'black') +
    scale_alpha_discrete(range=c(0.4,0.6),labels=c("Background","Target")) +
    scale_colour_manual(labels=c("Background","Target"),values=c("darkgrey","darkorange2")) +
    labs(x="Promoter Coverage (%)",
         y="Number of TFOs (abj1 - 1kb)",
         colour="LncRNA-gene\nDBRs",
         alpha="LncRNA-gene\nDBRs") +
    theme_cowplot() +
    theme(legend.title = element_text(hjust=0.5),
          plot.title = element_blank(),
          plot.margin=margin(l=0.1,b=0.1,t=-1,r=-1,unit="cm")) +
    border() 
  p1
  p2 <- p1 + rremove("legend")
  #p3 <- get_legend(p1)
  
  xplot <- ggplot(dat[dat$sum_rank_geneProp_TFO <= min_rank,],
                  aes(x=gene_prop*100,
                      fill = Background_factor,
                      alpha=Background_factor)) +
    geom_density() +
    scale_fill_manual(labels=c("Background","Target"),values=c("darkgrey","darkorange2")) +
    scale_alpha_discrete(range=c(0.2,0.4),labels=c("Background","Target")) +
    clean_theme() + 
    theme_cowplot() +
    labs(x="Promoter Coverage (%)",
         fill="LncRNA-gene\nDBRs",
         alpha="LncRNA-gene\nDBRs") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_text(hjust=0.5),
          plot.margin=margin(t=0.1,l=0.1,b=-1,r=0.1,unit="cm")) 
  
  xplot_legend <- get_legend(xplot)
  xplot <- xplot + rremove("legend")
  
  yplot <- ggplot(dat[dat$sum_rank_geneProp_TFO <= min_rank,],
                  aes(x=DBR_TFO_abj1_sum_1kbstandardized,
                      fill = Background_factor,
                      alpha=Background_factor)) +
    geom_density() +
    scale_fill_manual(labels=c("Background","Target"),values=c("darkgrey","darkorange2")) +
    scale_alpha_discrete(range=c(0.1,0.4),labels=c("Background","Target")) +
    rotate() +
    clean_theme() + 
    theme_cowplot() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin=margin(t=0.1,r=0.1,b=0.1,l=-1,unit="cm")) +
    rremove("legend")
  
  p5 <- plot_grid(xplot,xplot_legend,p2,yplot, ncol = 2, 
                  align = c("hv","none","hv","hv"), 
                  axis = c("tblr","none","tblr","tblr"),
                  rel_widths = c(3, 1), rel_heights = c(1, 3))
  title <- ggdraw() + draw_label(paste0(title," : Top ",Num_Cand," candidates"), fontface='bold')
  plot_grid(title, p5, ncol=1, rel_heights=c(0.1, 1))
}

DBREnrichment_viz2 <- function(dat,
                               tfo_thres,tfo_per,
                               tfo_adjusted=F,
                               gene_prop_thres,
                               gene_prop_per,
                               label_rank,
                               num_lnc = 10,
                               title="NULL",
                               save=F){
  
  dat$Background_factor <- as.factor(ifelse(dat$Background == T,"0","1"))
  if(tfo_adjusted == T){
    dat_temp <- dat[dat$DBR_TFO_abj1_sum_1kbstandardized >= tfo_thres & dat$gene_prop >= gene_prop_thres,]
    dat_temp$TFO_count <- dat_temp$DBR_TFO_abj1_sum_1kbstandardized
    dat <- dat %>%
      arrange(sum_rank_geneProp_TFO)
    dat$TFO_count <- dat$DBR_TFO_abj1_sum_1kbstandardized
    dat_labels <- match(unique(dat$lnc),dat$lnc)[1:num_lnc]
    #dat_labels <- which(dat$sum_rank_geneProp_TFO <= label_rank & dat$Background_factor == 1)
  }else{
    dat_temp <- dat[dat$DBR_TFO_abj1_sum >= tfo_thres & dat$gene_prop >= gene_prop_thres,]
    dat_temp$TFO_count <- dat_temp$DBR_TFO_abj1_sum
    dat <- dat %>%
      arrange(sum_rank_geneProp_TFOsum)
    dat$TFO_count <- dat$DBR_TFO_abj1_sum
    #dat_labels <- which(dat$sum_rank_geneProp_TFOsum <= label_rank & dat$Background_factor == 1)
    dat_labels <- match(unique(dat$lnc), dat$lnc)[1:num_lnc]
  }
  
  bg <- sum(dat_temp$Background == T)
  tar <- sum(dat_temp$Background == F)
  p1 <- ggplot(dat_temp,
               aes(x=gene_prop*100,
                   y=TFO_count,
                   colour=Background_factor,
                   alpha=Background_factor)) +
    geom_hline(aes(yintercept = tfo_thres)) +
    geom_vline(aes(xintercept = gene_prop_thres*100)) +
    geom_point() +
    annotate("text",x=gene_prop_thres*100+1,y = max(dat_temp$TFO_count),
             hjust=0,vjust=1,
             label = paste0("Targets: ",tar,"\n",
                            "Background: ",bg,"\n",
                            "Total: ",nrow(dat_temp),"\n",
                            "Top DBR for first ",num_lnc," lncRNAs")) +
    geom_label_repel(data=dat[dat_labels,],aes(label = paste0(lnc,"\n",gname),colour=Background_factor),
                     box.padding   = 0.25,
                     alpha = 0.75,
                     max.overlaps = 15,
                     min.segment.length = unit(0, 'lines'),
                     point.padding = 0.1,
                     size=3,
                     segment.color = 'black') +
    scale_alpha_discrete(range=c(0.4,0.6),labels=c("Background","Target")) +
    scale_colour_manual(labels=c("Background Genes","Sex-Bias Genes (Target)"),values=c("darkgrey","darkorange2")) +
    labs(x="Promoter Region Coverage (%)",
         y="TFO Count (Adjusted - 1kb)",
         colour="LncRNA-gene\nDBRs",
         alpha="LncRNA-gene\nDBRs") +
    theme_cowplot() +
    theme(legend.title = element_text(hjust=0.5),
          plot.title = element_blank(),
          plot.margin=margin(l=0.1,b=0.1,t=0,r=-1,unit="cm")) +
    border()
  p1
  p2 <- p1 + rremove("legend")
  #p3 <- get_legend(p1)
  
  xplot <- ggplot(dat_temp,
                  aes(x=gene_prop*100,
                      fill = Background_factor,
                      alpha = Background_factor)) +
    geom_density() +
    scale_fill_manual(labels=c("Background Genes","Sex-Bias Genes (Target)"),values=c("darkgrey","darkorange2")) +
    scale_alpha_discrete(range=c(0.2,0.4),labels=c("Background Genes","Sex-Bias Genes (Target)")) +
    clean_theme() + 
    theme_cowplot() +
    labs(x="Promoter Coverage (%)",
         fill="LncRNA-gene\nDBRs",
         alpha="LncRNA-gene\nDBRs") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_text(hjust=0.5),
          plot.margin=margin(t=0.1,l=0.1,b=-1,r=0.1,unit="cm")) 
  
  xplot_legend <- get_legend(xplot)
  xplot <- xplot + rremove("legend")
  
  yplot <- ggplot(dat_temp,
                  aes(x = TFO_count,
                      fill = Background_factor,
                      alpha = Background_factor)) +
    geom_density() +
    scale_fill_manual(labels=c("Background","Target"),values=c("darkgrey","darkorange2")) +
    scale_alpha_discrete(range=c(0.1,0.4),labels=c("Background","Target")) +
    rotate() +
    clean_theme() + 
    theme_cowplot() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin=margin(t=0.1,r=0.1,b=0.1,l=-1,unit="cm")) +
    rremove("legend")
  
  p5 <- plot_grid(xplot,xplot_legend,p2,yplot, ncol = 2, 
                  align = c("hv","none","hv","hv"), 
                  axis = c("tblr","none","tblr","tblr"),
                  rel_widths = c(3, 1), rel_heights = c(1, 3))
  
  if(tfo_adjusted == T){
    title <- ggdraw() + draw_label(paste0(title," : Scaled TFO (1kb std.) thres: ",round(tfo_thres,3)," - Gene Prop.: ",round(gene_prop_thres,3)), fontface='bold')
  }else{
    title <- ggdraw() + draw_label(paste0(title," : TFO thres ",round(tfo_thres,3)," - Gene Prop. ",round(gene_prop_thres,3)), fontface='bold')
  }
  
  p6 <- plot_grid(title, p5, ncol=1, rel_heights=c(0.1, 1))
  
  if(save == T){
    if(tfo_adjusted == T){
      ggsave(p6,filename = paste0(dir,"output/",project,"/sbDBRsummary_TopCandidateViz_TFOScaled_TFO_",tfo_per,"_geneProp_",gene_prop_per,"_Figure.png"),
            width=18,height = 7,dpi = 600)
    }else{
      ggsave(p6,filename = paste0(dir,"output/",project,"/sbDBRsummary_TopCandidateViz_TFONotScaled_TFO_",tfo_per,"_geneProp_",gene_prop_per,"_Figure.png"),
             width=18,height = 7,dpi = 600)
    }
  }
  return(p6)
}


#### DBR Threshold estimate figure ####

DBRThres_viz <- function(dat,criteria,rank,probs=c(0.01,0.05,0.10),sample_num=NULL,
                         criteria_label=NULL,rank_label=NULL,title_label=NULL){
  #rank = 'gene_prop_rank'
  #criteria = 'gene_prop'
  #Subsample all data to reduce plotting time
  if(is.null(sample_num)){dat_sample <- dat}else{
    dat_sample <- dat[sample(1:nrow(dat),sample_num,replace = F),]
  }
  dat_sample$x <- unlist(dat_sample[,rank])
  dat_sample$y <- unlist(dat_sample[,criteria])
  
  cols <- colorRampPalette(c("red3","steelblue2"))(length(probs))
  # Calculate quantile value for various probabilities
  thres <- NULL
  for(i in 1:length(probs)){
    thres <- c(thres,quantile(unlist(dat[,criteria]),probs = 1-probs[i]))
  }
  offset <- quantile(unlist(dat[,criteria]),probs = 0.01)
  
  p1 <- ggplot(dat_sample,aes(x=x,y=y)) 
  for(i in 1:length(probs)){
    p1 <- p1 + 
      geom_segment(x = 0, y = thres[i], xend = max(dat_sample$x)*0.50, yend = thres[i],
                   colour=cols[i],size=2,lineend = c('round')) +
      #geom_hline(yintercept=thres[i],colour=cols[i],size=2) +
      geom_text(label=paste0("Probability: ",probs[i]*100,"% - Thres. value: ",round(as.numeric(thres[i]),3)),
                x=max(dat_sample$x)*0.75,
                y=thres[i],
                check_overlap = T)
  }
  p1 <- p1 + geom_point() +
    labs(y=criteria_label,
         x=rank_label,
         title=title_label) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(thresholds=thres,
              plot=p1))
}

sbDBRsummary_viz <- function(dat_target,
                             dat_background,
                             tfo_thres,tfo_adjusted=F,
                             gene_prop_thres,
                             tfo_per,gene_prop_per,
                             save=F){
  
  ## Summarize lncRNA for top DBR candidates with SB genes
  if(tfo_adjusted == T){
    dat_target_candidates <- dat_target[dat_target$DBR_TFO_abj1_sum_1kbstandardized >= tfo_thres & dat_target$gene_prop >= gene_prop_thres,]
  }else{
    dat_target_candidates <- dat_target[dat_target$DBR_TFO_abj1_sum >= tfo_thres & dat_target$gene_prop >= gene_prop_thres,]
  }
  
  dat_target_GeneSummary <- dat_target_candidates %>%
    group_by(lnc,SB_Response_Summary,gname) %>%
    summarise(freq = n())
  unique(dat_target_GeneSummary$gname)
  top_lncs <- data.frame(table(dat_target_candidates$lnc,
                               dat_target_candidates$SB_Response_Summary))
  
  colnames(top_lncs) <- c("lnc","Promoter_SB","count")
  top_lncs <- left_join(top_lncs,lnclist,by=c("lnc"= "Lnc_Id"))
  
  ## Summarize lncRNA for top DBR candidates with background genes
  max_rank <- max(dat_target_candidates$sum_rank_geneProp_TFOsum) # max rank of candidate DBRs w/ SB genes
  top_background <- data.frame(table(dat_background$lnc[dat_background$sum_rank_geneProp_TFOsum <= max_rank]))
  colnames(top_background) <- c("lnc","count_b")
  top_lnc_alt <- top_lncs %>%
    group_by(lnc) %>%
    summarize(gene_f=count[Promoter_SB == "SexBiased_Female"],
              gene_m=count[Promoter_SB == "SexBiased_Male"])
  top_background_match <- left_join(top_lnc_alt,top_background)
  
  sum(top_background_match$gene_f) + sum(top_background_match$gene_m)
  sum(top_background_match$count_b)
  
  ### Panel 1 -  Count number of DBRs between lncRNAs and SB genes by sex
  # Order factor based on number of genes
  top_lncs$lnc <- as.factor(top_lncs$lnc)
  top_lncs$lnc <- factor(top_lncs$lnc,levels = unique(top_lncs$lnc[order(top_lncs$count,decreasing = T)]))
  top_lncs_alt <- top_lncs
  top_lncs_alt$Hypoxclass[is.na(top_lncs_alt$Hypoxclass)] <- "No Class"
  top_lncs_alt$Hypoxclass[top_lncs_alt$Hypoxclass == unique(top_lncs_alt$Hypoxclass)[6]] <-  "Female Class I/II"
  levels(top_lncs_alt$Promoter_SB) <- c("Gene - Female Bias","Gene - Male Bias")
  
  
  p1 <- ggplot(top_lncs_alt,aes(y=count,x=lnc,fill=Hypoxclass,colour=Sex_Bias)) +
    facet_grid(cols = vars(Promoter_SB)) +
    geom_col(width=0.8,size=0.8) +
    geom_label(aes(y = -18,x = lnc,label=lnc),colour="grey90", hjust = 0) +
    theme_cowplot() +
    scale_fill_manual(values = c("palevioletred1","palevioletred3","palevioletred4",
                                 "steelblue1","steelblue4",
                                 "grey25")) +
    scale_colour_manual(values = c("red3","navyblue"),
                        guide = guide_legend(override.aes = list(linetype = c(0,0),
                                                                 fill = c("red3","navyblue"),
                                                                 shape = c(2,2),
                                                                 colour = c("red3","navyblue")))) +
    labs(y="Number of Genes",
         x="",
         fill="LncRNA\nHypox Class",
         colour="LncRNA\nSex Bias") +
    theme(legend.title = element_text(hjust=0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_line(),
          strip.background = element_blank()) +
    coord_cartesian(ylim = c(0, 100),xlim = c(0, 100),clip = "off") +
    rotate()
  p1
  #p1_nolegend <- p1 + rremove("legend")
  #p1_legend <- get_legend(p1)
  
  # p1 <- ggplot(top_lncs_alt,aes(y=count,x=lnc,fill=Hypoxclass,colour=Sex_Bias)) +
  #   facet_grid(cols = vars(Promoter_SB)) +
  #   geom_col(width=0.8,size=0.8) +
  #   theme_cowplot() +
  #   scale_fill_manual(values = c("palevioletred1","palevioletred3","palevioletred4",
  #                                "steelblue1","steelblue4",
  #                                "grey50")) +
  #   scale_colour_manual(values = c("red3","navyblue"),
  #                       guide = guide_legend(override.aes = list(linetype = c(0,0),
  #                                                                fill = c("red3","navyblue"),
  #                                                                shape = c(2,2),
  #                                                                colour = c("red3","navyblue")))) +
  #   labs(y="Number of Genes",
  #        x="",
  #        fill="LncRNA\nHypox Class",
  #        colour="LncRNA\nSex Bias") +
  #   theme(legend.title = element_text(hjust=0.5),
  #         panel.grid.major.y = element_line(),
  #         strip.background = element_blank()) +
  #   rotate()
  # p1
  ### Calculate proportion of DBRS between lncRNAs and either malebiased,female biased, or backgrouund genes
  top_background_match$count_b[is.na(top_background_match$count_b)] <- 0
  top_background_match$count_sum <- top_background_match$gene_f+top_background_match$gene_m+top_background_match$count_b
  top_background_match$per_f <- top_background_match$gene_f/top_background_match$count_sum
  top_background_match$per_m <- top_background_match$gene_m/top_background_match$count_sum
  top_background_match$per_b <- top_background_match$count_b/top_background_match$count_sum
  top_background_match$gene_f_scale <- top_background_match$gene_f/length(glist$Gene[glist$SB_Response_Summary == "SexBiased_Female"])
  top_background_match$gene_m_scale <- top_background_match$gene_m/length(glist$Gene[glist$SB_Response_Summary == "SexBiased_Male"])
  top_background_match$gene_b_scale <- top_background_match$count_b/length(glist$Gene[glist$Background == T])
  top_background_match$scale_sum <- top_background_match$gene_f_scale + top_background_match$gene_m_scale + top_background_match$gene_b_scale
  top_background_match$scale_f <- top_background_match$gene_f_scale/top_background_match$scale_sum
  top_background_match$scale_m <- top_background_match$gene_m_scale/top_background_match$scale_sum
  top_background_match$scale_b <- top_background_match$gene_b_scale/top_background_match$scale_sum
  
  top_background_match_rv <- rbind(data.frame(lnc=top_background_match$lnc,
                                              bias="Female",
                                              per=top_background_match$per_f,
                                              scale=top_background_match$scale_f),
                                   data.frame(lnc=top_background_match$lnc,
                                              bias="Male",
                                              per=top_background_match$per_m,
                                              scale=top_background_match$scale_m),
                                   data.frame(lnc=top_background_match$lnc,
                                              bias="Background",
                                              per=top_background_match$per_b,
                                              scale=top_background_match$scale_b))
  top_background_match_rv$lnc <- as.factor(top_background_match_rv$lnc)
  top_background_match_rv$lnc <- factor(top_background_match_rv$lnc,levels = unique(top_lncs$lnc[order(top_lncs$count,decreasing = T)]))
  p_col <- viridisLite::viridis(2)
  # Proportion of each gene list unscaled
  p2 <- ggplot(top_background_match_rv,aes(x=lnc,y=per,fill=bias)) +
    geom_bar(position="stack", stat="identity",size=0.8,alpha=0.8) +
    #scale_fill_manual(values = c("grey80","palevioletred3","steelblue3")) +
    scale_fill_manual(values = c("grey80",p_col)) +
    theme_cowplot() +
    labs(y="Proportion",
         fill="Gene\nSex Bias") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.title = element_text(hjust=0.5),
          panel.grid.major.y = element_line()) +
    rotate()
  p2_nolegend <- p2 + rremove('legend')
  p2_legend <- get_legend(p2)
  
  # Proportion of each gene list scaled
  p3 <- ggplot(top_background_match_rv,aes(x=lnc,y=scale,fill=bias)) +
    geom_bar(position="stack", stat="identity",size=0.8,alpha=0.8) +
    scale_fill_manual(values = c("grey80",p_col)) +
    theme_cowplot() +
    labs(y="Proportion (Scaled)",
         fill="Gene\nSex Bias") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          panel.grid.major.y = element_line()) +
    rotate()
  
  p4 <- plot_grid(p1_legend,p2_legend,nrow=2)
  p5 <- plot_grid(p1_nolegend,p2_nolegend,p3,p4,ncol = 4,align = "hv",axis = c("tb","tb","tb"),rel_widths = c(5,2,2,1))
  title <- ggdraw() + draw_label(paste0("Top Candidate sex-biased DBR summary - TFO: ",round(tfo_thres,2)," Gene Prop: ",round(gene_prop_thres,2)), fontface='bold')
  p6 <- plot_grid(title, p5, ncol=1, rel_heights=c(0.1, 4))
  if(save == T){
    ggsave(p6,filename = paste0(dir,"output/",project,"/sbDBRsummary_TFO_",tfo_per,"_geneProp_",gene_prop_per,"_Figure.png"),
           width=22,height=9,dpi = 600)
  }
  return(p6)
}
