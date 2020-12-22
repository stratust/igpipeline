library(dplyr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(writexl)
library(argparser, quietly = TRUE)
library(ggplot2)
library(svglite)
library(reshape2)
library(scales)
library(ggpubr)

get.unified.df = function ( ig.df ) {
  
  ig.df.na.rows.removed = ig.df[rowSums(is.na(ig.df)) != ncol(ig.df), ]
  
  seqid.ref.column = grep('SEQUENCE_ID',names(ig.df),value=F)
  names(seqid.ref.column) = c('IGA','IGM','IGG','VK','VL')
  HC.ref.column = seqid.ref.column[1:3]
  LC.ref.column = seqid.ref.column[4:5]
  
  # Unify heavy chain
  HC.df = NULL
  isotype <- NULL 
  for (i in names(HC.ref.column)) {
    ref.column = HC.ref.column[i]
    seqid.vector = ig.df.na.rows.removed[,ref.column]
    column.range =  ref.column:(ref.column + 11)
    
    if (is.null(HC.df) & sum(!is.na(seqid.vector)) > 0){
      HC.df = ig.df.na.rows.removed[,column.range]
      isotype <- rep(i,nrow(HC.df))
    } else if (!is.null(HC.df) & sum(!is.na(seqid.vector)) > 0) {
      selected.rows = which(!is.na(seqid.vector))
      HC.df[selected.rows,] =  ig.df.na.rows.removed[selected.rows,column.range]
      isotype[selected.rows] = i
    }
    
  }
  
  # Unify light chain
  LC.df = NULL
  light.chain = NULL
  for (i in names(LC.ref.column)) {
    ref.column = LC.ref.column[i]
    seqid.vector = ig.df.na.rows.removed[,ref.column]
    column.range =  ref.column:(ref.column + 11)
    
    if (is.null(LC.df) & sum(!is.na(seqid.vector)) > 0){
      LC.df = ig.df.na.rows.removed[,column.range]
      light.chain = rep(i, nrow(LC.df))
    } else if (!is.null(LC.df) & sum(!is.na(seqid.vector)) > 0) {
      selected.rows = which(!is.na(seqid.vector))
      LC.df[selected.rows,] =  ig.df.na.rows.removed[selected.rows,column.range]
      light.chain[selected.rows] = i
    }
    
  }
  
  unified.df = cbind(data.frame(ISOTYPE=isotype,LIGHTCHAIN=light.chain),ig.df.na.rows.removed[,1:3],HC.df,LC.df)
  unified.df = unified.df %>% mutate( CLUSTERCLASS = ifelse( num_seqs > 1, 'clone', 'singlet' ) )
  
  return(unified.df)
  
}

parse_excel <- function(file) {
  xl.obj = as.data.frame( read_excel( file, sheet = "PROPER", skip = 1, col_names = T ) )
  
  unified.df = get.unified.df( ig.df = xl.obj )
  
  return(unified.df)
}

get_v_gene_plot <- function(df) {
  
  order.df <- df %>% group_by(V_CALL) %>% summarise(max = max(frequency), min = min(frequency)) %>% arrange(max, min, desc(V_CALL))
  df$V_CALL <- factor(df$V_CALL, levels <- order.df$V_CALL)
  
  calculate.binomial.test <- function(data) {

    t1.row <- data %>% filter(TIMEPOINT == 'Baseline')
    t2.row <- data %>% filter(TIMEPOINT == '6mo')
    output.df <- data.frame(group1 = 'Baseline', group2 = '6mo')
    
    output.df$frequency <- max(data$frequency)
    
    n.t1 <- 0
    total.t1 <- 0 
    prob <- 0
    if (nrow(t1.row) > 0){
      n.t1 <- t1.row$n
      total.t1 <- t1.row$total
      prob = n.t1/total.t1
    }
    
    n.t2 <- 0
    total.t2 <- 0 
    if (nrow(t2.row) > 0 & nrow(t1.row) > 0 ){
      n.t2 <- t2.row$n
      total.t2 <- t2.row$total
      result = binom.test(n.t2, total.t2, p = prob)
      output.df$p = result$p.value
    } else{
      output.df$p <- NA
    }
    
    return(output.df)
  }
  
  pvalue.df <-  df %>% group_by(TIMEPOINT) %>% mutate(total=sum(n)) %>% group_by(V_CALL) %>% group_modify(~ calculate.binomial.test(.x) )
  pvalue.df$p.adj <- p.adjust(pvalue.df$p, method = 'fdr')
  pvalue.df <- pvalue.df %>% mutate(
    p.formatted = case_when(
      p.adj < 0.05 ~ '*',
      p.adj < 0.01 ~ '**',
      p.adj < 0.001 ~ '***'
      #is.na(p.adj) ~ NA
    )
  )
  
  xmax <- ceiling(max(df$frequency) * 100)/100
  
  # get title
  plot.title <- 'All patients'

  
  ggplot(df, aes(y = V_CALL, x = frequency)) + 
    geom_bar(aes(fill = TIMEPOINT), position = position_dodge(preserve = "single"), stat="identity", color = 'black', width=0.8) + 
    scale_x_continuous(labels = percent, position = "top") +
    coord_cartesian(expand = FALSE, xlim = c(0,xmax)) +
    theme_classic() + 
    ggtitle(plot.title) +
    theme(
      panel.grid.major.x = element_line(colour="grey", size=0.5, linetype= 'dashed'),
      axis.text = element_text(color = 'black'),
      plot.title = element_text(hjust = 0.5)
    ) + #coord_flip(expand = FALSE) 
    geom_text(data=pvalue.df, aes(label=p.formatted), nudge_x = 0.005)
  
}


get_v_gene_plot_by_patient <- function(df) {
  
  #ggplot(df, aes(y = V_CALL, x = n, fill = ISOTYPE)) + geom_bar(position="dodge", stat="identity") + facet_wrap(sampleid ~ ISOTYPE, scales = 'free')
  
  order.df <- df %>% group_by(V_CALL) %>% summarise(max = max(frequency), min = min(frequency)) %>% arrange(max, min, desc(V_CALL))
  df$V_CALL <- factor(df$V_CALL, levels <- order.df$V_CALL)
  
  calculate.binomial.test <- function(data) {
    t1.row <- data %>% filter(TIMEPOINT == 'Baseline')
    t2.row <- data %>% filter(TIMEPOINT == '6mo')
    output.df <- data %>% dplyr::select(sampleid) %>% distinct() %>% mutate(group1 = 'Baseline', group2 = '6mo')
    
    output.df$frequency <- max(data$frequency)
    
    n.t1 <- 0
    total.t1 <- 0 
    prob <- 0
    if (nrow(t1.row) > 0){
        n.t1 <- t1.row$n
        total.t1 <- t1.row$total
        prob = n.t1/total.t1
    }
    
    n.t2 <- 0
    total.t2 <- 0 
    if (nrow(t2.row) > 0 & nrow(t1.row) > 0 ){
      n.t2 <- t2.row$n
      total.t2 <- t2.row$total
      result = binom.test(n.t2, total.t2, p = prob)
      output.df$p = result$p.value
    } else{
      output.df$p <- NA
    }

    return(output.df)
  }
  
  pvalue.df <-  df %>% group_by(TIMEPOINT) %>% mutate(total=sum(n)) %>% group_by(V_CALL) %>% group_modify(~ calculate.binomial.test(.x) )
  pvalue.df$p.adj <- p.adjust(pvalue.df$p, method = 'fdr')
  pvalue.df <- pvalue.df %>% mutate(
    p.formatted = case_when(
      p.adj < 0.05 ~ '*',
      p.adj < 0.01 ~ '**',
      p.adj < 0.001 ~ '***'
      #is.na(p.adj) ~ NA
    )
  )
  
  xmax <- ceiling(max(df$frequency) * 100)/100
  
  # get title
  samples <- unique(df$sampleid)
  plot.title <- NULL
  if (length(samples) > 1){
    plot.title <- paste0(samples, collapse = ', ')  
  } else {
    plot.title <- samples
  }
  
  
  ggplot(df, aes(y = V_CALL, x = frequency)) + 
    geom_bar(aes(fill = TIMEPOINT), position = position_dodge(preserve = "single"), stat="identity", color = 'black', width=0.8) + 
    scale_x_continuous(labels = percent, position = "top") +
    coord_cartesian(expand = FALSE, xlim = c(0,xmax)) +
    theme_classic() + 
    ggtitle(plot.title) +
    theme(
      panel.grid.major.x = element_line(colour="grey", size=0.5, linetype= 'dashed'),
      axis.text = element_text(color = 'black'),
      plot.title = element_text(hjust = 0.5)
    ) + #coord_flip(expand = FALSE) 
    geom_text(data=pvalue.df, aes(label=p.formatted), nudge_x = 0.005)

}

v_gene_usage <- function(file, output_prefix) {
  patient.df <- parse_excel(file)
  patient.df <- patient.df %>% 
    mutate( TIMEPOINT = case_when(
        ISOTYPE == 'IGG' ~ 'Baseline',
        ISOTYPE == 'IGM' ~ '6mo'
      )
    )
  

  # freqq by patient
  vh.by.patient.df <- patient.df %>%
    mutate( V_CALL...71  = gsub('\\*.*','', V_CALL...71)) %>%
    group_by(sampleid, TIMEPOINT, V_CALL...71) %>% 
    summarise(n = n()) %>% 
    rename( V_CALL = V_CALL...71) %>%
    mutate(frequency = n/sum(n))

  
  vk.by.patient.df <- patient.df %>% 
    filter( LIGHTCHAIN == 'VK') %>%
    mutate( V_CALL...197  = gsub('\\*.*','', V_CALL...197)) %>%
    group_by(sampleid, TIMEPOINT, V_CALL...197) %>% 
    summarise(n = n()) %>% 
    rename( V_CALL = V_CALL...197) %>%
    mutate(frequency = n/sum(n))
  
  vl.by.patient.df <- patient.df %>% 
    filter( LIGHTCHAIN == 'VL') %>% 
    mutate( V_CALL...197  = gsub('\\*.*','', V_CALL...197)) %>%
    group_by(sampleid, TIMEPOINT, V_CALL...197) %>% 
    summarise(n = n()) %>% 
    rename( V_CALL = V_CALL...197) %>%
    mutate(frequency = n/sum(n))
  
  # frequency by chain
  vh.df <- patient.df %>%
    mutate( V_CALL...71  = gsub('\\*.*','', V_CALL...71)) %>%
    group_by(TIMEPOINT, V_CALL...71) %>% 
    summarise(n = n()) %>% 
    rename( V_CALL = V_CALL...71) %>%
    mutate(frequency = n/sum(n))
  
  vk.df <- patient.df %>% 
    filter( LIGHTCHAIN == 'VK') %>%
    mutate( V_CALL...197  = gsub('\\*.*','', V_CALL...197)) %>%
    group_by( TIMEPOINT, V_CALL...197) %>% 
    summarise(n = n()) %>% 
    rename( V_CALL = V_CALL...197) %>%
    mutate(frequency = n/sum(n))
  
  vl.df <- patient.df %>% 
    filter( LIGHTCHAIN == 'VL') %>% 
    mutate( V_CALL...197  = gsub('\\*.*','', V_CALL...197)) %>%
    group_by(TIMEPOINT, V_CALL...197) %>% 
    summarise(n = n()) %>% 
    rename( V_CALL = V_CALL...197) %>%
    mutate(frequency = n/sum(n))
  
  all.patients.plots <- list()
  all.patients.plots[['VH']] <- get_v_gene_plot(vh.df)
  all.patients.plots[['VL'] ]<- get_v_gene_plot(vl.df)
  all.patients.plots[['VK']] <- get_v_gene_plot(vk.df)
  
  for (chain in names(all.patients.plots)) {
    pdf(file = paste0(output_prefix, '/all_patients_',chain,'_distribution.pdf'))
    print(all.patients.plots[[chain]])
    dev.off()
  }
  
  
  vh.list <- split(vh.by.patient.df, vh.by.patient.df$sampleid)

  vh.plots <- lapply(vh.list,get_v_gene_plot_by_patient) 

  for (patient in names(vh.plots)) {
    pdf(file = paste0(output_prefix, '/', patient, '_VH_distribution.pdf'))
    print(vh.plots[[patient]])
    dev.off()
    #svglite(file = paste0(output_prefix, '/', patient, '_VH_distribution.svg'))
    #print(vh.plots[[patient]])
    #dev.off()
  }
  
  vl.list <- split(vl.by.patient.df, vl.by.patient.df$sampleid)
  
  vl.plots <- lapply(vl.list,get_v_gene_plot_by_patient) 
  
  for (patient in names(vl.plots)) {
    pdf(file = paste0(output_prefix, '/', patient, '_VL_distribution.pdf'))
    print(vl.plots[[patient]])
    dev.off()
    #svglite(file = paste0(output_prefix, '/', patient, '_VL_distribution.svg'))
    #print(vl.plots[[patient]])
    #dev.off()
  }
  
  vk.list <- split(vk.by.patient.df, vk.by.patient.df$sampleid)
  
  vk.plots <- lapply(vk.list, get_v_gene_plot_by_patient) 
  
  for (patient in names(vk.plots)) {
    pdf(file = paste0(output_prefix, '/', patient, '_VK_distribution.pdf'))
    print(vk.plots[[patient]])
    dev.off()
    #svglite(file = paste0(output_prefix, '/', patient, '_VK_distribution.svg'))
    #print(vk.plots[[patient]])
    #dev.off()
  }
  
}

main <- function() {
  # Create a parser
  p <- arg_parser("Parse IgParse.pl output and create V gene usage plots and tables")
  
  # Add command line arguments
  p <- add_argument(p, "--input_file", help = "The excel file")
  p <- add_argument(p, "--output_prefix", help = "output prefix")
  
  #Parse the command line arguments
  argv <- parse_args(p)
  argv$input_file <- list.files(pattern = ".*all.*_clonal.xlsx", path = '/home/stratus/GitHub/igpipeline2_timepoint//results/PATIENTS/igpairs/', full.names = T)
  argv$output_prefix <- '/home/stratus/TMP/'
  
  # check if input and codefile arguments are not NA
  if (is.na(argv$input_file) || is.na(argv$output_prefix)) {
    print(p)
    quit()
  }
  
  
  # For strict files COV
  v_gene_usage(argv$input_file, output_prefix = argv$output_prefix)
}


# Call main function
main()
