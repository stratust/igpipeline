#################################################### Libraries ####################################################
library( readxl )
library( dplyr )
library( ggplot2 )
library( ggpubr )
library( reshape2 )
library( argparser )
library( ggbeeswarm )
library( Peptides )
library( writexl )
library( parallel )
#library( plotly )


p <- arg_parser("SHM analysis")

p <- add_argument(p, "--input", help="directory where the excel files are located")

p <- add_argument(p, "--bcelldbtsv", help="TSV file containing the CDR3 of the B cell database")
p <- add_argument(p, "--rdsbcelldb", help="RDS file containing CDR3 score of the B cell database")

p <- add_argument(p, "--output", help="output to export the images", default="output.txt")

argv <- parse_args(p)

#argv$input = "/home/stratus/GitHub/igpipeline2_tbev/results/PATIENTS/igpairs/"
#argv$rdsbcelldb = "/home/stratus/GitHub/igpipeline2/database/bcelldb.rds"

#argv$output = "/home/stratus/TMP"

#################################################### Function ####################################################
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

get.ggplot.jitter = function( data, x, fill , ylab = '# of mutations' ,legend = FALSE ){
  class = as.character(unique(data[[x]]))
  
  my_comparisons <- NULL
  if (length(class) > 1){
    my_comparisons = combn(class, 2, simplify = F)
  }

  plot = ggplot(data=data, aes_string(x=x, y='value')) +
    geom_quasirandom(aes_string(fill=fill, text='SEQUENCE_ID'), color='black', groupOnX=TRUE, shape = 21, size = 1.5 ) + 
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="red", size = 0.4) +
    coord_cartesian( ylim=c(0, max(data$value) + (length(my_comparisons) * max(data$value) * 0.3  ) ), expand=T) +
    ylab(ylab) +
    xlab("") +
    theme_bw() + 
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5))
  
  if (legend == FALSE){
    plot = plot + theme(legend.position = "none")
  }
  
  if (!is.null(my_comparisons)){
    plot = plot +
      stat_compare_means(comparisons = my_comparisons,
                         label = "p.signif",
                         hide.ns = F) #+ # Add pairwise comparisons p-value
      #stat_compare_means(label.y =  max(data$value) + (length(my_comparisons) * max(data$value) * 0.15) )  # Add global p-value
  }
  
  return(plot)  
}

get.shm.jitter.plot = function ( df.melt ) {
  data.plot = subset( df.melt, variable %in% c("VH mutations (nt)","VL mutations (nt)","VH + VL mutations (nt)") )
  data.plot$variable= as.character(data.plot$variable)
  
  plots <- list()
  
  plots[['isotype_by_sample']] = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE', fill = 'ISOTYPE') + 
    facet_grid(patient ~ variable) +
    scale_fill_manual(values = c('IGA'='black','IGG'='black','IGM'='black'))
  
  plots[['isotype']] = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE', fill = 'ISOTYPE') + 
    facet_grid(. ~ variable) + ggtitle(paste0('SHM (',paste0(unique(df.melt$patient),collapse=", "),')')) +
    scale_fill_manual(values = c('IGA'='black','IGG'='black','IGM'='black'))
  
  
  plots[['cluster_by_sample']] = get.ggplot.jitter(data = data.plot, x = 'CLUSTERCLASS', fill = 'CLUSTERCLASS') + 
    facet_grid(patient ~ variable) +
    scale_fill_manual(values = c('clone'='black','singlet'='white'))
  
  plots[['cluster']]  = get.ggplot.jitter(data = data.plot, x = 'CLUSTERCLASS', fill = 'CLUSTERCLASS') +
    facet_grid(. ~ variable) + ggtitle(paste0('SHM (',paste0(unique(df.melt$patient),collapse=", "),')')) +
    scale_fill_manual(values = c('clone'='black','singlet'='white'))
  
  
  plots[['isotype_cluster_by_sample']]  = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE_CLUSTERCLASS', fill = 'CLUSTERCLASS') +
    facet_grid(patient ~ variable)  +
    scale_fill_manual(values = c('clone'='black','singlet'='white'))
  
  plots[['isotype_cluster']]   = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE_CLUSTERCLASS', fill = 'CLUSTERCLASS') +
    facet_grid(. ~ variable) + ggtitle(paste0('SHM (',paste0(unique(df.melt$patient),collapse=", "),')')) +
    scale_fill_manual(values = c('clone'='black','singlet'='white'))
 
  n.samples = length(unique(df.melt$patient))
  
  plots_size = list(
    'isotype'= c(5,6)
    , 'cluster'= c(5,6)
    , 'isotype_cluster'= c(5,6)    
    , 'isotype_by_sample' = c(2.85 * n.samples, 6)
    , 'cluster_by_sample' = c(2.85 * n.samples, 6)
    , 'isotype_cluster_by_sample' = c(2.85 * n.samples, 6)
  )
  
  # add SHM to lists names
  names(plots) = paste0('SHM_',names(plots))
  names(plots_size) = paste0('SHM_',names(plots_size))
  
  return(list(plots=plots,size=plots_size))

}

get.cdr3.jitter.plot = function ( df.melt ) {

  data.plot = subset( df.melt, variable %in% c("CDR3 HC length (aa)","CDR3 LC length (aa)") )
  data.plot$variable= as.character(data.plot$variable)
  
  plots <- list()
  
  plots[['isotype_by_sample']] = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE', fill = 'ISOTYPE') + 
    facet_grid(patient ~ variable) +
    scale_fill_manual(values = c('IGA'='black','IGG'='black','IGM'='black')) + ylab('CDR3 length (aa)')
  
  plots[['isotype']] = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE', fill = 'ISOTYPE') + 
    facet_grid(. ~ variable) + ggtitle(paste0('CDR3 Length (',paste0(unique(df.melt$patient),collapse=", "),')')) +
    scale_fill_manual(values = c('IGA'='black','IGG'='black','IGM'='black')) + ylab('CDR3 length (aa)')
  
  
  plots[['cluster_by_sample']] = get.ggplot.jitter(data = data.plot, x = 'CLUSTERCLASS', fill = 'CLUSTERCLASS') + 
    facet_grid(patient ~ variable) +
    scale_fill_manual(values = c('clone'='black','singlet'='white')) + ylab('CDR3 length (aa)')
  
  plots[['cluster']]  = get.ggplot.jitter(data = data.plot, x = 'CLUSTERCLASS', fill = 'CLUSTERCLASS') +
    facet_grid(. ~ variable) + ggtitle(paste0('CDR3 Length (',paste0(unique(df.melt$patient),collapse=", "),')')) +
    scale_fill_manual(values = c('clone'='black','singlet'='white')) + ylab('CDR3 aa length')
  
  
  plots[['isotype_cluster_by_sample']] = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE_CLUSTERCLASS', fill = 'CLUSTERCLASS') +
    facet_grid(patient ~ variable)  +
    scale_fill_manual(values = c('clone'='black','singlet'='white')) + ylab('CDR3 length (aa)')
  
  plots[['isotype_cluster']] = get.ggplot.jitter(data = data.plot, x = 'ISOTYPE_CLUSTERCLASS', fill = 'CLUSTERCLASS') +
    facet_grid(. ~ variable) + ggtitle(paste0('CDR3 Length (',paste0(unique(df.melt$patient),collapse=", "),')')) +
    scale_fill_manual(values = c('clone'='black','singlet'='white')) + ylab('CDR3 length (aa)')
 
  n.samples = length(unique(df.melt$patient))
  
  plots_size = list(
    'isotype'= c(5,8.5)
    , 'cluster'= c(5,8.5)
    , 'isotype_cluster'= c(5,8.5)    
    , 'isotype_by_sample' = c(2.85 * n.samples, 8.5)
    , 'cluster_by_sample' = c(2.85 * n.samples, 8.5)
    , 'isotype_cluster_by_sample' = c(2.85 * n.samples, 8.5)
  )
  
  # add SHM to lists names
  names(plots) = paste0('CDR3_',names(plots))
  names(plots_size) = paste0('CDR3_',names(plots_size))
  
  return(list(plots=plots,size=plots_size))
  
}

create.summary.list <- function ( unified.df ){
  summary.list <- list()
  
  seqid.ref.column = grep('SEQUENCE_ID\\.\\.',names(unified.df),value=T)
  nt_mismatch.ref.column = grep('nt_mismatches_V_region\\.\\.',names(unified.df),value=T)
  aa_mismatch.ref.column = grep('v_region_aa_mismatches\\.\\.',names(unified.df),value=T)
  cdr3aalength.ref.column = grep('cdr3_aa_length\\.\\.',names(unified.df),value=T)
  
  # Rlang unevaluated expressions
  summarise.expression <- quos(
    "NUM SEQS" = as.list(n()),
    "MEAN VH mutations (nt)" = mean( !!as.name(nt_mismatch.ref.column[1]) ),
    "MEDIAN VH mutations (nt)" = median( !!as.name(nt_mismatch.ref.column[1]) ),
    "MIN VH mutations (nt)"  = min( !!as.name(nt_mismatch.ref.column[1]) ),
    "MAX VH mutations (nt)"  = max( !!as.name(nt_mismatch.ref.column[1]) ),
    "MEAN VL mutations (nt)" = mean( !!as.name(nt_mismatch.ref.column[2]) ),
    "MEDIAN VL mutations (nt)" = median( !!as.name(nt_mismatch.ref.column[2]) ),
    "MIN VL mutations (nt)"  = min( !!as.name(nt_mismatch.ref.column[2]) ),
    "MAX VL mutations (nt)"  = max( !!as.name(nt_mismatch.ref.column[2]) ),
    "MEAN VH+VL mutations (nt)" = mean( !!as.name(nt_mismatch.ref.column[1]) + !!as.name(nt_mismatch.ref.column[2]) ),
    "MEDIAN VH+VL mutations (nt)" = median( !!as.name(nt_mismatch.ref.column[1]) + !!as.name(nt_mismatch.ref.column[2])  ),
    "MIN VH+VL mutations (nt)"  = min( !!as.name(nt_mismatch.ref.column[1]) + !!as.name(nt_mismatch.ref.column[2])  ),
    "MAX VH+VL mutations (nt)"  = max( !!as.name(nt_mismatch.ref.column[1]) + !!as.name(nt_mismatch.ref.column[2])  ),
    
    "MEAN VH mutations (aa)" = mean( !!as.name(aa_mismatch.ref.column[1]) ),
    "MEDIAN VH mutations (aa)" = median( !!as.name(aa_mismatch.ref.column[1]) ),
    "MIN VH mutations (aa)"  = min( !!as.name(aa_mismatch.ref.column[1]) ),
    "MAX VH mutations (aa)"  = max( !!as.name(aa_mismatch.ref.column[1]) ),
    "MEAN VL mutations (aa)" = mean( !!as.name(aa_mismatch.ref.column[2]) ),
    "MEDIAN VL mutations (aa)" = median( !!as.name(aa_mismatch.ref.column[2]) ),
    "MIN VL mutations (aa)"  = min( !!as.name(aa_mismatch.ref.column[2]) ),
    "MAX VL mutations (aa)"   = max( !!as.name(aa_mismatch.ref.column[2]) ),
    "MEAN VH+VL mutations (aa)" = mean( !!as.name(aa_mismatch.ref.column[1]) + !!as.name(aa_mismatch.ref.column[2]) ),
    "MEDIAN VH+VL mutations (aa)" = median( !!as.name(aa_mismatch.ref.column[1]) + !!as.name(aa_mismatch.ref.column[2]) ),
    "MIN VH+VL mutations (aa)"  = min( !!as.name(aa_mismatch.ref.column[1]) + !!as.name(aa_mismatch.ref.column[2]) ),
    "MAX VH+VL mutations (aa)"  = max( !!as.name(aa_mismatch.ref.column[1]) + !!as.name(aa_mismatch.ref.column[2]) ),
    
    "MEAN CDR3 HC length (nt)" =  mean( !!as.name(cdr3aalength.ref.column[1]) * 3 ),
    "MEDIAN CDR3 HC length (nt)" =  median( !!as.name(cdr3aalength.ref.column[1]) * 3 ),
    "MIN CDR3 HC length (nt)" =  min( !!as.name(cdr3aalength.ref.column[1]) * 3),
    "MAX CDR3 HC length (nt)" =  max( !!as.name(cdr3aalength.ref.column[1]) * 3),
    "MEAN CDR3 LC length (nt)" =  mean( !!as.name(cdr3aalength.ref.column[2]) * 3),
    "MEDIAN CDR3 LC length (nt)" =  median( !!as.name(cdr3aalength.ref.column[2]) * 3 ),
    "MIN CDR3 LC length (nt)" =  min( !!as.name(cdr3aalength.ref.column[2]) * 3 ),
    "MAX CDR3 LC length (nt)" =  max( !!as.name(cdr3aalength.ref.column[2]) * 3 ),
    
    "MEAN CDR3 HC length (aa)" =  mean( !!as.name(cdr3aalength.ref.column[1]) ),
    "MEDIAN CDR3 HC length (aa)" =  median( !!as.name(cdr3aalength.ref.column[1]) ),
    "MIN CDR3 HC length (aa)" =  min( !!as.name(cdr3aalength.ref.column[1]) ),
    "MAX CDR3 HC length (aa)" =  max( !!as.name(cdr3aalength.ref.column[1]) ),
    "MEAN CDR3 LC length (aa)" =  mean( !!as.name(cdr3aalength.ref.column[2]) ),
    "MEDIAN CDR3 LC length (aa)" =  median( !!as.name(cdr3aalength.ref.column[2]) ),
    "MIN CDR3 LC length (aa)" =  min( !!as.name(cdr3aalength.ref.column[2]) ),
    "MAX CDR3 LC length (aa)" =  max( !!as.name(cdr3aalength.ref.column[2]) )
  ) 
  
  summary.list[['SHM patient and isotype']] <- unified.df %>% group_by(sampleid,ISOTYPE) %>% summarise( !!!summarise.expression )
  summary.list[['SHM by patient']] <- unified.df %>% group_by(sampleid) %>% summarise( !!!summarise.expression )
  summary.list[['SHM by isotype']] <- unified.df %>% group_by(ISOTYPE) %>% summarise( !!!summarise.expression )
 
  return( summary.list )
}

# create.summary.list <- function ( unified.df ){
#   summary.list <- list()
#   
#   # Rlang unevaluated expressions
#   summarise.expression <- quos(
#     "NUM SEQS" = as.list(n()),
#     "MEAN VH mutations (nt)" = mean( nt_mismatches_V_region...12 ),
#     "MEDIAN VH mutations (nt)" = median( nt_mismatches_V_region...12 ),
#     "MIN VH mutations (nt)"  = min( nt_mismatches_V_region...12 ),
#     "MAX VH mutations (nt)"  = max( nt_mismatches_V_region...12 ),
#     "MEAN VL mutations (nt)" = mean( nt_mismatches_V_region...48 ),
#     "MEDIAN VL mutations (nt)" = median( nt_mismatches_V_region...48 ),
#     "MIN VL mutations (nt)"  = min( nt_mismatches_V_region...48 ),
#     "MAX VL mutations (nt)"  = max( nt_mismatches_V_region...48 ),
#     "MEAN VH+VL mutations (nt)" = mean( nt_mismatches_V_region...12 + nt_mismatches_V_region...48 ),
#     "MEDIAN VH+VL mutations (nt)" = median( nt_mismatches_V_region...12 + nt_mismatches_V_region...48  ),
#     "MIN VH+VL mutations (nt)"  = min( nt_mismatches_V_region...12 + nt_mismatches_V_region...48  ),
#     "MAX VH+VL mutations (nt)"  = max( nt_mismatches_V_region...12 + nt_mismatches_V_region...48  ),
#     
#     "MEAN VH mutations (aa)" = mean( v_region_aa_mismatches...13 ),
#     "MEDIAN VH mutations (aa)" = median( v_region_aa_mismatches...13 ),
#     "MIN VH mutations (aa)"  = min( v_region_aa_mismatches...13 ),
#     "MAX VH mutations (aa)"  = max( v_region_aa_mismatches...13 ),
#     "MEAN VL mutations (aa)" = mean( v_region_aa_mismatches...49 ),
#     "MEDIAN VL mutations (aa)" = median( v_region_aa_mismatches...49 ),
#     "MIN VL mutations (aa)"  = min( v_region_aa_mismatches...49 ),
#     "MAX VL mutations (aa)"   = max( v_region_aa_mismatches...49 ),
#     "MEAN VH+VL mutations (aa)" = mean( v_region_aa_mismatches...13 + v_region_aa_mismatches...49 ),
#     "MEDIAN VH+VL mutations (aa)" = median( v_region_aa_mismatches...13 + v_region_aa_mismatches...49 ),
#     "MIN VH+VL mutations (aa)"  = min( v_region_aa_mismatches...13 + v_region_aa_mismatches...49 ),
#     "MAX VH+VL mutations (aa)"  = max( v_region_aa_mismatches...13 + v_region_aa_mismatches...49 ),
#     
#     "MEAN CDR3 HC length (nt)" =  mean( cdr3_aa_length...15 * 3 ),
#     "MEDIAN CDR3 HC length (nt)" =  median( cdr3_aa_length...15 * 3 ),
#     "MIN CDR3 HC length (nt)" =  min( cdr3_aa_length...15 * 3),
#     "MAX CDR3 HC length (nt)" =  max( cdr3_aa_length...15 * 3),
#     "MEAN CDR3 LC length (nt)" =  mean( cdr3_aa_length...51 * 3),
#     "MEDIAN CDR3 LC length (nt)" =  median( cdr3_aa_length...51 * 3 ),
#     "MIN CDR3 LC length (nt)" =  min( cdr3_aa_length...51 * 3 ),
#     "MAX CDR3 LC length (nt)" =  max( cdr3_aa_length...51 * 3 ),
#     
#     "MEAN CDR3 HC length (aa)" =  mean( cdr3_aa_length...15 ),
#     "MEDIAN CDR3 HC length (aa)" =  median( cdr3_aa_length...15 ),
#     "MIN CDR3 HC length (aa)" =  min( cdr3_aa_length...15 ),
#     "MAX CDR3 HC length (aa)" =  max( cdr3_aa_length...15 ),
#     "MEAN CDR3 LC length (aa)" =  mean( cdr3_aa_length...51 ),
#     "MEDIAN CDR3 LC length (aa)" =  median( cdr3_aa_length...51 ),
#     "MIN CDR3 LC length (aa)" =  min( cdr3_aa_length...51 ),
#     "MAX CDR3 LC length (aa)" =  max( cdr3_aa_length...51 )
#   ) 
#   
#   summary.list[['SHM patient and isotype']] <- unified.df %>% group_by(sampleid,ISOTYPE) %>% summarise( !!!summarise.expression )
#   summary.list[['SHM by patient']] <- unified.df %>% group_by(sampleid) %>% summarise( !!!summarise.expression )
#   summary.list[['SHM by isotype']] <- unified.df %>% group_by(ISOTYPE) %>% summarise( !!!summarise.expression )
#   
#   return( summary.list )
# }

plot.shm.and.cdr3.length = function( excel.file ) {
  project = gsub( "(\\S+?)_.*","\\1", basename(excel.file))
  cluster_type =  gsub( ".*_(\\S+?)_selected_columns.*","\\1", basename(excel.file))
  xl.obj = as.data.frame( read_excel( excel.file, sheet = "PROPER", skip = 1, col_names = T ) )
  
  unified.df = get.unified.df( ig.df = xl.obj )
  # get cdr3aa cols
  seqid.ref.column = grep('SEQUENCE_ID\\.\\.',names(unified.df),value=F)
  nt_mismatch.ref.column = grep('nt_mismatches_V_region\\.\\.',names(unified.df),value=F)
  cdr3aalength.ref.column = grep('cdr3_aa_length\\.\\.',names(unified.df),value=F)
  
  df = tibble(SEQUENCE_ID = unified.df[,seqid.ref.column[1]],
              patient = unified.df$sampleid,
              CLUSTERCLASS = unified.df$CLUSTERCLASS,
              ISOTYPE = unified.df$ISOTYPE,
              LIGHTCHAIN = unified.df$LIGHTCHAIN,
              "VH mutations (nt)" = unified.df[,nt_mismatch.ref.column[1]],
              "VL mutations (nt)" = unified.df[,nt_mismatch.ref.column[2]],
              "VH + VL mutations (nt)" = ( unified.df[,nt_mismatch.ref.column[1]] + unified.df[,nt_mismatch.ref.column[2]] ),
              "CDR3 HC length (aa)" = unified.df[,cdr3aalength.ref.column[1]],
              "CDR3 LC length (aa)" = unified.df[,cdr3aalength.ref.column[2]])
  
  df.melt = reshape2::melt( df )
  
  df.melt$ISOTYPE <- factor(df.melt$ISOTYPE,levels = c('IGM','IGG','IGA'))
  df.melt$ISOTYPE <- factor(df.melt$ISOTYPE,levels = c('IGG','IGM','IGA'))
  df.melt$ISOTYPE_CLUSTERCLASS = paste0(df.melt$ISOTYPE,'\n',df.melt$CLUSTERCLASS)
  df.melt$ISOTYPE_CLUSTERCLASS <- factor(df.melt$ISOTYPE_CLUSTERCLASS,  levels = c('IGM\nclone','IGM\nsinglet','IGG\nclone','IGG\nsinglet','IGA\nclone', 'IGA\nsinglet'))
  df.melt$ISOTYPE_CLUSTERCLASS <- factor(df.melt$ISOTYPE_CLUSTERCLASS,  levels = c('IGG\nclone','IGG\nsinglet','IGM\nclone','IGM\nsinglet','IGA\nclone', 'IGA\nsinglet'))
  
  # Get shm plots
  shm.list = get.shm.jitter.plot( df.melt )
  # Get CDR3 plots
  cdr3.list = get.cdr3.jitter.plot( df.melt )
  
  dir.create( output.dir, recursive = T, showWarnings = F )
  
  #write.table( file = paste0(output.dir,'/', patient, "_SHM_CDR3_table.txt"), x = shm.df, quote = F, row.names = F, col.names = T, sep = "\t" )
  
  plots = c(shm.list$plots,cdr3.list$plots)
  plots_size = c(shm.list$size,cdr3.list$size)
  
  for (n in names(plots)){
    svglite::svglite(file = paste0(output.dir, '/',project,'_',cluster_type,'_',n,"_jitter_plot.svg"), height = plots_size[[n]][1], width = plots_size[[n]][2])
    print(plots[[n]])
    dev.off()
    pdf(file = paste0(output.dir, '/',project,'_',cluster_type,'_',n,"_jitter_plot.pdf"), height = plots_size[[n]][1], width = plots_size[[n]][2])
    print(plots[[n]])
    dev.off()
  }
  
  shm.summary.list <- create.summary.list( unified.df )
  writexl::write_xlsx(shm.summary.list, path = paste0(output.dir, '/',project,'_',cluster_type,'_summary.xlsx'))
  
}

 
# Hydrophobicity plots
# -------------------------------------------------------------------------------

check.hydrophobicity = function ( seq.vec ) {
  gravy.score = hydrophobicity( as.character( gsub("X","", seq.vec )  ), scale = "Guy")
  # remove gravy.score == 0
  return( gravy.score )
}

get.boxplot.plot = function ( db.and.covid ) {
  x <- 'ISOTYPE'
  my_comparisons = combn(as.character(unique(db.and.covid[[x]])), 2, simplify = F)
  
  # Create outlier function
  check_outlier <- function(df, group_info){
    coef=1.5
    quantiles <- quantile(df$CDR3_IGH_GRAVY_SCORE, probs=c(0.25,0.75) )
    IQR <- quantiles[2] - quantiles[1]
    res <- (df$CDR3_IGH_GRAVY_SCORE < ( quantiles[1]- coef*IQR )) | (df$CDR3_IGH_GRAVY_SCORE > ( quantiles[2]+ coef*IQR ))
    df$outlier <- NA
    df$outlier[res] <- 'outlier'
    return(df)
  }
  
  # Create box plots by sample
  plot.by.sample <- list()
  samples  <- unique(db.and.covid$sampleid)
  samples <- grep(pattern='BCR\nDatabase',x = samples,invert = T, value = T)
  # reorder samples
  samples <- samples[order(samples)]
  for (sample in samples) {
    db.and.covid.sample <- db.and.covid %>% filter(sampleid %in% c(sample, 'BCR\nDatabase'))
    # Apply with data.table "by" method
    db.and.covid.outliers.sample = db.and.covid.sample %>% group_by(!!as.name(x)) %>% group_modify(check_outlier) %>% filter(!is.na(outlier)) %>% distinct()
    plot.by.sample[[sample]] <- ggplot( data = db.and.covid.sample, aes_string(x=x , y='CDR3_IGH_GRAVY_SCORE', fill='ISOTYPE')) +
      stat_boxplot( geom ='errorbar', linetype=1, width=0.3 ) +  #whiskers
      geom_boxplot(color = "black", outlier.shape = NA, width=0.4) +
      geom_point(data=db.and.covid.outliers.sample) +
      #coord_cartesian( ylim=c( round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2), round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2)) ) +
      ylab("GRAVY score") +
      scale_y_continuous(breaks=seq(round(min(db.and.covid$CDR3_IGH_GRAVY_SCORE), 2),  round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2), length.out = 10)) +
      xlab("") +
      ggtitle( paste0(sample,"\n Hydrophobicity score\n HC CDR3" )) +
      theme_bw() + 
      theme(legend.position = "none") +
      theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
      stat_compare_means(comparisons = my_comparisons,
                         label = "p.signif",
                         hide.ns = F) #+ # Add pairwise comparisons p-value
    #stat_compare_means(label.y =  max(xl.obj$CDR3_IGH_GRAVY_SCORE))  # Add global p-value
  }
  
  # arrange sample plots in 3 columns
  boxplot.score.by.sample <- ggarrange(plotlist = plot.by.sample, ncol = 3)
  
  # Create global boxplot
  db.and.covid.outliers = db.and.covid %>% group_by(!!as.name(x)) %>% group_modify(check_outlier) %>% filter(!is.na(outlier)) %>% distinct()
  boxplot.score <- ggplot( data = db.and.covid, aes_string(x=x , y='CDR3_IGH_GRAVY_SCORE', fill='ISOTYPE')) +
    stat_boxplot( geom ='errorbar', linetype=1, width=0.3 ) +  #whiskers
    geom_boxplot(color = "black", outlier.shape = NA, width=0.4) +
    geom_point(data=db.and.covid.outliers) +
    #coord_cartesian( ylim=c( round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2), round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2)) ) +
    ylab("GRAVY score") +
    scale_y_continuous(breaks=seq(round(min(db.and.covid$CDR3_IGH_GRAVY_SCORE), 2),  round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2), length.out = 10)) +
    xlab("") +
    ggtitle( paste0(paste0(samples,collapse = ', '),"\n Hydrophobicity score\n HC CDR3" )) +
    theme_bw() + 
    theme(legend.position = "none") +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       hide.ns = F) #+ # Add pairwise comparisons p-value
  
  plots = list(
    'hydrophobicity_score' = boxplot.score,
    'hydrophobicity_score_by_sample' = boxplot.score.by.sample
  )
  
  plots_size = list(
    'hydrophobicity_score' = c(8.5,11/3),
    'hydrophobicity_score_by_sample' = c(8.5, 11)
  )
  
  return(list(plots=plots,plots_size=plots_size))
  
}

plot.hidrophobicity = function( excel.file ) {
  xl.obj = NULL
  
  project = gsub( "(\\S+?)_.*","\\1", basename(excel.file))
  cluster_type =  gsub( ".*_(\\S+?)_selected_columns.*","\\1", basename(excel.file))
  obj = read_excel( excel.file, sheet = "PROPER", skip = 1, col_names = T, .name_repair = "universal")
  
  xl.obj = get.unified.df(obj)
  
  # remove sequences with empyt CDR3aa
  # get cdr3aa cols
  cdr3aa.ref.column = grep('cdr3_aa\\.\\.',names(xl.obj),value=F)
  xl.obj = xl.obj[which(!is.na(xl.obj[,cdr3aa.ref.column[1]]) & !is.na(xl.obj[,cdr3aa.ref.column[1]])),]

    # calculate GRAVY SCORE
  xl.obj$CDR3_IGH_GRAVY_SCORE = check.hydrophobicity(xl.obj[,cdr3aa.ref.column[1]])
  
  
  # Part 2: Calculating GRAVY scores for B cell database
  # x = read.table("~/Downloads/public-bcell-dataset-IGHV-cdr3aa.tsv")
  #x = read.table(argv$bcelldbtsv)
  #colnames(x) = c("SEQUENCE_ID...3", "cdr3_aa...13")
  #cat('Please wait... this analysis will take ~25 minutes to finish!')
  
  #gravy.score.db = apply(X = x, 1, check.hydrophobicity )
  #gravy.score.db.df = get.cdr3.score.df( gravy.score.db )
  
  #gravy.score.db.df$patient = "CDR3 from BCR Database"
  #gravy.score.db.df = gravy.score.db.df %>% rename(CDR3_IGH_GRAVY_SCORE =  "gravy.score.all.cdr3")
  
  # save RDS
  #saveRDS(gravy.score.db.df, file = paste0( argv$output, "/bcelldb.rds") )
  
  #read RDS
  gravy.score.db.df = readRDS(argv$rdsbcelldb)

  #set.seed(12345)
  #gravy.score.db.df <- gravy.score.db.df[sample(1:nrow(gravy.score.db.df),1000),]
  gravy.score.db.df$sampleid <- "BCR\nDatabase"
  gravy.score.db.df$patient <- "BCR\nDatabase"
  gravy.score.db.df$ISOTYPE <- "BCR\nDatabase"
  
  
  # Part 3: Concatenating database sequences data frame with our sequences database
  xl.obj.aux = xl.obj %>% select(sampleid, ISOTYPE, CDR3_IGH_GRAVY_SCORE) %>% mutate(patient='patients')

  #db.and.covid = plyr::rbind.fill(xl.obj.aux, gravy.score.db.df[sample(1:nrow(gravy.score.db.df),nrow(gravy.score.db.df)),])
  db.and.covid = rbind(xl.obj.aux, gravy.score.db.df)
  db.and.covid$sampleid = factor(db.and.covid$sampleid, levels=c("BCR\nDatabase", unique(xl.obj.aux$sampleid)))
  db.and.covid$patient = factor(db.and.covid$patient, levels=c("BCR\nDatabase", unique(xl.obj.aux$patient)))
  db.and.covid$ISOTYPE = factor(db.and.covid$ISOTYPE, levels=c("BCR\nDatabase", c('IGM','IGG','IGA')))

  
  p.list = get.boxplot.plot( db.and.covid )

  plots = p.list$plots
  plots_size = p.list$plots_size
  
  for (n in names(plots)){
    #svglite::svglite(file = paste0(output.dir, '/',project,'_',cluster_type,'_',n,"_boxplot.svg"), height = plots_size[[n]][1], width = plots_size[[n]][2])
    #print(plots[[n]])
    #dev.off()
    pdf(file = paste0(output.dir, '/',project,'_',cluster_type,'_',n,"_boxplot.pdf"), height = plots_size[[n]][1], width = plots_size[[n]][2])
    print(plots[[n]])
    dev.off()
  }
}


#################################################### Execution ####################################################
# excel.files = list.files( path = "~/Documents/Rockefeller/Davide/davide_analysis/covid_tables", full.names = TRUE, pattern = "xlsx" )
excel.files <- list.files( path = argv$input, full.names = TRUE, pattern = ".*_all_samples.*selected.*xlsx" )

output.dir <- argv$output

# Jitter plot separately per patient
mclapply(excel.files, plot.shm.and.cdr3.length)


# Jitter plot hydrophobicity
mclapply(excel.files, plot.hidrophobicity)
