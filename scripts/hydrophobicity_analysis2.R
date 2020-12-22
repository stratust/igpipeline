#########################################################################################
####################################### Libraries ####################################### 
#########################################################################################

library( Peptides )
library( readxl )
library( dplyr )
library( ggpubr )
library( argparser )

p <- arg_parser("SHM analysis")

p <- add_argument(p, "--input", help="Excel file with all samples")

p <- add_argument(p, "--bcelldbtsv", help="TSV file containing the CDR3 of the B cell database")
p <- add_argument(p, "--rdsbcelldb", help="RDS file containing CDR3 score of the B cell database")

p <- add_argument(p, "--output", help="output to export the images", default="output.txt")

argv <- parse_args(p)
#argv$input = "/home/stratus/GitHub/igpipeline2_tbev/results/PATIENTS/igpairs/PATIENTS_all_samples_clonal_selected_columns.xlsx"
#argv$rdsbcelldb = "/home/stratus/GitHub/igpipeline2/database/bcelldb.rds"
#argv$output = "/home/stratus/GitHub/igpipeline2/results/PATIENTS/SHM"
#argv$output = "/home/stratus/TMP/"
#########################################################################################
####################################### Functions ####################################### 
#########################################################################################

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
    
    if (is.null(HC.df) && !is.na(seqid.vector)){
      HC.df = ig.df.na.rows.removed[,column.range]
      isotype <- rep(i,nrow(HC.df))
    } else if (!is.null(HC.df) && !is.na(seqid.vector)) {
      selected.rows = !is.na(seqid.vector)
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
    
    if (is.null(LC.df) && !is.na(seqid.vector)){
      LC.df = ig.df.na.rows.removed[,column.range]
      light.chain = rep(i, nrow(LC.df))
    } else if (!is.null(LC.df) && !is.na(seqid.vector)) {
      selected.rows = !is.na(seqid.vector)
      LC.df[selected.rows,] =  ig.df.na.rows.removed[selected.rows,column.range]
      light.chain[selected.rows] = i
    }
    
  }
  
  unified.df = cbind(data.frame(ISOTYPE=isotype,LIGHTCHAIN=light.chain),ig.df.na.rows.removed[,1:3],HC.df,LC.df)
  
  return(unified.df)
  
}

check.hydrophobicity = function ( seq.vec ) {
  gravy.score = hydrophobicity( as.character( gsub("X","", seq.vec )  ), scale = "Guy")
  # remove gravy.score == 0
  return( gravy.score )
}

get.boxplot.plot = function ( xl.obj ) {
  
  xl.obj$patient = factor(xl.obj$patient, levels = unique( xl.obj$patient) )
  xl.obj$sampleid = factor(xl.obj$sampleid, levels = unique( xl.obj$sampleid) )

  # boxplot.score = ggplot( data = xl.obj, aes(x=patient , y=CDR3_IGH_GRAVY_SCORE)  ) +
  #   geom_boxplot( aes (fill = patient), color = "black" ) +
  #   coord_cartesian( ylim=c( round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2), round( max(xl.obj$CDR3_IGH_GRAVY_SCORE) , 2)) ) +
  #   ylab("GRAVY score") +
  #   scale_y_continuous(breaks=c( seq(round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2),  round( max(xl.obj$CDR3_IGH_GRAVY_SCORE) , 2), 0.1))) +
  #   xlab("") +
  #   ggtitle( "Hydrophobicity score using HC CDR3" ) +
  #   theme_classic() + 
  #   theme(legend.position = "none")
  
  my_comparisons = combn(as.character(unique(db.and.covid$sampleid)), 2, simplify = F)
  
  boxplot.score = ggplot( data = db.and.covid, aes(x=sampleid , y=CDR3_IGH_GRAVY_SCORE)  ) +
    geom_boxplot( aes (fill = patient), color = "black" ) +
    #coord_cartesian( ylim=c( round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2), round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2)) ) +
    ylab("GRAVY score") +
    #scale_y_continuous(breaks=c( seq(round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2),  round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2), length.out = 5))) +
    xlab("") +
    ggtitle( "Hydrophobicity score using heavy chain CDR3" ) +
    theme_classic() + 
    theme(legend.position = "none") #+
    #stat_compare_means(comparisons = my_comparisons,
    #                   label = "p.signif",
    #                   hide.ns = F) #+ # Add pairwise comparisons p-value
    #stat_compare_means(label.y =  max(xl.obj$CDR3_IGH_GRAVY_SCORE))  # Add global p-value

  return(boxplot.score)
  
}


#########################################################################################
####################################### Execution ####################################### 
#########################################################################################

# Part 1: Calculating GRAVY scores for our sequences
# excel.files = list.files( path = "~/Documents/Rockefeller/Davide/davide_analysis/covid_tables", full.names = TRUE, pattern = "xlsx" )
excel.file = argv$input
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
gravy.score.db.df$sampleid <- "CDR3 from BCR Database"

# Part 3: Concatenating database sequences data frame with our sequences database
xl.obj.aux = xl.obj %>% select(sampleid, CDR3_IGH_GRAVY_SCORE) %>% mutate(patient='patients')
set.seed(12345)
db.and.covid = plyr::rbind.fill(xl.obj.aux, gravy.score.db.df[sample(1:nrow(gravy.score.db.df),1000000),])
#db.and.covid[ which(db.and.covid$patient != "CDR3 from BCR Database"), "patient"] = "patients"

boxplot.score = get.boxplot.plot( xl.obj = db.and.covid )
boxplot.score

# Part 4: Statistical tests and boxplot

# Shapiro-Wilk test was used to determine whether the GRAVY scores are normally distributed.

# Shapiro-Wilk test applied to the database sequences - 5000 sequences randomly selected
#shapiro.test( subset(db.and.covid, patient == "CDR3 from BCR Database")[sample(1:nrow(subset(db.and.covid, patient == "CDR3 from BCR Database")), 5000, replace=F), "CDR3_IGH_GRAVY_SCORE"] )

# Shapiro-Wilk test applied to our sequences only
#shapiro.test( subset(xl.obj, patient != "CDR3 from BCR Database")$CDR3_IGH_GRAVY_SCORE )

t = boxplot.score #+ stat_compare_means(paired = FALSE, method = "wilcox.test", label.y = 0.5)
pdf(file = paste0( argv$output, "/all_patients_collapsed_and_db_boxplot_2.pdf"), height = 12, width = 18 )
print(t)
dev.off()
