#########################################################################################
####################################### Libraries ####################################### 
#########################################################################################

library( Peptides )
library( readxl )
library( dplyr )
library( ggpubr )
library( argparser )

p <- arg_parser("SHM analysis")

p <- add_argument(p, "--input", help="directory where the excel files are located")

p <- add_argument(p, "--bcelldbtsv", help="TSV file containing the CDR3 of the B cell database")

p <- add_argument(p, "--output", help="output to export the images", default="output.txt")

argv <- parse_args(p)

#########################################################################################
####################################### Functions ####################################### 
#########################################################################################

get.df.information = function ( xl.obj ) {
  
  xl.obj = xl.obj[rowSums(is.na(xl.obj)) != ncol(xl.obj), ]
  xl.obj[ is.na(xl.obj)] = ""
  
  xl.obj = xl.obj[, c("cluster_id", "num_seqs",
                      "SEQUENCE_ID...3",
                      "translation_of_corrected_input_from_start...8", "cdr3_aa...13",
                      "translation_of_corrected_input_from_start...20", "cdr3_aa...25",
                      "translation_of_corrected_input_from_start...32", "cdr3_aa...37"
                      )]
  
  return( xl.obj )
}

check.hydrophobicity = function ( seq ) {

  gravy.score = NULL
  if (seq["SEQUENCE_ID...3"] != "COV047_P4_IgG_51-P1369") {
    
      gravy.score = hydrophobicity( as.character( gsub("X","", seq["cdr3_aa...13"] )  ), scale = "Guy")
  }
  
  return( paste0(seq["SEQUENCE_ID...3"], " SCORE: ", gravy.score) )
  
}

get.boxplot.plot = function ( xl.obj ) {
  
  xl.obj$patient = factor(xl.obj$patient, levels = unique( xl.obj$patient) )

  # boxplot.score = ggplot( data = xl.obj, aes(x=patient , y=CDR3_IGH_GRAVY_SCORE)  ) +
  #   geom_boxplot( aes (fill = patient), color = "black" ) +
  #   coord_cartesian( ylim=c( round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2), round( max(xl.obj$CDR3_IGH_GRAVY_SCORE) , 2)) ) +
  #   ylab("GRAVY score") +
  #   scale_y_continuous(breaks=c( seq(round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2),  round( max(xl.obj$CDR3_IGH_GRAVY_SCORE) , 2), 0.1))) +
  #   xlab("") +
  #   ggtitle( "Hydrophobicity score using HC CDR3" ) +
  #   theme_classic() + 
  #   theme(legend.position = "none")
  
  boxplot.score = ggplot( data = db.and.covid, aes(x=patient , y=CDR3_IGH_GRAVY_SCORE)  ) +
    geom_boxplot( aes (fill = patient), color = "black" ) +
    coord_cartesian( ylim=c( round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2), round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2)) ) +
    ylab("GRAVY score") +
    scale_y_continuous(breaks=c( seq(round(min(xl.obj$CDR3_IGH_GRAVY_SCORE), 2),  round( max(db.and.covid$CDR3_IGH_GRAVY_SCORE) , 2), length.out = 5))) +
    xlab("") +
    ggtitle( "Hydrophobicity score using heavy chain CDR3" ) +
    theme_classic() + 
    theme(legend.position = "none")

  return(boxplot.score)
  
}

get.cdr3.score.df = function ( gravy.score.all.cdr3 ) {
  
  names.gravy.score.all.cdr3 = gsub("^(\\S+) SC.*","\\1", gravy.score.all.cdr3 )
  gravy.score.all.cdr3 = as.numeric( gsub(".*SCORE\\:\\s+(\\S+)","\\1", gravy.score.all.cdr3 ) )
  names(gravy.score.all.cdr3) = names.gravy.score.all.cdr3
  gravy.score.all.cdr3 = as.data.frame(gravy.score.all.cdr3)
  gravy.score.all.cdr3 = na.omit(gravy.score.all.cdr3) # REMOVING NA DUE TO THE NON-EXISTING VALUE FOR COV047_P4_IgG_51-P1369
  
  return(gravy.score.all.cdr3)
  
}


#########################################################################################
####################################### Execution ####################################### 
#########################################################################################

# Part 1: Calculating GRAVY scores for our sequences
# excel.files = list.files( path = "~/Documents/Rockefeller/Davide/davide_analysis/covid_tables", full.names = TRUE, pattern = "xlsx" )
excel.files = list.files( path = argv$input, full.names = TRUE, pattern = "*.selected.*xlsx" )
xl.obj = NULL
for (excel.file in excel.files) {
  
  patient = gsub( "(\\S+?)_.*","\\1", basename(excel.file))
  obj = as.data.frame( read_excel( excel.file, sheet = "PROPER", skip = 1, col_names = T ) )
  
  obj = get.df.information( xl.obj = obj )
  
  xl.obj = rbind(xl.obj, obj)
  
}

gravy.score.all.cdr3 = apply(X = xl.obj, 1, check.hydrophobicity)
gravy.score.all.cdr3 = get.cdr3.score.df( gravy.score.all.cdr3 )

xl.obj = xl.obj %>% inner_join( gravy.score.all.cdr3 %>% tibble::rownames_to_column('sequence.id'), by = c("SEQUENCE_ID...3" = "sequence.id") )
xl.obj = xl.obj %>% rename(CDR3_IGH_GRAVY_SCORE = gravy.score.all.cdr3)
xl.obj$patient = gsub("^(COV.*?)_.*", "\\1", xl.obj$SEQUENCE_ID...3)

# Part 2: Calculating GRAVY scores for B cell database
# x = read.table("~/Downloads/public-bcell-dataset-IGHV-cdr3aa.tsv")
x = read.table(argv$bcelldbtsv)
colnames(x) = c("SEQUENCE_ID...3", "cdr3_aa...13")
cat('Please wait... this analysis will take ~25 minutes to finish!')

gravy.score.db = apply(X = x, 1, check.hydrophobicity )
gravy.score.db.df = get.cdr3.score.df( gravy.score.db )
gravy.score.db.df$patient = "CDR3 from BCR Database"
gravy.score.db.df = gravy.score.db.df %>% rename(CDR3_IGH_GRAVY_SCORE =  "gravy.score.all.cdr3")

# Part 3: Concatenating database sequences data frame with our sequences database
db.and.covid = plyr::rbind.fill(xl.obj, gravy.score.db.df)
db.and.covid[ which(db.and.covid$patient != "CDR3 from BCR Database"), "patient"] = "COVID-19 patients"

boxplot.score = get.boxplot.plot( xl.obj = db.and.covid )

# Part 4: Statistical tests and boxplot

# Shapiro-Wilk test was used to determine whether the GRAVY scores are normally distributed.

# Shapiro-Wilk test applied to the database sequences - 5000 sequences randomly selected
shapiro.test( subset(db.and.covid, patient == "CDR3 from BCR Database")[sample(1:nrow(subset(db.and.covid, patient == "CDR3 from BCR Database")), 5000, replace=F), "CDR3_IGH_GRAVY_SCORE"] )

# Shapiro-Wilk test applied to our sequences only
shapiro.test( subset(xl.obj, patient != "CDR3 from BCR Database")$CDR3_IGH_GRAVY_SCORE )

t = boxplot.score + stat_compare_means(paired = FALSE, method = "wilcox.test", label.y = 0.5)
pdf(file = paste0( argv$output, "/all_patients_collapsed_and_db_boxplot_2.pdf"), height = 6, width = 9 )
print(t)
dev.off()
