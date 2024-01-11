library(dplyr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(writexl)
library(argparser, quietly = TRUE)
library(ggplot2)
library(svglite)
library(reshape2)

COLORS <- c(
  '#FF0000',
  '#FF8000',
  '#FFFF00',
  '#80FF00',
  '#00FF00',
  '#00FF80',
  '#00FFFF',
  '#0080FF',
  '#0000FF',
  '#7F00FF',
  '#FF00FF',
  '#FF007F'
)

parse_excel <- function(file) {
  data <- read_excel(file, skip = 1)
  # remove NAs
  data <- data[complete.cases(data[, c(1, 2)]), ]
  data <- data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::rename(IGA=SEQUENCE_ID...4, IGM=SEQUENCE_ID...67, IGG=SEQUENCE_ID...130) %>%
    select(sampleid, cluster_id, IGA, IGM, IGG) %>%  gather(isotype, value, IGA:IGG)
  data <- data[complete.cases(data[, c(1, 4)]), ]
  data$isotype <- paste(data$sampleid, data$isotype, sep='\n')
  data <- data %>% dplyr::group_by(cluster_id,isotype) %>% tally(name='n_seqs')
}

add_start_end <- function(data){
  #sort data by new cluster size
  data <- data[ order( data[['n_seqs']], data[['shared']], -as.numeric( data[['cluster_id']] ), decreasing = T), ]
  n_seqs <- data$n_seqs
  # add one to the first element
  #n_seqs[1] = n_seqs[1] + 1 
  data$end <- cumsum(n_seqs)
  data$start <- c(0, head(data$end, -1))
  return(data)
}


# Data must have all isotypes sequences for one patient
index_clone_colors <- function(data) {

  data.single <- data %>% filter(n_seqs == 1)
  data.clones <- data %>% filter(n_seqs > 1)
  
  # sort and keep the largest cluster_id
  #data.clones.unique <-  data.clones %>% arrange(desc(n_seqs)) %>% distinct(cluster_id, .keep_all = TRUE) #keep all columns
  data.clones.unique <-  data.clones %>% arrange(desc(n_seqs)) %>% distinct(cluster_id) # keep only cluster_id
  
  n_clones <- nrow(data.clones.unique)
  
  # get number of color replicates (n_clones/colors)
  n_replicates <- ceiling(n_clones/length(COLORS))
  
  colors <- rep(COLORS,n_replicates)[1:n_clones]
  names(colors) <- data.clones.unique$cluster_id
  
  # index color and cluster_id
  for (i in 1:nrow(data)){
    cluster_id <- data[i, 'cluster_id'] 
    data[i,'color'] = colors[as.character(cluster_id)]
  }
  
  data[which(is.na(data$color)),'color'] <- 'white'
  
  # find shared clones/singles
  data.shared  = data %>% group_by(cluster_id) %>% summarise(isotypes=length(isotype)) %>% filter(isotypes > 1)
  data$shared <- 0
  data[which(data$cluster_id %in% data.shared$cluster_id),'shared'] <- 1
  # make all white
  data$color_shared <- 'white'
  # make all clones gray
  data$color_shared[which(data$n_seqs > 1)] <- 'grey'
  
  # get shared sequences position in the dataframe
  shared_index <- which(data$cluster_id %in% data.shared$cluster_id)
  data$color_shared[shared_index] <- data$color[shared_index]

  return(data) 
}


# Make the plot
piechart_plot <- function(isotype, isotypes.list, shared = FALSE){
  data <- isotypes.list[[isotype]]
  fill_color <- NULL
  if (shared == TRUE){
    fill_color <- data$color_shared   
  }
  else {
    fill_color <- data$color
  }
  
  names(fill_color) <- data$cluster_id
  
  single_start <- data %>% filter(n_seqs == 1) %>% pull(start)
  single_end <- data %>% filter(n_seqs == 1) %>% pull(end)
  
  total_seq <- sum(data$n_seqs)
  total_seq_clones <- sum(data %>% filter(n_seqs > 1) %>% pull(n_seqs))
  percent_seq_clones <- round(total_seq_clones/total_seq * 100)
  
  # Get shared singles
  shared.singles <- data.frame()
  if (shared == TRUE){
    shared.singles <- data %>% filter(n_seqs == 1 & shared == 1)
  }

  single.df <- data.frame(ymax=max(single_end), ymin=min(single_start) + nrow(shared.singles))
  
  p <- ggplot(data, aes(ymax=end, ymin=start, xmax=4, xmin=3)) +
    geom_rect(aes(color = as.factor(cluster_id), fill=as.factor(cluster_id)), color='black', size = 1)

  if (total_seq > total_seq_clones + nrow(shared.singles)) {
    p <- p + geom_rect(data=single.df, aes(ymax = ymax, ymin = ymin, xmax=4, xmin = 3), color='black', fill = 'white', size=1) # single white chunk
  }

    # Draw black bar representing percent of clones
    clones.df <- filter(data, n_seqs > 1)
    if (nrow(clones.df) > 0) {
        p <- p + geom_rect(data = clones.df,  aes(xmin=4.1,xmax=4.2), fill = 'black')
    }
    p <- p + coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2.3, 4.4)) + # Try to remove that to see how to make a pie chart
    annotate(geom = 'text', x = 2.3, y = 0.5, label = sum(data$n_seqs), family = 'sans', fontface = 'bold', size = 16) + # Add number at the center
    annotate(geom = 'text', x = 4.4, y = 0, label = paste0(percent_seq_clones,'%'), family = 'sans', fontface = 'bold', size = 10) + # Add percentage number on top of black circle
    ggtitle(isotype) +
    theme_void() +
    scale_fill_manual(values = fill_color) +
    #scale_color_manual(values = border_color) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25, family = 'Helvetica', face ='bold'), 
      legend.position="none",
      text=element_text(family="Helvetica")
      ) # format plot title

    return(p)
}


create_piecharts <- function(files = NULL, filename = NULL) {
  
  patient.df <- parse_excel(files)
  patient.df$sampleid = gsub('\n.*','',patient.df$isotype)
  
  # split by patient/sample
  sampleid.list <- split(patient.df,patient.df$sampleid)
  sampleid.list <- lapply(sampleid.list, index_clone_colors)
  
  all_isotypes.df <- bind_rows(sampleid.list)
  
  # split by patient+isotype
  isotypes.list <- split(all_isotypes.df,all_isotypes.df$isotype)
  isotypes.list <- lapply(isotypes.list, add_start_end)
  
  all_isotypes.df <- bind_rows(isotypes.list)
  
  summary.df <- all_isotypes.df %>% mutate(type = case_when( n_seqs  > 1 ~ 'n_clones', n_seqs == 1 ~ 'n_singles'  )) %>% group_by(isotype,type) %>% tally()
  summary.shared.df <- all_isotypes.df %>% filter(shared == 1) %>% mutate(type = case_when( n_seqs  > 1 ~ 'n_shared_clones', n_seqs == 1 ~ 'n_shared_singles'  )) %>% group_by(isotype,type) %>% tally()
  summary.all.df <- dcast(rbind(summary.df,summary.shared.df), isotype ~ type)
  summary.all.df[is.na(summary.all.df)] <- 0
  summary.all.df$isotype <- gsub('\n','_',summary.all.df$isotype)

  isotypes.plot.list.shared <- lapply(names(isotypes.list), piechart_plot, isotypes.list = isotypes.list, shared = TRUE)
  names(isotypes.plot.list.shared) <- names(isotypes.list)
  
  for (i in names(isotypes.plot.list.shared)) {
    isotype = gsub('\n', '_', i, perl = TRUE)
    pdf(file = paste0(filename, '_', isotype, '_shared.pdf'))
    print(isotypes.plot.list.shared[[i]])
    dev.off()
    svglite(file = paste0(filename, '_', isotype, '_shared.svg'))
    print(isotypes.plot.list.shared[[i]])
    dev.off()
  }
  
  isotypes.plot.list.not.shared <- lapply(names(isotypes.list), piechart_plot, isotypes.list = isotypes.list, shared = FALSE)
  names(isotypes.plot.list.not.shared) <- names(isotypes.list)
  
  for (i in names(isotypes.plot.list.not.shared)) {
    isotype = gsub('\n', '_', i, perl = TRUE)
    pdf(file = paste0(filename, '_', isotype, '.pdf'))
    print(isotypes.plot.list.not.shared[[i]])
    dev.off()
    
    fonts <- list(sans = "Helvetica", serif = "Times New Roman")
    svglite(file = paste0(filename, '_', isotype, '.svg'), width = 7, height = 7, pointsize = 12, )
    print(isotypes.plot.list.not.shared[[i]])
    dev.off()
  }
  
  # write the excel file
  write_xlsx(isotypes.list, path = paste0(filename,'_piechart_info.xlsx'))
  write_xlsx(summary.all.df, path = paste0(filename,'_piechart_summary.xlsx'))
  
  
}

main <- function() {
  # Create a parser
  p <- arg_parser("Parse IgParse.pl output and create piecharts")
  
  # Add command line arguments
  p <- add_argument(p, "--input_file", help = "The excel file")
  p <- add_argument(p, "--output_prefix", help = "output prefix")

  #Parse the command line arguments
  argv <- parse_args(p)
  #argv$input_file <- list.files(pattern = ".*all.*_clonal.xlsx", path = '/home/stratus/GitHub/igpipeline2_timepoint//results/PATIENTS/igpairs/', full.names = T)
  #argv$output_prefix <- '/home/stratus/TMP/out'
  
  # check if input and codefile arguments are not NA
  if (is.na(argv$input_file) || is.na(argv$output_prefix)) {
    print(p)
    quit()
  }
  

  # For strict files COV
  create_piecharts(argv$input_file, filename = argv$output_prefix)
}


# Call main function
main()




