#################################################### Libraries ####################################################
library( readxl )
library( ggplot2 )
library( ggpubr )
library( reshape2 )
library( argparser )

p <- arg_parser("SHM analysis")

p <- add_argument(p, "--input", help="directory where the excel files are located")

p <- add_argument(p, "--output", help="output to export the images", default="output.txt")

argv <- parse_args(p)

#################################################### Function ####################################################
get.SHM = function ( ig.df ) {
  
  ig.df.na.rows.removed = ig.df[rowSums(is.na(ig.df)) != ncol(ig.df), ]
  
  # Gathering Heavy-Chain information
  VH.nt.mutations.ref.column = "nt_mismatches_V_region...11"
  cdr3.HC.aa.length.ref.column = "cdr3_aa_length...14"
  VH.nt.mutations = ig.df.na.rows.removed[, VH.nt.mutations.ref.column ]
  cdr3.HC.aa.length = ig.df.na.rows.removed[, cdr3.HC.aa.length.ref.column ]
  
  
  # Gathering Light-Chain information
  VK.nt.mutations.ref.column = "nt_mismatches_V_region...23"
  VK.nt.mutations = ig.df.na.rows.removed[, VK.nt.mutations.ref.column ]
  VL.nt.mutations.ref.column = "nt_mismatches_V_region...35"
  VL.nt.mutations = ig.df.na.rows.removed[, VL.nt.mutations.ref.column ]
  VLC.nt.mutations = VL.nt.mutations
  VLC.nt.mutations[ !is.na(VK.nt.mutations) ] = VK.nt.mutations[!is.na(VK.nt.mutations)]
  
  cdr3.VK.aa.length.ref.column = "cdr3_aa_length...26"
  cdr3.VK.aa.length = ig.df.na.rows.removed[, cdr3.VK.aa.length.ref.column ]
  cdr3.VL.aa.length.ref.column = "cdr3_aa_length...38"
  cdr3.VL.aa.length = ig.df.na.rows.removed[, cdr3.VL.aa.length.ref.column ]
  cdr3.VLC.aa.length = cdr3.VL.aa.length
  cdr3.VLC.aa.length[ !is.na(cdr3.VK.aa.length) ] = cdr3.VK.aa.length[!is.na(cdr3.VK.aa.length)]
  cdr3.VLC.aa.length = cdr3.VLC.aa.length
  
  
  shm.df = data.frame(SEQUENCE_ID = as.character( ig.df.na.rows.removed$`SEQUENCE_ID...3`),
                      VH.nt.mutations = VH.nt.mutations,
                      VL.nt.mutations = VLC.nt.mutations,
                      SHM = ( VH.nt.mutations + VLC.nt.mutations ),
                      cdr3.HC.length = cdr3.HC.aa.length,
                      cdr3.LC.length = cdr3.VLC.aa.length)
  
  return(shm.df)
  
}

get.jitter.plot = function ( shm.df, patient, output ) {
  
  # shm.df.melt = reshape2::melt( shm.df )
  shm.df.melt = shm.df
  shm.df.melt$patient = patient
  
  first.plot = subset( shm.df.melt, variable %in% c("VH.nt.mutations","VL.nt.mutations") )
  first.plot$variable= as.character(first.plot$variable)
  
  mean.VH.nt.mutations = mean(subset( shm.df.melt, variable %in% c("VH.nt.mutations") )$value)
  mean.VL.nt.mutations = mean(subset( shm.df.melt, variable %in% c("VL.nt.mutations") )$value)
  
  nb.mutations.plot = ggplot() +
    geom_point(data = first.plot, aes(x=variable , y=value, fill="black"), 
               position = position_jitterdodge( jitter.width=0.3, dodge.width = 0.7), size=2)+
    geom_path(aes(x=c(0.7,1.3),y =c(mean.VH.nt.mutations, mean.VH.nt.mutations)), size=1, color = "black")+
    geom_path(aes(x=c(1.7,2.3),y=c(mean.VL.nt.mutations, mean.VL.nt.mutations)),size=1, color = "black")+
    coord_cartesian( ylim=c(-3, 85), expand=F) +
    ylab("# of mutations") +
    scale_y_continuous(breaks=c(seq(0, 85, 5))) +
    xlab("") +
    # ggtitle( paste0("Patient: ", patient) ) +
    theme_classic() + 
    theme(legend.position = "none")
  
  
  second.plot = subset( shm.df.melt, variable %in% c("cdr3.HC.length","cdr3.LC.length") )
  second.plot$variable= as.character(second.plot$variable)
  
  mean.cdr3.HC.length = mean(subset( shm.df.melt, variable %in% c("cdr3.HC.length") )$value)
  mean.cdr3.LC.length = mean(subset( shm.df.melt, variable %in% c("cdr3.LC.length") )$value)
  
  cdr3.plot = ggplot() +
    geom_point(data = second.plot, aes(x=variable , y=value, fill="black"),
               position = position_jitterdodge( jitter.width=0.3, dodge.width = 0.7), size=2)+
    geom_path(aes(x=c(0.7,1.3),y =c(mean.cdr3.HC.length, mean.cdr3.HC.length)), size=1, color = "black")+
    geom_path(aes(x=c(1.7,2.3),y=c(mean.cdr3.LC.length, mean.cdr3.LC.length)),size=1, color = "black")+
    coord_cartesian( ylim=c(-3, 35), expand = F) +
    scale_y_continuous(breaks=c(seq(0, 35, 5))) +
    ylab("CDR3 length") +
    xlab("") +
    theme_classic() +
    theme(legend.position="none")
  
  svglite::svglite(file = paste0(output.dir, '/',patient, "_jitter_plot.svg"), height = 6, width = 7 )
  figure = ggarrange(nb.mutations.plot, cdr3.plot, legend = "none")
  annotate_figure(figure, top = text_grob(paste0("Patient ", patient), color = "black", face = "bold", size = 14) )
  print( annotate_figure(figure, top = text_grob(paste0("Patient ", patient), color = "black", face = "bold", size = 14) ) )
  dev.off()
  
  pdf(file = paste0(output.dir,'/' ,patient, "_jitter_plot.pdf"), height = 6, width = 7 )
  figure = ggarrange(nb.mutations.plot, cdr3.plot, legend = "none")
  print( annotate_figure(figure, top = text_grob(paste0("Patient ", patient), color = "black", face = "bold", size = 14) ) )
  dev.off()
  
  return(figure)
  
}
#################################################### Execution ####################################################
# excel.files = list.files( path = "~/Documents/Rockefeller/Davide/davide_analysis/covid_tables", full.names = TRUE, pattern = "xlsx" )
excel.files = list.files( path = argv$input, full.names = TRUE, pattern = ".*selected.*xlsx" )

plots.list = list()
all.patients.shm = NULL

output.dir = argv$output 

# Jitter plot separately per patient
for (excel.file in excel.files) {
  
  patient = gsub( "(\\S+?)_.*","\\1", basename(excel.file))
  xl.obj = as.data.frame( read_excel( excel.file, sheet = "PROPER", skip = 1, col_names = T ) )
  
  shm.df = get.SHM( ig.df = xl.obj )
  shm.df.melt = reshape2::melt( shm.df )
  shm.df.melt$patient = patient
  
  # dir.create("~/Documents/Rockefeller/Davide/davide_analysis/results", recursive = T, showWarnings = F)
  dir.create( output.dir, recursive = T, showWarnings = F )
  write.table( file = paste0(output.dir,'/', patient, "_SHM_CDR3_table.txt"), x = shm.df, quote = F, row.names = F, col.names = T, sep = "\t" )
  
  plots.list[[excel.file]] = get.jitter.plot( shm.df = shm.df.melt, patient = patient, output = output.dir )
  
  all.patients.shm = rbind(all.patients.shm, shm.df.melt)
  
}

all.fig = ggarrange(plotlist = plots.list)
svglite::svglite(file = paste0(output.dir, "/all_patients_jitter_plot.svg"), height = 9, width = 16 )
print(all.fig)
dev.off()

# Jitter plot of all patients combined
figure = get.jitter.plot(shm.df = all.patients.shm, patient = paste0(unique(all.patients.shm$patient), collapse = ", "), output = output.dir )

# mean.VH.nt.mutations = mean(subset( all.patients.shm, variable %in% c("VH.nt.mutations") )$value)
# mean.VL.nt.mutations = mean(subset( all.patients.shm, variable %in% c("VL.nt.mutations") )$value)
# mean.cdr3.HC.length = mean(subset( all.patients.shm, variable %in% c("cdr3.HC.length") )$value)
# mean.cdr3.LC.length = mean(subset( all.patients.shm, variable %in% c("cdr3.LC.length") )$value)
# 
# values = data.frame(mean.VH.nt.mutations, mean.VL.nt.mutations, mean.cdr3.HC.length, mean.cdr3.LC.length)
# write.table(x = values, file = "~/Documents/Rockefeller/Davide/davide_analysis/results/all_patients_mean.txt", quote = F, col.names = T, row.names = F, sep = "\t")
