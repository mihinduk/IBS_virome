
sel_cols_plot_temp <- sel_cols_plot[-2]

pdf("./figure_output/circos plot correlations cohort only.pdf", width=6.5, height=6.8) 
par(mar=c(0,0,0,0))

# edit initialising parameters
circos.par(track.margin = c(0, 0.01), # adjust bottom and top margin
           track.height = 0,
           cell.padding=c(0,0,0,0))
#optionally add gap.after=0.01
#gap after also messes up distance between viruses and other elements, could be made into a vector

chordDiagramFromDataFrame(circos_data, grid.col = grid.col_input, col=sel_cols_plot_temp[factor(dataset_group)], 
                          scale=F, annotationTrack = "grid", 
                          preAllocateTracks = list(track.height = 0.4, track.margin = c(mm_h(0), 0)))
# max(strwidth(unlist(dimnames(circos_data))))
#track margin is the spaceing between the track and the label. 
#


circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  temp_test <- length(grep("s__", sector.name)) > 0
  #test against list of colors
  if(temp_test == TRUE) {
    circos.text(mean(xlim), ylim[1], sector.name, font=3, facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5), col= "black", cex=0.35)
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex=0.35)
  }
}, bg.border = NA)


#KO terms and metabolites do not have to be italics

#add legend for connections
y = -0.85
for (i in 1:length(names(type_list))) {
  text(0.85,y, names(type_list)[i], col=rev(sel_cols_plot_temp)[i], pos = 1, cex=0.7)
  y <- y - 0.06
}
text(0.85, 0.8*y, "connection types", pos = 1, cex=0.8, font=2)


#add legend presence absence
y = 1.1
for (i in 1:4) {
  text(0.85 ,y, c("unique to Healthy", "unique to IBS-C", "unique to IBS-D")[i], col=sel_colors[i], 
       pos = 1, cex=0.7)
  y <- y - 0.06
}

dev.off()

circos.clear()

