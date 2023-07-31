#' Circo plot of the connector clustering
#'
#'@description
#'	Generates the circos plot given two connector clustering
#'
#' @param dataFrom starting dataframe with one column ID and the second one the connector cluster membership
#' @param dataTo ending dataframe with one column ID and the second one the connector cluster membership
#' @param nameFrom name of the experiment to print above the starting cluster names 
#' @param nameTo name of the experiment to print above the ending cluster names 
#' @return
#'
#' @examples
#'
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'
#' @import circlize dplyr wesanderson viridis
#' @export

circo.generation = function(dataFrom,dataTo,nameFrom,nameTo,filename)
{
  # change the name of the columns 
  colnames(dataFrom) = c("ID","ClusterFrom")
  colnames(dataTo) = c("ID","ClusterTo")
  
  dataMerged = merge(dataFrom,dataTo,by = "ID")
  
  dataMerged2 = dataMerged %>% 
    tidyr::gather(key = "type",value = "cl", -ID) %>%
    dplyr::mutate(type = paste0(type,"_",cl))
  
  # color definition
  CLall = unique(dataMerged2$type)
  colAll = rep("white",length(CLall))
  names(colAll) = sort(CLall)
  colAll[grep("From",CLall)] = wes_palette("GrandBudapest1", n = length(grep("From",CLall)),type = "discrete")
  colAll[grep("To",CLall)] = viridis(length(grep("To",CLall)))
  
  colors = c("#d8b365","#5ab4ac")
  names(colors) = c(nameFrom,nameTo)
  Fromcol = colors[nameFrom]
  Tocol = colors[nameTo]
  
  # starting the definition of relations in the circo
  
  cl1 = CLall[grep(x = CLall,pattern = "To_")]
  cl2 = CLall[grep(x = CLall,pattern = "From_")]
  
  df = dataMerged %>% 
    mutate(ClusterFrom = paste0("ClusterFrom_",ClusterFrom),
           ClusterTo = paste0("ClusterTo_",ClusterTo)) %>%
    group_by(ClusterFrom, ClusterTo) %>% 
    tally()
  
  group = structure(gsub(paste0("_(",paste0(LETTERS,collapse = "|"),")$"),"", CLall), names = CLall)
  col = colAll[names(group)]
  
  # starting the circo:
  pdf(paste0(filename,".pdf"),width = 6,height = 6)
  
  circos.clear()
  circos.par(start.degree = -100,track.margin = c(mm_h(4), 0))
  
  chordDiagram(df,
               group = group, 
               big.gap = 20,
               small.gap = 1,
               annotationTrack = "grid", 
               grid.col = col,
               annotationTrackHeight = c(mm_h(3)),
               preAllocateTracks = list(
                 track.height = mm_h(4),
                 track.margin = c(mm_h(2), 0)
               )
  )
  
  # Adding the %
  circos.track(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    xplot = get.cell.meta.data("xplot")
    
    circos.lines(xlim, c(min(ylim), min(ylim)), lty = 3) # dotted line
    by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.2, 0.5)
    circos.lines(c(0,0),
                 c(min(ylim), min(ylim)- 0.15), lty = 1)
    for(p in seq(by, 1, by = by)) {
      circos.text(p*(xlim[2] - xlim[1]) + xlim[1], min(ylim) - 0.6, 
                  paste0(p*100, "%"), cex = 0.6,
                  adj = c(0.6, 0), niceFacing = TRUE )
      circos.lines(c(p*(xlim[2] - xlim[1]) + xlim[1],p*(xlim[2] - xlim[1]) + xlim[1]),
                   c(min(ylim), min(ylim)- 0.15), lty = 1)
    }
  }, bg.border = NA)
  # Adding the name of the clusters
  circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index = gsub("(ClusterFrom|ClusterTo)_","", get.cell.meta.data("sector.index"))
    #sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim),col = "white", sector.index, cex = 0.6, niceFacing = TRUE)
  }, bg.border = NA)
  # Adding the name of the two sets
  highlight.sector(cl1,
                   track.index = 1,
                   col = Tocol, 
                   text = nameTo,
                   cex = 0.8, text.col = "white", niceFacing = TRUE)
  highlight.sector(cl2, 
                   track.index = 1,
                   col = Fromcol, 
                   text = nameFrom, cex = 0.8, text.col = "white", niceFacing = TRUE)
  
  circos.clear() 
  
  dev.off()
  
}