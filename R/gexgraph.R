gexgraph <-
function (genome.test, scale, results.graph, W, cytobands)
{
   ###############################################################################
   #                                                                             #
   #                              GExGraph                                       #
   #                                                                             #
   #        INPUT: Genome.map                                                    #
   #        OUTPUT: graphics in pdf format                                       #
   #                                                                             #
   ###############################################################################

    write("Graph", file = "")
    lev = as.data.frame(genome.test)
    genome.test[as.numeric(genome.test[, "binome.test"]) > 0.05, "binome.test"] = ""
    genome.test[as.numeric(genome.test[, "Wilcoxon"]) > 0.05, "Wilcoxon"] = ""
    for (C in levels(lev[, "chr"]))
    {
        #temporary data matrix of each chr
        genome.temp = genome.test[genome.test[, "chr"] == C, ]
        
        #cytobands matrix
        cytobands_tmp = cytobands[cytobands[,"chr"]==C, c("cytoband", "debut", "fin")]       
        temp = cytobands_tmp
        cytobands_tmp = cbind(as.matrix(as.numeric(temp[,"debut"])), as.matrix(as.numeric(temp[,"fin"])))
        rownames(cytobands_tmp) = temp[,"cytoband"]
        colnames(cytobands_tmp) = c("debut", "fin")
        
        #maximum value of x = maximum value of total genes =>width of the graphic
        max = max(as.numeric(genome.temp[, "total"]))
        #maximum value of y = number of unit for the current chr = height of the graphic
        may = max(as.numeric(genome.temp[, "unit"]))
        
        write(paste("graph ", genome.temp[1, "chr"], sep = ""), file = "")
        
        colours = character()
        COL.ENS = "green"
        COL.list = "orange"
        COL.Haz = "black"
        COL.up = "red"
        COL.down = "blue"
        COL.regions = "red"
        COL.binome.test = "cyan"
        COL.Optmized = "orange"
        COL.Wilcoxon = "green"
        colours = c(COL.ENS, COL.list, COL.Haz, COL.up, COL.down, COL.regions, COL.binome.test, COL.Optmized, COL.Wilcoxon)
        
        y = as.matrix(genome.temp[, "unit"])
        
        #min and max of y axis
        min_up = min(as.numeric(genome.temp[, "up"]))
        min_down = min(as.numeric(genome.temp[, "down"]))
        min = min(min_down, min_up) - 1
        
        #replace the "1" which marks the regions of interest by the min value to draw the line of interesting region under the curves
        genome.temp[genome.temp[, "regions"] == 1, "regions"] = min
        
        #pvals selection <=0.05   
        genome.temp[genome.temp[, "binome.test"] != "", "binome.test"] = min - 0.5
        genome.temp[genome.temp[, "Optmized"] != "", "Optmized"] = min - 1
        genome.temp[genome.temp[, "Wilcoxon"] != "", "Wilcoxon"] = min - 1.5
        
        min = min - 1
        
        x = as.matrix(genome.temp[, c("total", "list", "hazard", "up", "down", "regions", "binome.test", "Optmized", "Wilcoxon")])
        
        pdf(file = paste(results.graph, "Chr", C, ".pdf", sep = ""), width = (may*0.8), height = (abs(min-2)+max)/2.5, family = "Helvetica", title = paste("chr", C, sep = ""))
        #pdf(file = paste(results.graph, "Chr", C, ".pdf", sep = ""), width = (may*1.5), height = (abs(min-2)+max)/2.5, family = "Helvetica", title = paste("chr", C, sep = ""))
        
        matplot(y, x, col = colours, type = "o", pch = 20, axes = FALSE, ylim = c(min, max), ylab = "", xlab = "")
        
        title(main = paste("Chr", C), cex.main = 3, font.main = 2, col.main = "black", font.lab = 2, col.lab = "black", xlab = "", ylab = "")
        title(xlab = paste(scale, " pb", sep = ""), line = 3.5, cex.lab = 3)
        title(ylab = "Gene frequency", line = 1.5, cex.lab = 3)
        axis(side = 2, at = min:max, tick = TRUE, lab = c(rep(min:max, 1)), cex.axis = 1)
        axis(side = 1, at = seq(1, (may + 1), by = 5), tick = TRUE, lab = c(seq(0, may, by = 5)), cex.axis = 1)
        
        #cytobands
        pair=0
        for(K in 1:nrow(cytobands_tmp))
        {
          if(pair == 0)
          {
            rect((cytobands_tmp[K, "debut"]/scale), (min-0.75), (cytobands_tmp[K, "fin"]/scale), (min - 1.25), col="black") 
            text(((((cytobands_tmp[K, "fin"]/scale) - (cytobands_tmp[K, "debut"]/scale))/2) + (cytobands_tmp[K, "debut"]/scale)) , (min-1), rownames(cytobands_tmp)[K], cex = 1, col="white")
            pair = 1
          } else
          {
            rect((cytobands_tmp[K, "debut"]/scale), (min-0.75), (cytobands_tmp[K, "fin"]/scale), (min - 1.25), col="white") 
            text(((((cytobands_tmp[K, "fin"]/scale) - (cytobands_tmp[K, "debut"]/scale))/2) + (cytobands_tmp[K, "debut"]/scale)) , (min-1), rownames(cytobands_tmp)[K], cex = 1, col="black")
           pair = 0
          }
        }
        
        #  Legend of curves
        leg.txt = character()
        leg.txt <- c(leg.txt, "Curves")
        leg.txt <- c(leg.txt, "")
        leg.txt <- c(leg.txt, "Ensembl genome")
        leg.txt <- c(leg.txt, "Tested list")
        leg.txt <- c(leg.txt, "Hazard estimation")
        leg.txt <- c(leg.txt, "Upregulated genes")
        leg.txt <- c(leg.txt, "Downregulated genes")
        COL.blank = "white"
        colours = c(COL.Haz, COL.blank, COL.ENS, COL.list, COL.Haz, COL.up, COL.down)
        legend(1, max, leg.txt, pch = "  -----", col = colours, cex = 3)
        
        #  Legend of statistical results
        leg.txt = character()
        leg.txt <- c(leg.txt, "Statistucal results")
        leg.txt <- c(leg.txt, "")
        leg.txt <- c(leg.txt, "Units where Tested list > Hazard estimation")
        leg.txt <- c(leg.txt, "Statistically validated")
        leg.txt <- c(leg.txt, "Optimized regions")
        leg.txt <- c(leg.txt, "Significant")
        colours = c(COL.Haz, COL.blank, COL.regions, COL.binome.test, COL.Optmized, COL.Wilcoxon)
        legend(15, max, leg.txt, pch = "  ----", col = colours, cex = 3)
        dev.off()
    }
}

