gexmap <-
function (genome = "homosapiens", scale = "", source = "", res = "", isGO = FALSE, isMAP = TRUE, lim_chi = 5, global_test_choice = 4, pcorrd = 1, pcorrv = 1)
{
    write("###############################################################################", file = "")
    write("#                                                                             #", file = "")
    write("#                              GExMap 1.1                                     #", file = "")
    write("#                                                                             #", file = "")
    write("###############################################################################", file = "")

	###############################################################################
	#                                                                             #
	#                              INPUT PARAMETERS                               #
	#                                                                             #
	#    scale = width of units (Mbp)                                             #
	#    path = Source folder                                                     #
	#    res = main results folder                                                #
	#                                                                             #
	###############################################################################
 
    library("multtest")
    if (scale == "")
    {
        scale = 1e+06
        write("Default scale: 1 unit = 1 000 000 bp", file = "")
    } else
    {
        write(paste("Scale unite =", scale, sep = ""), file = "")
    }
    #Default source folder
    if (source == "")
    {
        source = paste(getwd(), "/library/GExMap/data/", sep = "")
        write(paste("Default source folder : ", source, sep = ""), file = "")
    } else
    {
        write(paste("Custom source folder :", source, sep = ""), file = "")
    }
    #Results folder
    if (res == "")
    {
        res = paste(getwd(), "/GExMap", format(Sys.time(), "(%H-%M-%S) %a%d %b%Y"), sep = "")
        write(paste("Default results folder : ", res, sep = ""), file = "")
    } else
    {
        write(paste("Custom results folder :", res, sep = ""), file = "")
    }
    dir.create(res)
    res = paste(res, "/", sep = "")
    
   ###############################################################################
   #                                                                             #
   #                               UPLOAD LISTE                                  #
   #                                                                             #
   ###############################################################################
  
    write("Upload of list to analyze", file = "")
    if (source == "test")
    {
        data(list)
    } else
    {
        list = read.table(file = file.choose(), sep = "\t", header = FALSE, row.names = 1)
    }
    names = c(rownames(list)[1], as.vector(list[1, ]))
    col = as.matrix(rownames(list))
    col = as.matrix(col[2:nrow(col), ])
    list = cbind(col, as.matrix(list[2:nrow(list), ]))
    colnames(list) = names
    
    write(paste("Loading the ", genome, " genome Rdata file", sep = ""), file = "")
    if (source == "test")
    {
        data(data)
    } else
    {
        data = gexload.data(source, genome)
    }
    
    list.ENS = gexload.corr(data, list, scale, source, res)

   ######################################################################################################
   #                                                                                                    #
   #                             FORMATTING GENOME TABLE                                                #
   #                                                                                                    #
   #   chr = chromosome                                                                                 #
   #   unit = scale unit (million of bp...)                                                             #
   #   total = Nbre of genes by unit for the ENSEMBL genome                                             #
   #   list = Nbr of genes by unit for the tested list                                                  #
   #   hasard = Nbr of genes by unit for the tested gene list expected by chance                        #
   #                                                                                                    #                                                                           
   #   genome = matrix of all data for graphics                                                         #
   #   ENS = ENSEMBL gross data matrix for the current chromosome                                       #
   #   genome.temp = matrix of the number of ENSEMBL genes of each unit for the current chromosome      #   
   #   liste.chr = matrix of the number of genes from the list of each unit for the current chromosome  #
   #   list.ENS = all data about the gene list                                                          #
   #                                                                                                    #
   ######################################################################################################

    if (isMAP)
    {
        genome = matrix(ncol = 6, nrow = 0, dimnames = NULL)
        colnames(genome) = c("chr", "unit", "total", "list", "up", "down")
        
        # matrix of number of ENSEMBL genes by unit
        for (I in levels(factor(data[, "chr"])))   # chromosome par chromosome
        {
            write(paste("Treatment of Chr ", I, sep = ""), file = "")
            list.chr = list.ENS[list.ENS[, "chr"] == I, ]      #extracts list data for the current chromosome
            
            #If the list contains only one gene the defaul class is "character"
            if (!is.matrix(list.chr))
            {
                list.chr = t(as.matrix(list.chr))
            } else   # matrix of only one gene
            {
                list.chr = as.matrix(list.chr[order(list.chr[, "m"]), ])      #sorts conlumn by m
            }
            if (nrow(list.chr) > 5)   #If there is at least five gene on the current chromosome
            {
                genome.temp = gex.mapping(data, list.chr, scale, I, list.ENS)
            } else 
            {
                genome.temp = matrix(ncol = ncol(genome), nrow = 0)
            }
            genome = rbind(genome, genome.temp)
        }
        
        #cytobands matrix creation
        cytobands = matrix(0, ncol=4, nrow=0)
        colnames(cytobands) = c("chr", "cytoband", "debut", "fin")
        for (I in levels(as.data.frame(data)[,"chr"]))
        {
          chr_tmp = as.data.frame(data[data[,"chr"]==I,])
          for (J in levels(factor(chr_tmp[, "cytoband"])))
          {
             cytobands = rbind(cytobands, matrix(0, ncol=4, nrow=1))
             cytobands[nrow(cytobands), "chr"] = I
             cytobands[nrow(cytobands), "cytoband"] = J
             cytobands[nrow(cytobands), "debut"] = min(as.matrix(chr_tmp[chr_tmp[,"cytoband"]==J, "m"]))
             cytobands[nrow(cytobands), "fin"] = max(as.matrix(chr_tmp[chr_tmp[,"cytoband"]==J, "m"]))
          }
        }
        
        #TESTS using two computation of the hazard: global then specific for each one of the two test packages
        #Resulting files are placed in separated folders
       
        #Estimation of hazard global then specific
        
        chance.compute = matrix(0, ncol = 2, nrow = 2)
        colnames(chance.compute) = c("nhazard", "genome.graph")
        rownames(chance.compute) = c("specific method", "global method")
        chance.compute[1, 1] = sum(as.numeric(genome[genome[, "list"] != 0, "total"]))   #number of ENSEMBL genes presents in units where there is at least one gene frome the gene list
        chance.compute[2, 1] = sum(as.numeric(genome[, "total"]))
        chance.compute[1, 2] = paste(res, "specific method/", sep = "")
        chance.compute[2, 2] = paste(res, "global method/", sep = "")
        for (W in 1:2)
        {
            write(paste("Hazard computation ", rownames(chance.compute)[W], sep = ""), file = "")
            hazard = matrix(0, ncol = 1, nrow = nrow(genome))
            colnames(hazard) = "hazard"
            nhazard = chance.compute[W, 1]
            results.graph = chance.compute[W, 2]
            for (H in 1:nrow(genome))
            {
                hazard[H, "hazard"] = (sum(as.numeric(genome[, "list"]))/as.numeric(nhazard)) * as.numeric(genome[H, "total"])  #computes number of genes expected by chance
            }
            
            if (!dir.create(results.graph))
            {
                dir.create(results.graph)
            }
            
            genome.test = cbind(genome, hazard)
            
            #statistic method
            genome.test = gextest(genome.test, nhazard, results.graph, lim_chi, global_test_choice, pcorrd, pcorrv)
            write.table(genome.test, row.names = FALSE, sep = "\t", file = file.path(paste(results.graph, "genome.txt", sep = "")))
            gexgraph(genome.test, scale, results.graph, W, cytobands)
            write(paste("Results have been saved in the folder: ", results.graph, sep = ""), file = "")
        }
    }
    #GO
    if (isGO)
    {
        gexgo(list.ENS, source, res)
    }
}

