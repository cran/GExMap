gex.mapping <-
function (data, list.chr, scale, I, list.ENS)
{
   ###############################################################################
   #                                                                             #
   #                                 MAPPING                                     #
   #                                                                             #
   #                           List, Genome, Hazard                              #
   #                                                                             #
   ###############################################################################
  
    ENS = cbind(as.matrix(as.numeric(data[data[, "chr"] == I, "chr"])), as.matrix(as.numeric(data[data[, "chr"] == I, "m"]))) # extracts data for the current chr
    colnames(ENS) = c("chr", "m")
    ENS = as.matrix(ENS[order(ENS[, "m"]), ])       # sort m by asc
    size = max(ENS[, "m"])                          # max valure of m
    size = floor(size/scale) + 1                    # scaling, heigth = number of unit for the chr
    
    genome.temp = matrix(0, ncol = 3, nrow = size)   # matrix of the number of ENSEMBL genes by unit)
    genome.temp[, 1] = seq(1, size)
    chr = matrix(I, ncol = 1, nrow = size)
    genome.temp = cbind(chr, genome.temp)
    colnames(genome.temp) = c("chr", "scale", "ensembl", "list")
    
    # treatment of the matrix list.chr
    for (J in 1:nrow(ENS))                          # for each gene of the current chr
    {
        genome.temp[floor(as.integer(ENS[J, "m"])/scale) + 1, "ensembl"] = as.numeric(genome.temp[floor(as.integer(ENS[J, "m"])/scale) + 1, "ensembl"]) + 1
    }
    
    for (J in 1:nrow(list.chr))
    {
        genome.temp[floor(as.integer(list.chr[J, "m"])/scale) + 1, "list"] = as.numeric(genome.temp[floor(as.integer(list.chr[J, "m"])/scale) + 1, "list"]) + 1
    }
    up = matrix(0, ncol = 1, nrow = nrow(genome.temp)) #upregulated genes
    colnames(up) = "up"
    down = matrix(0, ncol = 1, nrow = nrow(genome.temp)) #downregulated genes
    colnames(down) = "down"
    for (J in 1:nrow(list.chr))                     # round of current chromosome
    {
        if (as.numeric(list.chr[J, "expression"] > 0))
        {
            up[floor(as.integer(list.chr[J, "m"])/scale) + 1, "up"] = up[floor(as.integer(list.chr[J, "m"])/scale) + 1, "up"] - 1
        }
        if (as.numeric(list.chr[J, "expression"] < 0))
        {
            down[floor(as.integer(list.chr[J, "m"])/scale) + 1, "down"] = down[floor(as.integer(list.chr[J, "m"])/scale) + 1, "down"] - 1
        }
    }
    genome.temp = cbind(genome.temp, up, down)
    return(genome.temp)
}

