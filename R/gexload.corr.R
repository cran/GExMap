gexload.corr <-
function (data, list, scale, source, res)
{
######################################
# upload of correspondances file
######################################

    mtype = unlist(strsplit(colnames(list)[1], ","))[2]
    ttype = unlist(strsplit(colnames(list)[1], ","))[1]
    if (length(agrep(",", colnames(list)[1], ignore.case = TRUE)) != 0)
    {
        write(paste("Loading the ", ttype, " ", mtype, " file", sep = ""), file = "")
        Rdata.dir = file.path(paste(source, ttype, "-", mtype, ".Rdata", sep = ""))
    } else
    {
        write(paste("Loading the ", ttype, ".Rdata file", sep = ""), file = "")
        Rdata.dir = file.path(paste(source, ttype, ".Rdata", sep = ""))
    }

    if (source == "test")
    {
        data(corr)
    } else
    {
        if (!file.exists(Rdata.dir))
        {
            write(sprintf("ERROR: file %s does not exist in the default directory", source), file = "")
            write(paste("Load manually the ", ttype, "-", mtype,".Rdata file", sep = ""), file = "")
            load(file = file.choose())
        } else
        {
            load(file = Rdata.dir)
            write("Rdata file Loaded")
        }
    }

    ######################################
    # Matrix of correspondences
    ######################################
    
    absents_corr = matrix(0, ncol = 1, nrow = 0)  #matrix of the identifiers not found in the data.Rdata matrix
    write("Looking for correspondences", file = "")
    #deb = format(Sys.time(), "%H:%M:%S")
    id = as.matrix(unlist(apply(list, 1, function(x) as.matrix(corr[corr[, "probes"] == x[1], "ensembl"]))))
    #write(paste("debut:", deb, ", fin:", format(Sys.time(), "%H:%M:%S"), sep=""), file="")
    id = cbind(id, as.matrix(rep(1, nrow(id))))
    
    write("Correspondences OK", file = "")
    if (nrow(absents_corr) != 0)
    {
      write.table(absents_corr, row.names = FALSE, sep = "\t", file = paste(res, "absents_corr.txt", sep = ""))
      write(paste(nrow(absents_corr), " identifiers unknown with no ENSEMBL correspondences have been placed in the absents_corr.txt file", sep = ""), file = "")
    }

    # Suppression of the ENSG00000000000 identifiers which corresponds to unidentified probes
    id = as.matrix(id[id[, 1] != "ENSG00000000000", ])
    write(paste("Nbre of identified probes=", nrow(id), sep = ""), file = "")
       
   ###############################################################################
   #                                                                             #
   #                               VALIDATION                                    #
   #                                                                             #
   ###############################################################################

    # looking for localisation information
    write("Collecting localisation data", file = "")
    tbl_scale = matrix(0, ncol = 1, nrow = nrow(id))
    absents_data = matrix(0, ncol = 1, nrow = 0)
    id_map = matrix(0, ncol=ncol(data), nrow=0)
    id_map = t(as.matrix(apply(id, 1, function(x) data[data[,"ensembl"]==x[1],])))
    id_map =  cbind(id_map, as.matrix(rep(1, nrow(id_map))))
    colnames(id_map) = c("chr", "cytoband", "ensembl", "name", "m", "GO", "expression")        
    id_map = id_map[,c("chr", "m", "cytoband", "name", "ensembl", "GO", "expression")]
    
    write("Localisation data OK", file = "")
    write("Saving data", file = "")
    write.table(id_map, row.names = FALSE, sep = "\t", file = file.path(paste(res, "report.txt", sep = "")))
    
    return(id_map)
}

