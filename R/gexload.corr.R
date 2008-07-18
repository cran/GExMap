
gexload.corr<-function(data, list, scale, source, res)
{
  ######################################
  # upload of correspondances file
  ######################################
  
  mtype = unlist(strsplit(colnames(list)[1], ","))[2] #column name = information about microarray technology and / or type
  ttype = unlist(strsplit(colnames(list)[1], ","))[1]
  if(length(agrep(",", colnames(list)[1], ignore.case = TRUE))!=0)  #virgule suivit du type de puce
  { 
    write(paste("Loading the ", ttype, " ", mtype, " file", sep=""), file="");
    Rdata.dir = file.path(paste(source, ttype, "-", mtype, ".Rdata", sep=""));
  }else{
    write(paste("Loading the ", ttype, ".Rdata file", sep=""), file="");
    Rdata.dir = file.path(paste(source, ttype, ".Rdata", sep=""));
  }
  
  if(source=="test")
  {
     data(corr)
  }else{
	  if(!file.exists(Rdata.dir)) 
	  {
	    write(sprintf("ERROR: file %s does not exist in the default directory", source), file="");
	    write(paste("Load manually the ", ttype, "-", mtype, ".Rdata file", sep="") ,file="")
	    load(file=file.choose());
	  } else
	  {
	    load(file=Rdata.dir)
	    write("Rdata file Loaded")
	  }
	}
  
  ######################################
  # Matrix of correspondences
  ######################################
  idt=0
  absents_corr = matrix(0, ncol=1, nrow=0)  #matrix of the identifiers not found in the data.Rdata matrix
  write("Looking for correspondences", file="");
  
  id = matrix(0, ncol=2, nrow=nrow(list))
  for(n in 1:nrow(list))
  {
    temp = as.matrix(corr[corr[,"probes"]==list[n, 1], "ensembl"])
    if(nrow(temp)==0)
    {
      #if temp is empty, the id frome the user list is unknown in the corr file
      absents_corr = rbind(absents_corr, list[n, 1])
    }else{
      #if the id is presnt in the corr file: OK
      id[n,1] = corr[corr[,"probes"]==list[n, 1], "ensembl"]
      id[n,2] = as.numeric(list[n, "expression"])
    }
    write(paste(n, " in ", nrow(list), sep=""), file="")
  }
  write("Correspondences OK", file="");
  if(nrow(absents_corr)!=0)
  {
    write.table(absents_corr, row.names=FALSE, sep="\t", file = paste(res, "absents_corr.txt", sep=""));
    write(paste(nrow(absents_corr), " identifiers unknown with no ENSEMBL correspondences have been placed in the absents_corr.txt file", sep=""), file="")
  }
  idt=1
  # Suppression of the ENSG00000000000 identifiers which corresponds to unidentified probes
  id = as.matrix(id[id[,1]!="ENSG00000000000",])
  write(paste("Nbre of identified probes=", nrow(id) , sep=""), file="");
  
  ###############################################################################
  #                                                                             #
  #                               VALIDATION                                    #
  #                                                                             #
  ###############################################################################
  
  if(idt!=0)   #Known identifier type
  {
    ################################################## looking for localisation information
    write("Collecting localisation data", file="");
    id_map = matrix(0, ncol=6, nrow=nrow(id))
    tbl_scale = matrix(0, ncol=1, nrow=nrow(id))
    colnames(id_map)=c("chr", "m", "cytoband", "name", "ensembl", "GO")
    absents_data = matrix(0, ncol=1, nrow=0)
    for(n in 1:nrow(id))
    {
      temp = as.matrix(data[data[,"ensembl"]==id[n,1], c("chr", "m", "cytoband", "name", "ensembl", "GO")])
      if(nrow(temp)==0)
      {
        absents_data = rbind(absents_data, id[n,1])
      }else{
        id_map[n,] = temp
        tbl_scale[n,1] = round(as.numeric(id_map[n, 2])/scale)
      }
    }
    if(nrow(absents_data)!=0)
    {
      write.table(absents_data, row.names=FALSE, sep="\t", file = paste(res, "absents_data.txt", sep=""));
      write(paste(nrow(absents_data), " identifiers unknown with no data have been placed in the absents_data.txt file", sep=""), file="")
    }
    id_map = cbind(id_map, id[,2])
    colnames(id_map)=c("chr", "m", "cytoband", "name", "ensembl", "GO", "expression")
    
    write("Localisation data OK", file="");
    write("Gestion of the duplicates", file="");
    id_map = as.matrix(id_map[order(id_map[,"ensembl"]),])
    r=as.matrix(levels(factor(id_map[,"ensembl"])))                   #Extraction of ENSEMBL id and filtering duplicates
    write(paste(nrow(r), " Probes after identification and duplicates filtering", sep=""), file="");
    map = matrix(0, ncol=ncol(id_map), nrow = nrow(r))
    colnames(map) = colnames(id_map)
    for(J in 1:nrow(r))
    {
      temp = matrix(0, ncol=ncol(id_map), nrow=0)
      temp = rbind(temp, id_map[id_map[,"ensembl"]==r[J,1],])
      map[J,] = temp[1,]
      #print(J)
    }
    id_map=map
    write("Saving data", file="");
    write.table(id_map, row.names=FALSE, sep="\t", file = file.path(paste(res, "report.txt", sep="")));
    
    return(id_map)
  } else   #Unknown identifier type
  {
    write("Unknown identifier type", file="");
    return(FALSE)
  }
}
