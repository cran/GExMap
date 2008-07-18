
gexgo<-function(list.ENS, source, res)
{
  write(paste("Analyzing GO identifiers", sep=""), file="")
  data(go)
  
  #bind two columns to the GO data matrix
  go = cbind(go, matrix(0, ncol=1, nrow=nrow(go))) #one column of 0 for the scores of each GO id
  go = cbind(go, matrix("", ncol=1, nrow=nrow(go))) #one column for the list of genes relatives to the GO id
  colnames(go)=c("id", "ontology", "type", "score", "genes")
  gengo = list.ENS[list.ENS[,"GO"]!="0",]   #filters only genes with GO id
  write(paste(nrow(gengo), " genes have GO identifiers", sep=""), file="")
  
  #comput GO id scrores and create corresponding genes lists
  for(I in 1:nrow(gengo))
  {
    write(paste(I, " of ", nrow(gengo), sep=""), file="")
    gotemp = as.matrix(unlist(strsplit(gengo[I,6], ","))) #convert GO id lists (strings with GO id separated by ",") in matrix
    for(J in 1:nrow(gotemp))
    {
      go[go[,1]==gotemp[J,],4] = as.numeric(go[go[,1]==gotemp[J,],4]) +1  #compute scores
      go[go[,1]==gotemp[J,],5] = paste(go[go[,1]==gotemp[J,],5], paste(gengo[I,4], "(", gengo[I, 5], ")", sep=""), sep=",") #concatenates genes names
    }
  }
  
  resume<- file(paste(res, "resume.txt", sep=""), "w+")
  cat(paste("Size of input list: ", nrow(list.ENS), sep=""), file=resume)
  cat(paste("\n", nrow(gengo), " genes have GO identifiers", sep=""), file=resume)
  
  #filters GO id with score=0
  resgo = go[go[,4]!="0",]
  resgoP = resgo[resgo[,3]=="P",]
  resgoF = resgo[resgo[,3]=="F",]
  resgoC = resgo[resgo[,3]=="C",]
  
  orderP=as.matrix(order(as.matrix(as.numeric(resgoP[,4])), decreasing = TRUE))  #sort the Biological Processes by score
  tP = as.matrix(unlist(strsplit(resgoP[orderP[1,],5], ",")))  #matrix of the longuest list of genes for a GO identifier 
  tP = as.matrix(tP[2:nrow(tP),]) #remove the first row which is empty
  detailP = matrix("", ncol=nrow(orderP), nrow=nrow(tP))  #matrix of all genes for each GO id
  sumP = sum(as.numeric(resgoP[,4]))   #Sum of genes for this GO type
  resP = matrix(0, ncol=4, nrow=0) 
  prop = matrix(0, ncol=2, nrow=nrow(orderP))   #part of each GO id in the GO type
  for(I in 1:nrow(orderP))
  {
    resP = rbind(resP, resgoP[orderP[I,],1:4])
    prop[I,1] = round((as.numeric(resP[resP[,1]== resgoP[orderP[I,],1], 4])/sumP)*100, 2)
    prop[I,2] = paste(resP[I,2], " (", prop[I,1], "%,", resP[I,4], ")", sep="")
    t = as.matrix(unlist(strsplit(resgoP[orderP[I,],5], ",")))
    temp = matrix("", ncol=1, nrow=(nrow(detailP)-nrow(t)+1))
    t = rbind(as.matrix(t[2:nrow(t),]), temp)
    detailP[,I] = t
  }
  resP = cbind(resP, prop)
  colnames(detailP) = as.vector(resP[,1])
  
  orderF=as.matrix(order(as.matrix(as.numeric(resgoF[,4])), decreasing = TRUE))  #sort the Molecular Function by score
  tF = as.matrix(unlist(strsplit(resgoF[orderF[1,],5], ",")))  #matrix of the longuest list of genes for a GO identifier 
  tF = as.matrix(tF[2:nrow(tF),]) #remove the first row which is empty
  detailF = matrix("", ncol=nrow(orderF), nrow=nrow(tF))  #matrix of all genes for each GO id
  sumF = sum(as.numeric(resgoF[,4]))   #Sum of genes for this GO type
  resF = matrix(0, ncol=4, nrow=0) 
  prop = matrix(0, ncol=2, nrow=nrow(orderF))   #part of each GO id in the GO type
  for(I in 1:nrow(orderF))
  {
    resF = rbind(resF, resgoF[orderF[I,],1:4])
    prop[I,1] = round((as.numeric(resF[resF[,1]== resgoF[orderF[I,],1], 4])/sumF)*100, 2)
    prop[I,2] = paste(resF[I,2], " (", prop[I,1],"%, ", resF[I,4], ")", sep="")
    t = as.matrix(unlist(strsplit(resgoF[orderF[I,],5], ",")))
    temp = matrix("", ncol=1, nrow=(nrow(detailF)-nrow(t)+1))
    t = rbind(as.matrix(t[2:nrow(t),]), temp)
    detailF[,I] = t
  }
  resF = cbind(resF, prop)
  colnames(detailF) = as.vector(resF[,1])
  
  orderC =as.matrix(order(as.matrix(as.numeric(resgoC[,4])), decreasing = TRUE))  #sort the Molecular Function by score
  tC = as.matrix(unlist(strsplit(resgoC[orderC[1,],5], ",")))  #matrix of the longuest list of genes for a GO identifier 
  tC = as.matrix(tC[2:nrow(tC),]) #remove the first row which is empty
  detailC = matrix("", ncol=nrow(orderC), nrow=nrow(tC))  #matrix of all genes for each GO id
  sumC = sum(as.numeric(resgoC[,4]))   #Sum of genes for this GO type
  resC = matrix(0, ncol=4, nrow=0) 
  prop = matrix(0, ncol=2, nrow=nrow(orderC))   #part of each GO id in the GO type
  for(I in 1:nrow(orderC))
  {
    resC = rbind(resC, resgoC[orderC[I,],1:4])
    prop[I,1] = round((as.numeric(resC[resC[,1]== resgoC[orderC[I,],1], 4])/sumC)*100, 2)
    prop[I,2] = paste(resC[I,2], " (", prop[I,1],"%,", resC[I,4], ")", sep="")
    t = as.matrix(unlist(strsplit(resgoC[orderC[I,],5], ",")))
    temp = matrix("", ncol=1, nrow=(nrow(detailC)-nrow(t)+1))
    t = rbind(as.matrix(t[2:nrow(t),]), temp)
    detailC[,I] = t
  }
  resC = cbind(resC, prop)
  colnames(detailC) = as.vector(resC[,1])
  
  res = paste(res, "GO.results", sep="")
  write(paste("Go results folder : ", res, sep=""), file="");
  if(!dir.create(res))
  {
    dir.create(res)
  }
  
  resP = cbind(resP[,1:2], resP[,4:6])
  colnames(resP) = c("ID GO", "Biological Process", "Score", "Proportion", "labels")
  resF = cbind(resF[,1:2], resF[,4:6])
  colnames(resF) = c("ID GO", "Molecular Functions", "Score", "Proportion", "labels")
  resC = cbind(resC[,1:2], resC[,4:6])
  colnames(resC) = c("ID GO", "Cellular Component", "Score", "Proportion", "labels")
  
  write.table(resP, row.names=FALSE, sep="\t", file = file.path(paste(res, "/GO.BP.txt", sep="")));
  write.table(detailP, row.names=FALSE, sep="\t", file = file.path(paste(res, "/detail.BP.txt", sep="")));
  write.table(resC, row.names=FALSE, sep="\t", file = file.path(paste(res, "/GO.CC.txt", sep="")));
  write.table(detailC, row.names=FALSE, sep="\t", file = file.path(paste(res, "/detail.CC.txt", sep="")));
  write.table(resF, row.names=FALSE, sep="\t", file = file.path(paste(res, "/GO.MF.txt", sep="")));
  write.table(detailF, row.names=FALSE, sep="\t", file = file.path(paste(res, "/detail.MF.txt", sep="")));
  
  lim=1   #minimum number genes linked with a GO to be mapped
  x = as.vector(as.numeric(resP[as.numeric(resP[,4])>=lim,4]))
  labels = as.vector(resP[as.numeric(resP[,4])>=lim,"labels"])
  pdf(file = file.path(paste(res, "/GO_BP.pdf", sep="")), width = 150, height = 80, family = "Helvetica", title = "Biological Processes")
  pie(x, labels, col=palette(rainbow(100)), cex=8)
  text(0,1, paste("Biological Processes (", sum(x), "%)", sep=""), cex=15)
  dev.off()
  
  write(paste(nrow(resgoP), " Biological Processes have have been identified in the list", sep=""), file="")
  cat("\n", paste(nrow(resgoP), " Biological Processes have have been identified in the list", sep=""), file=resume)
  cat("\n", paste(length(x), " GO id represent ", sum(x), "%", sep=""), file=resume)
  
  x = as.vector(as.numeric(resF[as.numeric(resF[,4])>=lim,4]))
  labels = as.vector(resF[as.numeric(resF[,4])>=lim,"labels"])
  pdf(file = file.path(paste(res, "/GO_MF.pdf", sep="")), width = 150, height = 80, family = "Helvetica", title = "Molecular Functions")
  pie(x, labels, col=palette(rainbow(100)), cex=8)
  text(0,1, paste("Molecular Functions, (", sum(x), "%)", sep=""), cex=15)
  dev.off()
  
  write(paste(nrow(resgoF), " Molecular Functions have have been identified in the list", sep=""), file="")
  cat("\n", paste(nrow(resgoF), " Molecular Functions have have been identified in the list", sep=""), file=resume)
  cat("\n", paste(length(x), " GO id represent ", sum(x), "%", sep=""), file=resume)
  
  x = as.vector(as.numeric(resC[as.numeric(resC[,4])>=lim,4]))
  labels = as.vector(resC[as.numeric(resC[,4])>=lim,"labels"])
  pdf(file = file.path(paste(res, "/GO_CC.pdf", sep="")), width = 150, height = 80, family = "Helvetica", title = "Cellular Component")
  pie(x, labels, col=palette(rainbow(100)), cex=8)
  text(0,1, paste("Cellular Component (", sum(x), "%)", sep=""), cex=15)
  dev.off()
  
  write(paste(nrow(resgoC), " Cellular Component have have been identified in the list", sep=""), file="")
  cat("\n", paste(nrow(resgoC), " Cellular Component have have been identified in the list", sep=""), file=resume)
  cat("\n", paste(length(x), " GO id represent ", sum(x), "%", sep=""), file=resume)
  
  close(resume)
  write(paste("GO analyze OK", sep=""), file="")
}
