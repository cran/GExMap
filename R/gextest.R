gextest<-function (genome.test, nhazard, results.graph, lim_chi, global_test_choice, pcorrd, pcorrv)
{
    #########################################################################################################
    #
    #	The test take into acount all the chr and compares reel effectives of our list to the presence probability issued from the ENSEMBL genome.
    #	Probabilities are computed chr by chr
    #	Sum of presence probabilities of one chr is 1
    #	The program tests each chr at progressives scaling: from 1 to 10 Mb by unit. 
    #
    #########################################################################################################


######## Global tests to select chromosomes
   
    #round of genome chr by chr

    lev = levels(factor(genome.test[, "chr"]))
    max = 20
    g_test = matrix(0, ncol = max, nrow = length(lev))   #matrix Nchr*11 will stock the results of the tests for each chr and for each scale
    rownames(g_test) = lev
    colnames(g_test) = seq(1, max, 1)
    
    #g_test_graph: only pval to be reported graphically
    g_test_graph = g_test
    
    global_tests = matrix(0, ncol = 5, nrow = length(lev))
    colnames(global_tests) = c("chr", "pval", "scale", "min", "W")
    global_tests[, 1] = as.matrix(lev)
    
    lev = as.data.frame(genome.test)
    
    local_test_results = matrix(0, ncol = 12, nrow = 0)
    colnames(local_test_results) = c("chr", "unit", "total", "list", "up", "down", "hazard", "theo", "regions", "binome.test", "Optmized", "Wilcoxon")
    
    write(paste("Global tests", sep = ""), file = "")
    for (I in levels(lev[, "chr"]))
    {
      write(paste("Chr ", I, sep = ""), file = "")
      gentest = genome.test[genome.test[, "chr"] == I, ]  # extract data of crh
      if(nrow(gentest)!=0)
      {      
          #Wilcoxon global test comparaing real distribution (gentest[,"list"]) to expected distribution (gentest[,"hazard"])
          p = 1
          p = try(as.numeric(wilcox.test(as.matrix(as.numeric(gentest[, "list"])), as.matrix(as.numeric(gentest[, "hazard"])), paired = TRUE, exact = FALSE, correct = FALSE)[3]), silent=TRUE)
          if (p <= 0.05)
          {
              global_tests[global_tests[, "chr"] == I, 5] = p
          }
          
          nhazard = sum(as.numeric(gentest[, "total"]))    #N = toltal ENSEMBL genes of the chr
          theo = matrix(0, ncol = 1, nrow = nrow(gentest))     #matrix of probabilities of presence
          colnames(theo) = "theo"
          for (K in 1:nrow(gentest)) #cration of the matrix of presence probabilities
          {
              theo[K, ] = as.numeric(gentest[K, "total"])/nhazard
          }
          
          gentest = cbind(gentest, theo)
          for (N in 1:max)
          {
              # concatenates units
              if (N >= 2)
              {
                  if (floor(nrow(gentest)/N) * N == nrow(gentest))
                  {
                    x = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo_tmp = as.matrix(gentest[, "theo"])
                    x_tmp = as.matrix(gentest[, "list"])
                  } else     #si le compte n'est pas juste
                  {
                    x = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo_tmp = as.matrix(gentest[, "theo"])
                    theo_tmp[(floor(nrow(gentest)/N) * N), ] = sum(as.numeric(gentest[((floor(nrow(gentest)/N)) * N):nrow(gentest), "theo"]))
                    theo_tmp = as.matrix(theo_tmp[1:((floor(nrow(gentest)/N)) * N), ])
                    x_tmp = as.matrix(gentest[, "list"])
                    x_tmp[(floor(nrow(gentest)/N) * N), ] = sum(as.numeric(gentest[((floor(nrow(gentest)/N)) * N):nrow(gentest), "list"]))
                    x_tmp = as.matrix(x_tmp[1:((floor(nrow(gentest)/N)) * N), ])
                  }
                  for (K in 1:nrow(x)) #concatenates units
                  {
                    x[K, ] = sum(as.numeric(x_tmp[(K * N - N + 1):(K * N), ]))
                    theo[K, ] = sum(as.numeric(theo_tmp[(K * N - N + 1):(K * N), ]))
                  }
              } else
              {
                  x = as.matrix(as.numeric(gentest[, "list"]))
                  theo = as.matrix(as.numeric(gentest[, "theo"]))
              }
              #test
              p = suppressWarnings(try(as.numeric(chisq.test(x, p = theo)[3]), silent=TRUE))
              if (is.nan(p) | (p> 0.05))
              {
                  g_test[I, N] = ""
                  g_test_graph[I, N] = ""
              } else
              {
                  g_test_graph[I, N] = 1
                  g_test_graph[I, N] = suppressWarnings(try(as.numeric(chisq.test(x, p = theo)[3]), silent=TRUE))
                  g_test[I, N] = paste(g_test_graph[I, N], "(", min(x), ")", sep = "")
                  #test of the minimal number of genes by unit
                  if ((min(x) >= lim_chi) & (global_tests[global_tests[, 1] == I, 2] == 0))
                  {
                    global_tests[global_tests[, 1] == I, 2] = g_test[I, N]
                    global_tests[global_tests[, 1] == I, 3] = N
                    global_tests[global_tests[, 1] == I, 4] = min(x)
                  }
              }
          }
          
          #interesting region detection by a local test if the chr is interesting
          local_test = matrix(0, ncol = (ncol(genome.test) + 3), nrow = 0)
          colnames(local_test) = c(colnames(genome.test), "theo", "regions", "binome.test")
          
          # user choice of CHI/Wilcoxon/CHI-Wilcoxon
          pass_to_local = FALSE
          # 1- At least CHI is OK
          if ((global_tests[global_tests[, 1] == I, 2] != 0) & (global_test_choice == 1))
          {
              pass_to_local = TRUE
              write(paste("1- CHI is OK", sep = ""), file = "")
          }
          # 2- At least Wilcoxon is OK
          if ((global_tests[global_tests[, 1] == I, 5] != 0) & (global_test_choice == 2))
          {
              pass_to_local = TRUE
              write(paste("2- At least Wilcoxon is OK", sep = ""), file = "")
          }
           # 3- CHI & Wilcoxon are OK
          if ((global_tests[global_tests[, 1] == I, 5] != 0) & (global_tests[global_tests[, 1] == I, 2] != 0) & (global_test_choice == 3))
          {
              pass_to_local = TRUE
              write(paste("3- CHI & Wilcoxon are OK", sep = ""), file = "")
          }
          # 4- CHI OR Wilcoxon is OK
          if ((global_tests[global_tests[, 1] == I, 5] != 0) | (global_tests[global_tests[, 1] == I, 2] != 0) & (global_test_choice == 4)) {
              pass_to_local = TRUE
              write(paste("4- CHI OR Wilcoxon are OK", sep = ""), file = "")
          }
          
          if (pass_to_local)
          {
              #matrix of detected regions of interest
              regions = matrix(0, ncol = 1, nrow = nrow(gentest))
              colnames(regions) = "regions"
              
              #matrix of pval form test of each unit of regions of interest
              test = matrix("", ncol = 1, nrow = nrow(gentest))
              colnames(test) = "test"
              
              totlist = sum(as.numeric(genome.test[, "list"]))   # total number of genes of the tested list
              for (U in 1:(nrow(gentest)))
              {
                  if ((as.numeric(gentest[U, "list"]) > as.numeric(gentest[U, "hazard"])))
                  {
                    regions[U, 1] = 1
                    test[U, ] = 1
                    test[U, ] = try(as.numeric(binom.test(as.numeric(gentest[U, "list"]), totlist, as.numeric(gentest[U, "total"])/nhazard)[3]), silent=TRUE)
                  }
                  else
                  {
                    regions[U, 1] = 0
                    test[U, ] = 1
                  }
              }
              
              #pval correction for binomial test	
              write(paste("global pval correction", sep = ""), file = "")
              library(multtest)
              procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
              # correct results
              test = as.matrix(as.numeric(test))
              pval.adjusted <- mt.rawp2adjp(test, procs)
              pval.adjusted <- pval.adjusted$adjp[order(pval.adjusted$index), ]
              test = pval.adjusted[, pcorrd]
              
              gentest = cbind(gentest, regions, test)
              local_test = rbind(local_test, gentest)
  
        			#########################################################################################################
        			#
       			  #        Validation & optimization
        			#
       	 			#########################################################################################################
              
              write(paste("Validation & optimization", sep = ""), file = "")
              GExRegions = as.matrix(local_test[, "regions"])
              #region extension
              for (U in 1:nrow(GExRegions))
              {
                  if (local_test[U, "hazard"] >= local_test[U, "list"])
                  {
                    GExRegions[U, ] = 0
                  }
              }
              for (U in 2:(nrow(GExRegions) - 1))
              {
                  if ((as.numeric(local_test[U, "binome.test"]) > 0.05))
                  {
                    GExRegions[U, ] = 0
                  }
              }
              GExRegionst = GExRegions
              
              for (U in 2:(nrow(GExRegions) - 1))
              {
                  if (GExRegions[U, ] == 1)
                  {
                    GExRegionst[U - 1, ] = 1
                    GExRegionst[U + 1, ] = 1
                  }
              }
              GExRegions = GExRegionst
              if (GExRegions[1, ] == 1)
              {
                  GExRegions[2, ] = 1
              }
              if (GExRegions[nrow(GExRegions), ] == 1)
              {
                  GExRegions[(nrow(GExRegions) - 1), ] = 1
              }
              
              genome.test.tmp = cbind(local_test[, c("chr", "unit", "list", "hazard")], GExRegions)
              colnames(genome.test.tmp) = c("chr", "unit", "list", "hazard", "regions")
              #matrix for pval of all tested genome
              test = matrix(1, ncol = 1, nrow = 0)
              colnames(test) = "test"
              
              # listing data localisations of interesting regions before binding in gexROI matrix
              for (C in levels(as.data.frame(genome.test.tmp)[, "chr"]))
              {
  				      #########################################################################################################
  				      #
  				      #        censing
  				      #
  				      #########################################################################################################
  				      
                  write(paste("Finding interesting regions for chr ", C, sep = ""), file = "")
                  gexChr = genome.test.tmp[genome.test.tmp[, "chr"] == C, ]
                  start = matrix(0, ncol = 3, nrow = (nrow(gexChr) + 1))
                  start[1:nrow(start) - 1, 1] = as.numeric(gexChr[, "regions"])
                  start[2:nrow(start), 2] = as.numeric(gexChr[, "regions"])
                  
                  start[, 3] = start[, 1] - start[, 2]
                  start = start[1:nrow(start) - 1, ]
                  
                  #correspondences matrix of statrs and units
                  regions = cbind(gexChr[, "unit"], start[, 3])
                  
                  #end shifted matrix (cf to annexes)
                  end = matrix(0, ncol = 3, nrow = (nrow(gexChr) + 1))
                  end[2:nrow(end), 1] = as.numeric(gexChr[, "regions"])
                  end[1:nrow(end) - 1, 2] = as.numeric(gexChr[, "regions"])
                  
                  end[, 3] = end[, 1] - end[, 2]
                  end = end[2:nrow(end), ]
                  
                  regions = cbind(regions, end[, 3])
                  
                  gexROI = matrix(0, ncol = 3, nrow = length(regions[regions[, 2] == 1, 1]))
                  gexROI[, 1] = as.numeric(regions[as.numeric(regions[, 2]) == 1, 1])
                  gexROI[, 2] = as.numeric(regions[regions[, 3] == 1, 1])
                  
                  gexROI[, 3] = gexROI[, 2] - gexROI[, 1]
                  
                  ########################################################################################################
    				      #
    				      #        TEST
    				      #
    				      #########################################################################################################
    				      
    				      #Matrix of test results
    				      #two columns: the first for the chromosome units, the second for the pvalues
                  pval = matrix("", ncol = 2, nrow = nrow(gexChr))
                  pval[, 1] = gexChr[, "unit"]
                  
                  if(nrow(gexROI)!=0)
                  {
                    for (I in 1:nrow(gexROI))
                    {
                      #Test of entire region
                      list = as.vector(as.numeric(gexChr[gexROI[I, 1]:gexROI[I, 2], "list"]))
                      hazard = as.vector(as.numeric(gexChr[gexROI[I, 1]:gexROI[I, 2], "hazard"]))
                      p = wilcox.test(list, hazard, paired = TRUE, exact = FALSE, correct = FALSE)
                      p = as.numeric(p[3])
                      #saving results
                      pval[gexROI[I, 1]:gexROI[I, 2], 2] = p
                      
                      				         
    				      #########################################################################################################
    				      #
    				      #        PROGRESSIVE OPTIMISATION 
    				      #
    				      #########################################################################################################
    				         
    				         #localization information of the region
                      l = gexROI[I, 2] - gexROI[I, 1]
                      start = gexROI[I, 1]
                      end = gexROI[I, 2]
                      
                      #While region width>3 and pval is non significative
                      while ((p > 0.05) && (l >= 3))
                      {
    				            #upstream reduction 
                        list = as.vector(as.numeric(gexChr[(start + 1):end, "list"]))
                        hazard = as.vector(as.numeric(gexChr[(start + 1):end, "hazard"]))
                        p = wilcox.test(list, hazard, paired = TRUE, exact = FALSE, correct = FALSE)
                        p1 = as.numeric(p[3])
                        
                        #downstream reduction 
                        list = as.vector(as.numeric(gexChr[start:(end - 1), "list"]))
                        hazard = as.vector(as.numeric(gexChr[start:(end - 1), "hazard"]))
                        p2 = 1
                        p2 = try(as.numeric(wilcox.test(list, hazard, paired = TRUE, exact = FALSE, correct = FALSE)[3]), silent=TRUE)
                        
                        #Choice of best pval
                        if (p1 > p2)
                        {
                          p = p2
                          pval[start:(end - 1), 2] = p
                          end = end - 1
                        }
                        if (p2 >= p1)
                        {
                          p = p1
                          pval[(start + 1):end, 2] = p
                          start = start + 1
                        }
                        l = l - 1
                      } #End of optimisation	
                    } #End of tests for one chromosome
                    #pval correction for local wolcoxon tests
                    write(paste("local pval correction", sep = ""), file = "")
                 }else{
                  write(paste("*****No Region of Interest found\n", sep=""), file="")
                 }
                 test = rbind(test, as.matrix(pval[, 2]))
                  test = as.matrix(as.numeric(test))
                  procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
                  # correct results
                  pval.adjusted <- mt.rawp2adjp(test, procs)
                  pval.adjusted <- pval.adjusted$adjp[order(pval.adjusted$index), ]
                  test = as.matrix(pval.adjusted[, pcorrv])
              }
              colnames(GExRegions) = "Optmized"
              colnames(test) = "Wilcoxon"
              local_test = cbind(local_test, GExRegions, test)
              
              local_test[local_test[, "regions"] == "0", "regions"] = ""
              local_test[local_test[, "Optmized"] == "0", "Optmized"] = ""
              
              local_test_results = rbind(local_test_results, local_test)
          }
      }
    }
    pdf(file = file.path(paste(results.graph, "test_global.pdf", sep = "")))
    matplot(t(g_test_graph), pch = c(1:ncol(t(g_test_graph))), type = "p", col = c(1:ncol(t(g_test_graph))))
    legend(1, 0.045, pch = c(1:ncol(t(g_test_graph))), colnames(t(g_test_graph)), col = c(1:ncol(t(g_test_graph))), cex = 0.7)
    dev.off()
    write.table(g_test, row.names = TRUE, sep = "\t", file = file.path(paste(results.graph, "Chi_results.txt", sep = "")))
    write.table(global_tests, row.names = TRUE, sep = "\t", file = file.path(paste(results.graph, "global_tests_results.txt", sep = "")))
    return(local_test_results)
}
