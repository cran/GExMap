\name{gextest}
\alias{gextest}
\title{ Perform the statistical analysis step }
\description{
  This function selects the chromosoms by two global tests comparing an estimated hazard genomic distribution to the real distribution extracted from the user list.
  If the chromosom is selected, a local analysis is performed to localize and test regions of interest .
}
\usage{
gextest(genome.test, nhazard, results.graph, lim_chi, global_test_choice, pcorrd, pcorrv)
}
\arguments{
  \item{genome.test}{ All genomic data to be tested. }
  \item{nhazard}{ Total number of genes of the hazard distribution. }
  \item{results.graph}{ Data necessary to produces the graphic files. }
  \item{lim_chi}{ Limit of genes to perform the CHI squared global test (set to 5 by default). }
  \item{global_test_choice}{Choice of chromosomes selection according to the results of the global statistical tests (see \code{gextest()} function).\\
			1- At least CHI is OK\\
			2- At least Wilcoxon is OK\\
			3- CHI & Wilcoxon are OK\\
			4- CHI OR Wilcoxon is OK\\
			Default value is set at 4.}
  \item{pcorrd}{ Choice of pvalue correction method for the global tests:\\
     1- "Bonferroni"
     2- "Holm"
     3- "Hochberg"
     4- "SidakSS"
     5- "SidakSD"
     6- "BH"
     7- "BY"
  }
  \item{pcorrv}{ Choice of pvalue correction method for the local tests:\\
     1- "Bonferroni"
     2- "Holm"
     3- "Hochberg"
     4- "SidakSS"
     5- "SidakSD"
     6- "BH"
     7- "BY"
  }
}
\details{
  This function is included in the GExMap package and could not be used independently.
}
\value{
}
\references{  }
\author{ Nicolas Cagnard <nicolas.cagnard@inserm.fr> }
\note{ 
}
\seealso{GExMap }
\keyword{ misc }

