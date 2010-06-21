\name{gexmap}
\Rdversion{1.1}
\alias{gexmap}
\title{
A visual, intuitive, easy to use software giving access to a new type of information buried into your microarray data.
}
\description{
GExMap analyze genomic distribution and GO of genes lists produced by transcriptomic studies.
}
\usage{
gexmap(genome = "homosapiens", scale = "", source = "", res = "", isGO = FALSE, isMAP = TRUE, lim_chi = 5, global_test_choice = 4, pcorrd = 1, pcorrv = 1)
}
\arguments{
  \item{genome}{ Annotaions data and genomic localizations for all genes of a specific genome }
  \item{scale}{ The variable \code{scale} is the scale used to produce the graphics and to perform the statistical tests. Customizable in input in the gexmap(scale="") main function, the scale is set to 1 000 000 bp by default. }
  \item{source}{ Folder which would contain alla data source as *.Rdata files. }
  \item{res}{ Folder which will contain the results data files and folders. }
  \item{isGO}{ TRUE/FALSE, customizable input variable to perform or no the Gene Ontology analysis of the tested gene list. }
  \item{isMAP}{ TRUE/FALSE, customizable input variable to perform or no the graphical reoprts of genomic analysis in pdf files. }
  \item{lim_chi}{ Limit number of gene necessary to perform the CHI squared global test to statistically select the chromosomes. (used in the \code{gestest()} function).}
  \item{global_test_choice}{Choice of chromosomes selection according to the results of the global statistical tests (see \code{gextest()} function).\\
			1- At least CHI is OK\\
			2- At least Wilcoxon is OK\\
			3- CHI & Wilcoxon are OK\\
			4- CHI OR Wilcoxon is OK\\
			Default value is set at 4.}
  \item{pcorrd}{ Choice of pvalue correction method for the global tests (see \code{gextest()} function). }
  \item{pcorrv}{ Choice of pvalue correction method for the local tests (see \code{gextest()} function). }
}
\details{
GExMap analyze genomic distribution and GO of genes lists produced by transcriptomic studies.
Its
}
\references{
http://gexmap.voila.net/index.html
Research Article JPB/Vol.2.1 051- 059(2009)/January 2009 doi:10.4172/jpb.1000060
GExMap: An Intuitive Visual Tool to Detect and Analyze Genomic Distribution in Microarray-generated Lists of Differentially Expressed Genes
Nicolas Cagnard, Carlo Lucchesi, Gilles Chiocchia 
}
\author{
Nicolas Cagnard <nicolas.cagnard@yahoo.fr>
}
