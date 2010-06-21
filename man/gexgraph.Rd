\name{gexgraph}
\Rdversion{1.1}
\alias{gexgraph}
\title{
Produces graphics to illustrate genomic density and statistical test results
}
\description{
Function of the GExMap package which produces graphics in pdf files with data produces by the gexmap() main function.
}
\usage{
gexgraph(genome.test, scale, results.graph, W, cytobands)
}
\arguments{
  \item{genome.test}{ Matrix of all data of genomic density, localizations and tests results to be graphically reported. }
  \item{scale}{ The variable \code{scale} is the scale used to produce the graphics and to perform the statistical tests. Customizable in input in the gexmap(scale="") main function, the scale is set to 1 000 000 bp by default. }
  \item{results.graph}{ Folder path where will be saved all the graphics pdf files. }
  \item{W}{  \code{W} is the variable which permit to produecs graphics for each of the two methods of hazard estimation. 1 for the specific method, 2 for the global method. }
  \item{cytobands}{ Matrix of genomic localization to draw cytobands }
}
\details{
This function is included in the GExMap package and could not be used independently.
}
\value{
This function produces a map in a pdf file by chromosom statistically validated.
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
