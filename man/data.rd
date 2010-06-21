\name{data}
\Rdversion{1.1}
\alias{data}
\docType{data}
\title{
Genomic data set used as an example by GExMap
}
\description{
Text file tab separated with two columns.
 The first columns is for the probes identifiers  and the second for the variation code (1 or -1).
 The title of the first column must be formated as follow: [type of mocroarray], [microarray name].
 The title of the second column must be "expression".
}
\usage{data(data)}
\format{
Rdata file which contains a matrix data with all informations and annotations about a genome needed by GExMap.
}
\details{
All the data sources have been taked from the Enemble website (biomart).\\
The matrix columns are:\\ 
"chr": the chromosome names.\\
"cytoband": Cytoband names.\\
"ensembl": ensembl IDs.\\
"name": Genes names.\\
"m": The middle of the cytoband genomic localization.\\
"GO": All the GO terms about the gene.\\
}
\source{
http://www.ensembl.org/biomart/martview/
}
\references{
http://gexmap.voila.net/index.html
}
\examples{
data(data)
}
\keyword{datasets}
