% Part of the bcool package for R
% Distributed under the BSD license: see bcool/COPYING
% 


\name{rpoS}
\alias{rpoS}
\alias{rpoSalign}
\alias{labels}
\alias{tree7}
\title{rpoS Data from Spangenberg \emph{et al.} (2011)}
\description{
  In Spangenberg \emph{et al.} (2011) an RpoS (sigma38) alignment was scanned 
  searching for relevant columns, probably associated with pathogenicity.
  The labels (pathogenicity YES or No) were obtained from the NCBI. 
  The phylogenetic tree was reconstructed from the concatenated alignment 
  of 7 genes (K03070, K03073, K03076, K03087, K030106, K030110, K03217). 	
}
\usage{
data("rpoS")
}
\format{
  A table with the pathogenicity labels ("labels"), the phylogenetic tree 
  ("tree7") and the alignment ("rpoSalign") of the RpoS protein is provided.
}
\source{
  See below.
}
\references{
  L. Spangenberg, F. Battke, M. Grana, K. Nieselt and H. Naya (2011)
  Identifying associations between amino acid changes and meta information in alignments


}
\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}
\keyword{datasets}
