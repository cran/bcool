\name{RGBcolor.APPT.list}
\alias{RGBcolor.APPT.list}
\title{
Visualization method using the RGB color scale
}
\description{
Visualization method. Converts the amino acid properties into colors, using the RGB scale
}
\usage{
RGBcolor.APPT.list(x, which.prop, write.score, ...)
}
\arguments{
  \item{x}{
  	APPT.list object
}
  \item{which.prop}{
  	numeric, the amino acid properties to be considered
}
  \item{write.score}{
  	logical, if TRUE the scores are written into a file
}
  \item{\dots}{
}
}

\value{
A matrix containing colors according to the RGB scale
}

\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}

\examples{
data("rpoS")
data("aaindex")

# select the properties to consider from "aaindex"
prop<-c("JOND920102","KYTJ820101")

# instantiate the APPT.list object
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)
colu<-c(375,376)

# define the priors and run the model
# here the number of processors is one (count=1). 
# increase this number if possible 
\dontrun{
prior<-list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))
myAPPT.list<-MCMCglmm.APPT.list(myAPPT.list, ~ -1+pathogenicity,random.eff="spKEGG",nitt=1.5e2,burnin=5e1,prior,scale=FALSE,parallel=TRUE,which.columns=colu,maxngaps=10,nummin=2,count=1,pr=FALSE)

# function that return a matrix of RGB colors for the columns in the alignment
matRGB<-RGBcolor.APPT.list(myAPPT.list,which.prop=0,FALSE)
matRGB
}
}

