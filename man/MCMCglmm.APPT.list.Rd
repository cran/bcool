\name{MCMCglmm.APPT.list}
\alias{MCMCglmm.APPT.list}
\title{
Calculates the MCMCglmm for the APPT.list object
}
\description{
Calculates the MCMCglmm for the APPT.list object. The MCMCglmm is run for each column, integrated in a list of MCMCglmm and stored in the correspondent slot MCMCglmm.list of the APPT.list object
}


\usage{
MCMCglmm.APPT.list(x, class.var,random.eff,nitt=5.5e4,burnin=5e3,prior,scale=FALSE,parallel=FALSE,maxngaps=0,nummin=2,count=4,pr=FALSE,which.columns=NULL)

}
\arguments{
  \item{x}{
  	APPT.list object

}

\item{class.var}{
character, the class variable considered
}
  \item{random.eff}{
character, random effects considered
}
  \item{nitt}{
numeric, number of iterations for the MCMCglmm runs
}
  \item{burnin}{
numeric, burnin for the MCMCglmm run
}
  \item{prior}{
list, prior distributions for the MCMCglmm run
}
  \item{scale}{
logical, if TRUE scaling is performed in the MCMCglmm
}
  \item{parallel}{
logical, if TRUE the process runs in parallel
}

  \item{maxngaps}{
numeric, maximal number of gaps in oder to apply calculations to the column. If gaps>maxngaps, the column is not considered.
}
  \item{nummin}{
numeric, minimal number of properties considered
}
  \item{count}{
numeric, number of cores to run in parallel. The parallel argument must be set to TRUE.
}
  \item{pr}{
logical, if TRUE save random effects
}

  \item{which.columns}{
numeric, which columns are considered for the bootstrap
}

}
\value{
The APPT.list object with the slot MCMCglmm.list filled
}

\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}

\seealso{
  \code{\link{MCMCglmm}}
}

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
}
}

