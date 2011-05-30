\name{bootMCMCglmm.APPT.list}
\alias{bootMCMCglmm.APPT.list}
\title{
Bootstrap function for significance assessment
}
\description{
After all calculations, one might want to validate the significance of the chosen columns. A bootstrap method performed on columns selected by the user is applied permutating the observations in each iteration. In each step the scores are recalculated for each property and the integrated score is determined. Hence, an empiric random distribution of the score is determined. Score values above the bootstrap scores are significant.   
}
\usage{
bootMCMCglmm.APPT.list(x, class.var, boot, contrast, random.eff, nitt, burnin, prior, scale, parallel, what.prop, which.column, maxngaps, nummin, count, ...)
}
\arguments{
  \item{x}{
APPT.list object
}
  \item{class.var}{
character, the class considered
}
  \item{boot}{
numeric, number of iterations for the bootstrap
}
  \item{contrast}{
character, contrast considered
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
list, prior distribution for the MCMCglmm run
}
  \item{scale}{
logical, if TRUE scaling is performed in the MCMCglmm
}
  \item{parallel}{
logical, if TRUE the process runs in parallel
}
  \item{what.prop}{
numeric, which properties are considered
}
  \item{which.column}{
numeric, which columns are considered for the bootstrap
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
  \item{\dots}{
}
}

\value{
  \item{tabcor}{Table with dimensions |columns| x |properties|, containing the score calculated for each column and each property}
  \item{tabSize}{Table with dimensions |columns| x |properties|, containing the effect size for each column and each property}
  \item{resultado}{vector of length |column|, containing integrated score values}
  \item{resuCS}{vector of length |column|, containing chi-squared values}
  \item{liAAs}{Table containing the amino acid frequency in each of the classes, and the difference between the classes}
The same list as the summary.APPTlist function. Each entry might be smaller since in general the bootstrap is run only on a few potentially significant columns.
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
myAPPT.list.boot<-bootMCMCglmm.APPT.list(myAPPT.list, ~ -1+pathogenicity,boot=100,contrast=c("pathogenicityYes","pathogenicityNo"),random.eff="spKEGG",nitt=2.5e2,burnin=5e2,prior,scale=FALSE,parallel=TRUE,what.prop=c(1,2),which.column=c(376),maxngaps=10,nummin=2,count=1)
}

}




