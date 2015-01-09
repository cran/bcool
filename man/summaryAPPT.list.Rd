\name{summaryAPPT.list}
\alias{summaryAPPT.list}
\title{
Summarizes the results of the runs
}
\description{
This method summarizes the relevant information of the whole run. One obtaines the calculated scores for each column and property, as well as the integrated score for each column. Moreover, the effect size is also displayed for assessment of biological relevance. In addition, chi-square scores are also caculated, and the frequency of amino acid in the categories (e. g. pathogen/non pathogen) are shown.
}
\usage{
summaryAPPT.list(x, contrast, class.var, what.prop)
}
\arguments{
  \item{x}{
APPT.list
}
  \item{contrast}{
character vector with the contrasts
}
  \item{class.var}{
character, the variable to be analyzed
}
  \item{what.prop}{
numeric vector, with the amino acid properties to be considered 
}
 
}

\value{
  \item{tabcor}{Table with dimensions |columns| x |properties|, containing the score calculated for each column and each property}
  \item{tabSize}{Table with dimensions |columns| x |properties|, containing the effect size for each column and each property}
  \item{resultado}{vector of length |column|, containing integrated score values}
  \item{resuCS}{vector of length |column|, containing chi-squared values}
  \item{liAAs}{Table containing the amino acid frequency in each of the classes, and the difference between the classes}
}
\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}


\examples{
data("rpoS")

env <- new.env()
utils::data("aaindex", package = "seqinr", envir = env)
aaindex <- env$aaindex


# select the properties to consider from "aaindex"
prop<-c("JOND920102","KYTJ820101")

# instantiate the APPT.list object
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,
properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)
colu<-c(375,376)

# define the priors and run the model
# here the number of processors is one (count=1). 
# increase this number if possible 
\dontrun{
myAPPT.list<-MCMCglmm.APPT.list(myAPPT.list, ~ -1+pathogenicity,
random.eff="spKEGG",nitt=1.5e2,burnin=5e1,prior,scale=FALSE,
parallel=TRUE,which.columns=colu,maxngaps=10,
nummin=2,count=1,pr=FALSE)

# check what levels do you want to contrast in summary
levelsMCMCglmm(myAPPT.list)
	
# summarize results and statistics
my.summary<-summaryAPPT.list(myAPPT.list,contrast=c("pathogenicityYes","pathogenicityNo"),
class.var="pathogenicity",what.prop=1:2)
}
}
