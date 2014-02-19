### R code from vignette source 'bcool.Snw'

###################################################
### code chunk number 1: bcool.Snw:63-64
###################################################
options(width=50) 


###################################################
### code chunk number 2: bcool.Snw:90-92
###################################################
options(continue=" ")
set.seed(123456)


###################################################
### code chunk number 3: bcool.Snw:95-96
###################################################
suppressMessages(require("bcool"))


###################################################
### code chunk number 4: bcool.Snw:99-107
###################################################
library("bcool")
data("rpoS")

env <- new.env()
utils::data("aaindex", package = "seqinr", envir = env)
aaindex <- env$aaindex




###################################################
### code chunk number 5: bcool.Snw:118-121
###################################################
options(width=50) 
prop<-c("CHOC760102","KYTJ820101")
aaindex[prop]


###################################################
### code chunk number 6: bcool.Snw:128-135
###################################################
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,
properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)
head(columns(myAPPT.list))
head(pheno(myAPPT.list))
properties(myAPPT.list)
tree(myAPPT.list)



###################################################
### code chunk number 7: bcool.Snw:143-154
###################################################
my.cer<-cer.APPT.list(myAPPT.list,class.var="pathogenicity",
which.columns=NULL,nummin=2,maxngaps=10)
head(sort(my.cer))
tail(sort(my.cer))

my.anova<-anovaAPPT.list(myAPPT.list,class.var=~pathogenicity,
which.columns=NULL,nummin=2,maxngaps=10)
sum(my.anova$hm.signif==4,na.rm=TRUE)
which(my.anova$hm.signif==4)
head(my.anova)



###################################################
### code chunk number 8: bcool.Snw:160-172
###################################################
colu<-c(374:376)

# define the priors and run the model
# here the number of processors is one (count=1). 
# increase this number if possible 
prior<-list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))
myAPPT.list<-MCMCglmm.APPT.list(myAPPT.list, ~ -1+pathogenicity,
	random.eff="spKEGG",nitt=1.5e3,burnin=5e2,prior,scale=FALSE,
	parallel=TRUE,which.columns=colu,maxngaps=10,nummin=2,
	count=1,pr=FALSE)




###################################################
### code chunk number 9: bcool.Snw:182-198
###################################################
matGwk<-matrix(0,length(columns(myAPPT.list)),
length(properties(myAPPT.list)))
for (i in 1:length(properties(myAPPT.list))){
	for (j in 1:length(columns(myAPPT.list))){
		matGwk[j,i]<-geweke.diag(as.mcmc(
		multiMCMCglmm(myAPPT.list)[[i]][[j]]$Sol[,2]
		-multiMCMCglmm(myAPPT.list)[[i]][[j]]$Sol[,1]
		))$z
	}
}
colnames(matGwk)<-names(properties(myAPPT.list))
rownames(matGwk)<-columns(myAPPT.list)

matGwk
which(abs(matGwk)>2,TRUE)



###################################################
### code chunk number 10: bcool.Snw:204-217
###################################################
levelsMCMCglmm(myAPPT.list)
my.summary<-summaryAPPT.list(myAPPT.list,
contrast=c("pathogenicityYes","pathogenicityNo"),
class.var="pathogenicity",what.prop=NULL)
attributes(my.summary)

my.summary$tabcor
my.summary$tabSize
my.summary$SumTr
my.summary$ChiSq
my.summary$ChiSq.eff
my.summary$AAlist



###################################################
### code chunk number 11: bcool.Snw:226-237
###################################################
myAPPT.list.boot<-bootMCMCglmm.APPT.list(
	myAPPT.list, ~ -1+pathogenicity,boot=100,
	contrast=c("pathogenicityYes","pathogenicityNo"),
	random.eff="spKEGG",nitt=1.5e3,burnin=5e2,prior,
	scale=FALSE,parallel=TRUE,what.prop=c(1,2),
	which.column=c(376),maxngaps=10,nummin=2,count=1)

boxplot(myAPPT.list.boot$SumTr,col="orange",notch=TRUE,ylim=c(0,1))
par(new=TRUE)
plot(my.summary$SumTr["376"],ylim=c(0,1),pch=17,col="red",xaxt="n",yaxt="n",cex=2,ylab="SumTr")



###################################################
### code chunk number 12: bcool.Snw:246-261
###################################################
matH2<-matrix(0,length(columns(myAPPT.list)),
length(properties(myAPPT.list)))
for (i in 1:length(properties(myAPPT.list))){
	for (j in 1:length(columns(myAPPT.list))){
		matH2[j,i]<-median(
		multiMCMCglmm(myAPPT.list)[[i]][[j]]$VCV[,1]
		/(multiMCMCglmm(myAPPT.list)[[i]][[j]]$VCV[,1]+
		multiMCMCglmm(myAPPT.list)[[i]][[j]]$VCV[,2]))
	}
}
colnames(matH2)<-names(properties(myAPPT.list))
rownames(matH2)<-columns(myAPPT.list)

matH2



