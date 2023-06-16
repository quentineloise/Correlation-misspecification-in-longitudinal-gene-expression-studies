library(nlme)
library(reshape2)
library(Hmisc)
library(GEOquery)
library(ggplot2)

# Obtaining the raw data set from the GEOquery library

gse30550<-getGEO('GSE30550',GSEMatrix=TRUE)
des<-pData(phenoData(gse30550[[1]]))
eset<-exprs(gse30550[[1]])


#Cleaning the des (designFlu) dataset QE

sub<-substr(des$title,9,10)  # generate sub column with only the 9 to 10 characters of the orginal data column
time<-substr(des$characteristics_ch1.1,15,18) # generate time column with only the 17 to 18 characters of the original data column
cond<-substr(des$characteristics_ch1.3,15,19) # generate cond column with the 17 to 18 characters of the original data column
Array<-1:268 # generate array column with an index from 1:268
columnname<-rownames(des)
#time=as.numeric(time)
time2=time
des<-data.frame(sub,time,cond,Array,columnname,time2) # create the new des data set
des[,2]<- sub("line","-12",as.character(des[,2])) # replace ne with -1 (representing the baseline )
des[,3]<- sub("t","",as.character(des[,3])) #delete t from column 
des$time2<- sub("line","-12",as.character(des$time2)) #delete t from column 
des$time=as.numeric(des$time)
des$time2=as.numeric(des$time2)


#it is important before running any code to check to make sure that the following are true
#correct it if it is not.
#
# 1.  To be safe make sure that eset is a matrix of numbers and that 
#  the matrix has row name id's.  Sometimes the first column will have the row ids.
# 
# 2.  After adressing one, make sure that the expression values are on the log scale.
#    Just check if there are any values in the hundreds and thousands.  If so it needs
#    to be log2 transformed.
#
# 3. Make sure the colnames of eset match up with the des$columnname in des.
# 4. Make sure that cond,time, and sub are factors in des.  Time2 should be numeric.

#Checking #1
head(eset)
#first column is the row id's so lets fix it.
#rownames(eset)<-eset[,1]
#eset<-as.matrix(eset[,-1])


#Checking #2
summary(as.vector(eset))
#The values are already log2 transformed since the max is only 15.11. If it were not
#it would be something like 2^15=32768
#If you needed to fix simply overwrite eset,  eset<-log2(eset) or log2(eset-min(eset))
#if there are negative numbers.

#Checking #3
length(setdiff(colnames(eset),as.character(des$columnname)))
#This should be 0, if not, then there are cols in eset that do not match 
#in design or vice versa.  Need to fix by subsetting the two files so that they match.

#Checking #4
is.factor(des$cond)
is.factor(des$time)
is.factor(des$sub)
is.numeric(des$time2)

#fixing variable types
table(des$time)
des$time<-factor(des$time)
des$sub<-factor(des$sub)
des$integer<-as.numeric(factor(des$time))

#Sorting desing based on subject and time point.
#First ordering everything correctly.
#This is critical to do.
index<-order(des$columnname)
des<-des[index,]
index2<-order(colnames(eset))
eset<-eset[,index2]

#Ordering by subject and time point now
index3<-order(des$sub,des$time)
des<-des[index3,]
eset<-eset[,index3] # ?? NOT RUNNING NOT SURE WHAT IS THAT DOING


#Initializing result to store gene level correlation estimates as a function of d_jk
f.hat.final<-c()

for(i in 1:11961){
  des$y<-eset[i,]
  r.raw<- residuals(lm(y~ time+cond+time:cond, data = des))
  r.frame<-data.frame(sub=des$sub,time=des$time,resids=r.raw)
  
  
  r.wide<-acast(r.frame,sub~time,value.var="resids")
  r.cor<-rcorr(r.wide)$r
  r.tall<-data.frame(rows=rownames(r.cor)[row(r.cor)], vars=colnames(r.cor)[col(r.cor)],
                     values=c(r.cor))
  
  r.tall$d<-abs(as.numeric(as.character(r.tall$rows))-as.numeric(as.character(r.tall$vars)))
  
  f.hat<-aggregate(values~d,data=r.tall,mean)
  f.hat.final<-rbind(f.hat.final,f.hat)
}

#Computing emperical correlation estimates across all genes
f.bar<-aggregate(values~d,data=f.hat.final,mean)


#Plotting Figure 7  Emperical correlation function of the flu data set
f.bar<-f.bar[-1,]  #Removing d_jk=0
#graph without the curve
fig.7a<-ggplot(f.bar,aes(x=d,y=values))+geom_point()+labs(x=expression(d[jk]),y=expression(f(d[jk])))+
  ylim(0,0.75)+theme_bw()+theme(plot.title = element_text(hjust = 0.5))

#graph with the curve
fig.7b<-ggplot(f.bar,aes(x=d,y=values))+geom_point()+labs(x=expression(d[jk]),y=expression(f(d[jk])))+geom_smooth(span=1.5)+
  ylim(0,0.75)+theme_bw()+theme(plot.title = element_text(hjust = 0.5))



png("Figure7.png", width = 5, height = 3.75, units="in", res=2400)
#fig.7a
fig.7b
dev.off()






#Creating Figure 8.  Fitting sphercial correlation function to the empirical estimates 
#via nonlinear least squares.


#Obtaining least squares estimates using the spherical correlation function
f1<- values~ifelse(rho1>d,(1-rho2-nugget)*(1-1.5*d/rho1+0.5*(d/rho1)^3)+rho2,rho2) 
spherfit<-nls(f1,data=f.bar,start=list(rho1=30,rho2=0.3,nugget=0.6))
coef(spherfit)


#Creating modified spherical correlation function with parameters rho1, rho2, and nugget
spher<-function(d,rho1,rho2,nugget){
  result<-c()
  for (i in 1:length(d)){
    result[i]<-ifelse(rho1>d[i],(1-rho2-nugget)*(1-1.5*d[i]/rho1+0.5*(d[i]/rho1)^3)+rho2,rho2)
  }
  return(result)
}



png("Figure8.png", width = 5, height = 3.75, units="in", res=2400)
init.start<-list(rho=f.bar$values[1]^(1/5))
arfit<-nls(values~ rho^(d),data=f.bar,start=init.start)
csfit<-lm(values~ 1,data=f.bar)
plot(f.bar$d,f.bar$values,ylim=c(0,.75),pch=1,cex=0.4,ylab=expression(Correlation~~f(d[jk])),xlab=expression(d[jk]))
index<-0:120
lines(index,spher(index,rho1=10.4713871 ,rho2=0.32198904,nugget=0.2548449),lwd=3,lty="dotted")
lines(index,0.9659633^(index),lwd=3,lty="twodash")
abline(h=coef(csfit),lwd=3,lty="solid")
legend("topright", legend=c("CS", "AR","SPH"),lty=c("solid","twodash", "dotted"),cex=1,text.font=4,lwd=3)
dev.off()

#MSE estimates and AIC
summary(spherfit)$sigma^2
summary(arfit)$sigma^2
summary(csfit)$sigma^2

AIC(spherfit)
AIC(arfit)
AIC(csfit)

#############################################
#Table 1
#
#Applying the spherical fit to perform a genearlized linear model to determine difference over 
#time for Asymptomatic and Symptomatic groups
#
#

#Removing the -12 basline values 
index4<-which(des$time2==-12)
des<-des[-index4,]
eset<-eset[,-index4]
des$sub<-factor(des$sub)
des$integer<-as.numeric(factor(des$time))



djk<-as.vector(dist(unique(des$time2)))
fit.vals<-spher(djk,rho1=10.7705836 ,rho2=0.3107967,nugget=0.2532736)
custom.spher <- corSymm(value =fit.vals,
                        form = ~ integer | sub, fixed=TRUE)
custom.spher.init <- Initialize(custom.spher, data = des)
#corMatrix(custom.spher.init)  This line allows the user to see each subjects Correlation matrix as specified by the spherical correlation fit.




#On rare occasions, the gls model (using REML) may fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.

#The following loop performs the the three generalized least squares models using the SPH, AR, and CS
#structures.  The parameters for each correlation function were determined by the nonlinear least squares
#fits obtained in the code above when producing Figure 8.

#Since there are over 30 regression coefficients, we elected to not write specific contrasts to obtain the t-statistics provided by summary can obtain the t-statistics that correspond
#and p-values for the tests for mean changes (HrX vs Baseline) for each condition (symptomatic/asymptomatic). 
#To do this, we effectively run the model twice specifying the reference group for intercept to correspond to 
#the baseline symptomatic population for the first model, and the baseline asymptomatic population for the second.
#The coefficients corresponding to the time variable (no interaction) correspond to the HrX vs Baseline comparison for the reference group
#and their p-values are readily obtained using summary.


#Initializing Simulation Result objects
result.spher<-c()    #Object to store p-values for comparing HrX -Hr0 for Symp and Asymp groups using Spher
result.ar<-c()       #Object to store p-values for comparing HrX -Hr0 for Symp and Asymp groups using AR
result.cs<-c()

for (i in 1:11961){
  #des<-data.frame(time=rep(n.time,n.sub),sub=paste("S",rep(1:n.sub,each=length(n.time)),sep=""))
  #des$time2<-des$time
  #des$time<-factor(paste("T",des$time,sep=""))
  des$y<-eset[i,]
  #eset<-rbind(eset,des$y)
  des$cond<-factor(as.character(des$cond),levels=c("Symp","Asymp"))

  
  tryCatch({
    spher<-gls(model = y~ time+cond+time:cond, data = des, correlation = custom.spher.init,control = glsControl(opt = "optim")) 
    ar<-gls(model = y~ time+cond+time:cond, data = des, correlation = corExp(value=-1/log(coef(arfit)),form = ~time2|sub,fixed=T),control = glsControl(opt = "optim")) #value=-1/log(coef(arfit)),
    cs<-gls(model = y~ time+cond+time:cond, data = des, correlation = corCompSymm(value=coef(csfit),form = ~1|sub,fixed=T),control = glsControl(opt = "optim")) #value=coef(csfit),
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  des$cond<-factor(as.character(des$cond),levels=c("Asymp","Symp"))
  tryCatch({
    spher2<-gls(model = y~ time+cond+time:cond, data = des, correlation = custom.spher.init,control = glsControl(opt = "optim")) 
    ar2<-gls(model = y~ time+cond+time:cond, data = des, correlation = corExp(value=-1/log(coef(arfit)),form = ~time2|sub,fixed=T),control = glsControl(opt = "optim")) #value=-1/log(coef(arfit)),
    cs2<-gls(model = y~ time+cond+time:cond, data = des, correlation = corCompSymm(value=coef(csfit),form = ~1|sub,fixed=T),control = glsControl(opt = "optim")) #value=coef(csfit),
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  
  #Combining the HrX vs Hr0 p-value comparisons for Symptomatic (spher and ar) with the same comparisons for Asymptomatic (spher2 and ar2)
  result.spher<-rbind(result.spher,c(coef(summary(spher))[2:15,4],coef(summary(spher2))[2:15,4]))
  result.ar<-rbind(result.ar,c(coef(summary(ar))[2:15,4],coef(summary(ar2))[2:15,4]))
  result.cs<-rbind(result.cs,c(coef(summary(cs))[2:15,4],coef(summary(cs2))[2:15,4]))
  
}

adjp.spher<-apply(result.spher,2,p.adjust,method="fdr")
adjp.ar<-apply(result.ar,2,p.adjust,method="fdr")
adjp.cs<-apply(result.cs,2,p.adjust,method="fdr")

rejected.spher<-apply(adjp.spher<.1,2,sum)
rejected.ar<-apply(adjp.ar<.1,2,sum)
rejected.cs<-apply(adjp.cs<.1,2,sum)

cbind(rejected.spher,rejected.cs,rejected.ar)

table.1<-data.frame(Comparison=paste("HR",unique(des$time2)[-(1:1)],"-Base",sep=""),
                    SympSPHER=rejected.spher[1:15],
                    SympAR=rejected.ar[1:15],
                    SympCS=rejected.cs[1:15],
                    AsympSPH=rejected.spher[16:30],
                    AsympAR=rejected.ar[16:30],
                    AsympCS=rejected.cs[16:30])

View(table.1)

write.csv(table.1,file="Table1.csv")


#For curious readers and reviewers.  The significant genes in the Spher model are almost always contined within
#the gene lists of the AR model.  The Spher model just has additional significant genes.
counts.spher<- adjp.spher<.1
counts.ar<- adjp.ar<.1

common<-c()
for (i in 1:26){
  common[i]<-length(intersect(which(counts.spher[,i]==TRUE), which(counts.ar[,i]==TRUE)))
}

CommonGenes<-data.frame(Comparison=paste("HR",unique(des$time2)[-(1:2)],"-Base",sep=""),
                        Symp=common[1:14],
                        Asymp=common[15:28])

View(CommonGenes)


write.csv(CommonGenes,file="CommonGenes.csv")












