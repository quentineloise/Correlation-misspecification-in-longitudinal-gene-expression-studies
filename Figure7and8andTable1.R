###########################################################################

#(TO ONLY RUN ONCE) To install the GEOquery and also install BiocManager 
# and Biobase required for GEOquery (if not installed yet) 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("GEOquery")
###########################################################################

library(nlme)
library(reshape2)
library(Hmisc)
library(GEOquery)
library(ggplot2)

# Obtaining the raw data set from the GEOquery library
gse30550<-getGEO('GSE30550',GSEMatrix=TRUE)
des<-pData(phenoData(gse30550[[1]]))
eset<-exprs(gse30550[[1]])


#Extracting key sample info from experimental design file (des) dataset 
sub<-substr(des$title,9,10)  # generate subject column with only the 9 to 10 characters of the original data column
time<-substr(des$characteristics_ch1.1,15,18) # generate time column with only the 17 to 18 characters of the original data column
cond<-substr(des$characteristics_ch1.3,15,19) # generate symptom condition column with the 17 to 18 characters of the original data column
Array<-1:268 # generate array column with an index from 1:268
columnname<-rownames(des)
time2=time

#Creating new "clean" design file (des)
des<-data.frame(sub,time,cond,Array,columnname,time2) # create the new des data set
des[,2]<- sub("line","-12",as.character(des[,2])) # replace ne with -1 (representing the baseline )
des[,3]<- sub("t","",as.character(des[,3])) #delete t from column 
des$time2<- sub("line","-12",as.character(des$time2)) #delete t from column 
des$time=as.numeric(des$time)
des$time2=as.numeric(des$time2)


#It is important before analysis to check to make sure that the following are true
# and correct it if it is not.
#
# 1.  To be safe make sure that eset is a matrix of numbers and that 
#  the matrix has row name id's.  Sometimes the first column will have the row ids.
# 
# 2.  After addressing #1, make sure that the expression values are on the log scale.
#    Just check if there are any values in the thousands.  If so it needs
#    to be log2 transformed. Check documentation from GEO as well.
#
# 3. Make sure the colnames of eset match up with the des$columnname in des.
#
# 4. The order of the columns in eset should match the order of rows in design

#Checking #1
head(eset)  #eset is in the proper format

#Checking #2
summary(as.vector(eset))
#The values are already log2 transformed since the max is only 15.11. If it were not
#it would be something like 2^15=32768
#If you needed to fix simply overwrite eset,  eset<-log2(eset) or log2(eset-min(eset))
#if there are negative numbers.

#Checking #3
length(setdiff(colnames(eset),as.character(des$columnname)))
#This should be 0, if not, then there are cols in eset that do not match 
#in the design file or vice versa.  Need to fix by subsetting the two files so that they match.

#Checking #4.  This data set is already ordered correctly. 
head(cbind(colnames(eset),des$columnname))
cbind(colnames(eset),des$columnname)



#Final data cleaning processes
#1. Reorder design and eset so that they are ordered by subject and time point necessary
#   for computation of f.hat and f.bar (the correlogram visualization proposed)
#2. Convert time and subject to factors for linear modeling purposes.

#Addressing 1: (Samples are already ordered correctly but we do it here anyways for illustration)
index2<-order(des$sub,des$time)
des<-des[index2,]
eset<-eset[,index2]

#Addressing 2
des$time<-factor(des$time)
des$sub<-factor(des$sub)





#Computing correlation estimates and graphic
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
#graph without the loess curve (not in paper)
fig.7a<-ggplot(f.bar,aes(x=d,y=values))+geom_point()+labs(x=expression(d[jk]),y=expression(f(d[jk])))+
  ylim(0,0.75)+theme_bw()+theme(plot.title = element_text(hjust = 0.5))

#graph with loess curve
fig.7b<-ggplot(f.bar,aes(x=d,y=values))+geom_point()+labs(x=expression(d[jk]),y=expression(f(d[jk])))+geom_smooth(span=1.5)+
  ylim(0,0.75)+theme_bw()+theme(plot.title = element_text(hjust = 0.5))


png("Figure7.png", width = 5, height = 3.75, units="in", res=2400)
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
#Nonlinear fit for AR parameter estimate
arfit<-nls(values~ rho^(d),data=f.bar,start=init.start)
#Using intercept only model for CS parameter estimate
csfit<-lm(values~ 1,data=f.bar)
#Plotting the points
plot(f.bar$d,f.bar$values,ylim=c(0,.75),pch=1,cex=0.4,ylab=expression(Correlation~~f(d[jk])),xlab=expression(d[jk]))
index<-0:120
#Overlaying correlation function fits
lines(index,spher(index,rho1=10.7705836 ,rho2=0.3107967,nugget=.2532736),lwd=3,lty="dotted")
lines(index,0.9659633^(index),lwd=3,lty="twodash")
abline(h=coef(csfit),lwd=3,lty="solid")
#Adding legend
legend("topright", legend=c("CS", "AR","SPH"),lty=c("solid","twodash", "dotted"),cex=1,text.font=4,lwd=3)
dev.off()


#MSE estimates and AIC comparing correlation structure fits
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
#Since for this demonstration we want to compare each time point to baseline
#we will drop the data associated with HR12 so that the intercept of the model is associated with HR0
#This allows us to easily gather t-test and p-values for each comparison using summary() rather than
#writing individaul contrasts for the many many comparisons explored.

#Removing the -12 baseline values 
index3<-which(des$time2==-12)
des<-des[-index3,]
eset<-eset[,-index3]
des$sub<-factor(des$sub)
des$integer<-as.numeric(factor(des$time))



#Creating custom correlation function estimate to incorporate inside
#of gls()
djk<-as.vector(dist(unique(des$time2)))
fit.vals<-spher(djk,rho1=10.7705836 ,rho2=0.3107967,nugget=0.2532736)
custom.spher <- corSymm(value =fit.vals,
                        form = ~ integer | sub, fixed=TRUE)
custom.spher.init <- Initialize(custom.spher, data = des)




#On rare occasions, the gls model (using REML) may fail to converge, typically when fitting AR 
#or GAUS structures. We provide additional code to note them with a warning when they occur and move forward in the
#simulation loop.  We did not experience any issues for this particular data set.

#The following loop performs the  three generalized least squares models using the SPH, AR, and CS
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
  #Adding gene i to design file for modeling
  des$y<-eset[i,]
  
  #Making Symp group the reference group
  des$cond<-factor(as.character(des$cond),levels=c("Symp","Asymp"))
  
  tryCatch({
    spher<-gls(model = y~ time+cond+time:cond, data = des, correlation = custom.spher.init,control = glsControl(opt = "optim")) 
    ar<-gls(model = y~ time+cond+time:cond, data = des, correlation = corExp(value=-1/log(coef(arfit)),form = ~time2|sub,fixed=T),control = glsControl(opt = "optim")) #value=-1/log(coef(arfit)),
    cs<-gls(model = y~ time+cond+time:cond, data = des, correlation = corCompSymm(value=coef(csfit),form = ~1|sub,fixed=T),control = glsControl(opt = "optim")) #value=coef(csfit),
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  #Making Asymp group the reference group
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

#Computing FRD adjusted p-values for each comparison
adjp.spher<-apply(result.spher,2,p.adjust,method="fdr")
adjp.ar<-apply(result.ar,2,p.adjust,method="fdr")
adjp.cs<-apply(result.cs,2,p.adjust,method="fdr")

#Computing the number of rejected test using FDR adjustment (p-value<.1)
rejected.spher<-apply(adjp.spher<.1,2,sum)
rejected.ar<-apply(adjp.ar<.1,2,sum)
rejected.cs<-apply(adjp.cs<.1,2,sum)

#Creating final result table in paper
table.1<-data.frame(Comparison=paste("HR",unique(des$time2)[-(1:1)],"-Base",sep=""),
                    SympSPHER=rejected.spher[1:14],
                    SympAR=rejected.ar[1:14],
                    SympCS=rejected.cs[1:14],
                    AsympSPH=rejected.spher[15:28],
                    AsympAR=rejected.ar[15:28],
                    AsympCS=rejected.cs[15:28])

View(table.1)
write.csv(table.1,file="Table1.csv")












