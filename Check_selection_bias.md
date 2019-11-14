---
title: "Selection_bias_check"
author: "Cancan QI"
date: "11/14/2019"
output:
  html_document: 
    keep_md: true
  pdf_document: default
---



## R Markdown
load library

```r
library(foreign); library(data.table);library(MASS) 
```

```
## Warning: package 'foreign' was built under R version 3.5.2
```

```
## Warning: package 'data.table' was built under R version 3.5.2
```

```
## Warning: package 'MASS' was built under R version 3.5.2
```

```r
library(sandwich) ;library(lmtest) ;library(parallel) 
```

```
## Warning: package 'sandwich' was built under R version 3.5.2
```

```
## Warning: package 'lmtest' was built under R version 3.5.2
```

```
## Loading required package: zoo
```

```
## Warning: package 'zoo' was built under R version 3.5.2
```

```
## 
## Attaching package: 'zoo'
```

```
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```r
library(R.utils); library(matrixStats) ;library(plyr)
```

```
## Warning: package 'R.utils' was built under R version 3.5.2
```

```
## Loading required package: R.oo
```

```
## Loading required package: R.methodsS3
```

```
## R.methodsS3 v1.7.1 (2016-02-15) successfully loaded. See ?R.methodsS3 for help.
```

```
## R.oo v1.22.0 (2018-04-21) successfully loaded. See ?R.oo for help.
```

```
## 
## Attaching package: 'R.oo'
```

```
## The following objects are masked from 'package:methods':
## 
##     getClasses, getMethods
```

```
## The following objects are masked from 'package:base':
## 
##     attach, detach, gc, load, save
```

```
## R.utils v2.9.0 successfully loaded. See ?R.utils for help.
```

```
## 
## Attaching package: 'R.utils'
```

```
## The following object is masked from 'package:lmtest':
## 
##     reset
```

```
## The following object is masked from 'package:utils':
## 
##     timestamp
```

```
## The following objects are masked from 'package:base':
## 
##     cat, commandArgs, getOption, inherits, isOpen, parse, warnings
```

```
## 
## Attaching package: 'plyr'
```

```
## The following object is masked from 'package:matrixStats':
## 
##     count
```

```r
library(markdown)
```

```
## Warning: package 'markdown' was built under R version 3.5.2
```

load phenptype

```r
load("~/Documents/Projects/PIAMA/Allergic diseases/phenotype_origin/Phenotype_selection_bias.Rdata")
pheno16<-read.csv("~/Documents/Projects/PIAMA/Allergic diseases//phenotype_origin/Phenotype_brush16_PIAMA_add_allergy.csv")
df455$sample<-pheno16$Sample[match(df455$ID,pheno16$newid)]
```

load beta values of top CpGs

```r
load("~/Documents/Projects/PIAMA/Allergic diseases//data/trimedM_top_cpgs.Rdata")
load("~/Documents/Projects/PIAMA/Allergic diseases//data/cpgname_replication_noNA.Rdata")
cpg.list<-unique(c(allergy.rep,asthma.rep,rhinitis.rep))
length(cpg.list)
```

```
## [1] 68
```

**1.association between maternal educaiton and DNA methylation **

* make phenotype of maternal education

```r
pheno.ses<-df455[,c("sample","educmother","age16yrs","sex","COORD")]
pheno.ses$educmother<-as.factor(pheno.ses$educmother)
levels(pheno.ses$educmother)<-c("0","0","1")
pheno.ses<-na.omit(pheno.ses)
dim(pheno.ses)
```

```
## [1] 455   5
```
* match m values

```r
M_matrix<-t(tM.top[match(cpg.list,rownames(tM.top)),match(pheno.ses$sample,colnames(tM.top))])
```
* run logistic regression model

```r
GLMtest = function(meth_matrix,methcol,Y) {
  mod = glm(Y~meth_matrix[, methcol],family = "binomial")
  cf = summary(mod)$coefficients
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]  
}

system.time(ind.res <- mclapply(setNames(seq_len(ncol(M_matrix)), dimnames(M_matrix)[[2]]), GLMtest, meth_matrix=M_matrix, Y=pheno.ses[,2]))
```

```
##    user  system elapsed 
##   0.003   0.002   0.176
```

```r
all.results<-ldply(ind.res,rbind)
names(all.results)<-c("probeID","BETA","SE", "P_VAL")
```
* show the significant results

```r
all.results[which(all.results$P_VAL<0.05),]
```

```
## [1] probeID BETA    SE      P_VAL  
## <0 rows> (or 0-length row.names)
```

**2.association between maternal smoking and DNA methylation **

* make phenotype of maternal education

```r
pheno.ses<-df455[,c("sample","ampgsm","age16yrs","sex","COORD")]
pheno.ses$ampgsm<-as.factor(pheno.ses$ampgsm)
pheno.ses<-na.omit(pheno.ses)
```
* match M values

```r
M_matrix<-t(tM.top[match(cpg.list,rownames(tM.top)),match(pheno.ses$sample,colnames(tM.top))])
```
* run logistic regression model

```r
GLMtest = function(meth_matrix,methcol,Y) {
  mod = glm(Y~meth_matrix[, methcol],family = "binomial")
  cf = summary(mod)$coefficients
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]  
}

system.time(ind.res <- mclapply(setNames(seq_len(ncol(M_matrix)), dimnames(M_matrix)[[2]]), GLMtest, meth_matrix=M_matrix, Y=pheno.ses[,2]))
```

```
##    user  system elapsed 
##   0.001   0.003   0.147
```

```r
all.results<-ldply(ind.res,rbind)
names(all.results)<-c("probeID","BETA","SE", "P_VAL")
```
* show the significant results

```r
all.results[which(all.results$P_VAL<0.05),]
```

```
##       probeID       BETA        SE       P_VAL
## 33 cg23005227 -1.0666431 0.4470010 0.017022537
## 36 cg04891688 -1.1049380 0.4277655 0.009793222
## 65 cg18749617 -0.4186703 0.2103032 0.046503588
```
