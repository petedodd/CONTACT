## top level handler for CONTACT

## flags for sensitivity analyses
shell <- FALSE # whether running from shell script or not
if(shell){
  ## running from shell
  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  SA <- args[1]                  # none,base/lo/hi,tptru,hicoprev
  if(SA == 'none'){
    SA <- ''
  } 
} else { #set by hand
  rm(list=ls()) #clear all
  shell <- FALSE #whether running from shell script or not
  ##sensitivity analyses (mostly for PT):
  ## '' = basecase
  ## 'discr'='base'/'lo'/'hi'
  ## 'cdr' = making cdr higher for incidence
  ## 'txd' = making the completion influence tx/pt outcome
  sacases <- c('','lo','tptru','hicoprev', 'ctryeff','ugaattcsts', 'cdr')
  SA <- sacases[1]
}

# rm(list=ls())
library(here)
library(tidyverse)

## load other scripts
source(here('R/contact_tree.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/contact_functions.R'))      #functions for tree parameters

## number of reps
nreps <- 1e3
set.seed(1234)

## attributes to use
# tblevels <- c('TB+','TB-','noTB') #bac confirmable TB, bac unconfirmable TB, not TB
hivlevels <- c(0,1)
artlevels <- c(0,1)
agelevels <- c('0-4','5-14')
isoz <- c('CMR','UGA') #relevant countries


## --- life years and other outputs NOTE needs to be set FALSE on first run thru
LYSdone <- TRUE
if(!LYSdone){
  ## make discounted life-years if they haven't been done
  LYKc <- GetLifeYears(isolist=isoz,discount.rate=0.03,yearfrom=2021)
  LYKc0 <- GetLifeYears(isolist=isoz,discount.rate=0.00,yearfrom=2021)
  LYKc5 <- GetLifeYears(isolist=isoz,discount.rate=0.05,yearfrom=2021)
  LYKc <- merge(LYKc,LYKc0[,.(iso3,age,LYS0=LYS)],by=c('iso3','age'))
  LYKc <- merge(LYKc,LYKc5[,.(iso3,age,LYS5=LYS)],by=c('iso3','age'))
  LYK <- LYKc[,.(LYS=mean(LYS),LYS0=mean(LYS0),LYS5=mean(LYS5)),by=.(age)] #averaged life-years 4 generic tests
  save(LYKc,file=here('indata/LYKc.Rdata'))
  save(LYK,file=here('indata/LYK.Rdata'))
} else {
  load(file=here('indata/LYKc.Rdata'))
  load(file=here('indata/LYK.Rdata'))
}

if(SA %in% c('hi','lo')){
  LYKc[,LYS:=ifelse(SA=='lo', LYS0, 
                    ifelse(SA=='hi',LYS5, LYS))]  
}

## prior parameters
PD0 <- read.csv(here('indata/parameter_distributions.csv')) #read in
## parameters to be determined from cascade data
PD1 <- PD0[PD0$DISTRIBUTION=="",]
## the rest
PD0 <- PD0[PD0$DISTRIBUTION!="",]

drop <- PD0$NAME[grepl('sens|spec', PD0$NAME)]
PD0 <- PD0 %>% dplyr::filter(!NAME %in% drop)

# # proportions from cascade data (age aggregated)
PPA <- fread(here('indata/proportions.csv'))
# 
# # quick renaming some variables. TODO: rename in analysis code later
PPA[variable=='frac.clin.7d.dx', variable:='frac.clin.7d.clin.dx']
PPA[variable=='frac.noclin.7d.dx', variable:='frac.clin.7d.noclin.dx']

PPA <- melt(PPA, id.vars = c('country', 'variable'), variable.name = 'care_model')
PPA <- PPA[care_model=='INT', variable:=paste('int.', variable, sep = '')]
PPA <- PPA[care_model=='SOC', variable:=paste('soc.', variable, sep = '')]

PPA[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]

keep <- vrz[grepl('soc|int', vrz)]
PPA <- PPA %>% dplyr::filter(variable %in% keep)

ParmsTab2 <- PPA    # save for later use in creating parameter table for Appendix
fwrite(ParmsTab2,file=here('outdata/Parameters2.csv'))

PPA <- dcast(PPA, isoz~variable)

# load cohort data (fraction under 5 and HIV)
# load(here('indata/cohortdata.Rdata'))
# 
# popn <- setDT(popn)
# popn <- popn[country=='Cameroun', country:='Cameroon']
# popn <- popn[age_cat=='u5', variable:='F.u5']
# popn <- popn[age_cat=='o5', variable:='F.o5']
# popn <- popn[,.(country=as.character(country), age_cat=as.character(age_cat), variable, rando, prop)]
# popn <- popn[rando=='ITV', rando:='INT']
# popn <- dcast(popn, country+age_cat+variable~rando)
# 
# PP <- rbind(PP, popn)

# intervention effect estimates
INTE <- fread(gh('indata/pooled_effects.csv'))    #read cost data

# INTE[, c('mid', 'lo', 'hi') := tstrsplit('MEDIAN (IQR)', "", fixed=TRUE)][]

tmp <- INTE %>%
  extract(`MEDIAN (IQR)`, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") %>%
  separate(lo,c("lo","hi"),"-") %>%
  mutate_at(c("mid", "lo","hi"), as.numeric)

tmp <- setDT(tmp)
tmp1 <- getLNparms(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2,med=FALSE)
tmp[,DISTRIBUTION:=paste0("LN(",tmp1$mu,",",tmp1$sig,")")] #LN distributions

INTE[,DISTRIBUTION:=tmp$DISTRIBUTION]
INTE[DISTRIBUTION=='LN(NA,NA)', DISTRIBUTION:=NA]


ParmsTab3 <- copy(INTE)     # save for later use in creating parameter tables for Appendix
ParmsTab3[,DISTRIBUTION:=paste0("LN(",round(tmp1$mu,5),",",round(tmp1$sig,5),")")] #LN distributions
fwrite(ParmsTab3,file=here('outdata/Parameters3.csv'))

# # proportions from cascade data
PP <- fread(here('indata/proportions_age_cat.csv'))

# quick renaming some variables. TODO: rename in analysis code later
PP[variable=='frac.clin.7d.dx', variable:='frac.clin.7d.clin.dx']
PP[variable=='frac.noclin.7d.dx', variable:='frac.clin.7d.noclin.dx']

PP <- melt(PP, id.vars = c('country', 'age_cat', 'variable'), variable.name = 'care_model')
PP <- PP[care_model=='INT', variable:=paste('int.', variable, sep = '')]
PP <- PP[care_model=='SOC', variable:=paste('soc.', variable, sep = '')]

PP[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
PP[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]

keep <- vrz[grepl('soc|int', vrz)]
PP <- PP %>% dplyr::filter(variable %in% keep)

## PJD
PP[,.N,by=.(isoz,age,care_model)] #missing some UGA older ages
cmrlist <- PP[isoz=='CMR' & age=='5-14',unique(variable)]
(missed <- setdiff(cmrlist,PP[isoz=='UGA' & age=='5-14',unique(variable)]))
## [1] "int.frac.clin.7d.clin.dx"   "int.frac.clin.7d.noclin.dx"
## [3] "soc.frac.clin.7d.clin.dx"   "soc.frac.clin.7d.noclin.dx"
tmp <- PP[variable %in% missed & isoz=='UGA'] #u5 only
tmp <- copy(tmp)
tmp[,c('age_cat','age'):=.('o5','5-14')]
PP <- rbind(PP,tmp) #TODO BUG correction by making some missing rows

ParmsTab2a <- PP    # save for later use in creating parameter table for Appendix
fwrite(ParmsTab2a,file=here('outdata/Parameters2a.csv'))

PP <- dcast(PP, isoz+age~variable)
# 

# load enrollment data (number of contacts per index case & fraction under 5 years)
# also available but not currently being used: declared per household/enrolled, enrolled per household, frac declared u5, frac declared/enrolled hiv
library(readxl)
base_popn <- data.table(read_excel(here("indata/Baseline information.xlsx"), sheet = 'Sheet1', range = 'Z13:AB17'))
base_popn <- melt(base_popn, variable.name = 'isoz')
age_split <- base_popn[grepl('enrolled.', metric), .(isoz, variable=metric, value)]
age_split <- age_split[isoz=='CMR', variable:=paste('cmr.', variable, sep = '')]
age_split <- age_split[isoz=='UGA', variable:=paste('uga.', variable, sep = '')]

ParmsTab4 <- base_popn    # save for later use in creating parameter table for Appendix
fwrite(ParmsTab4,file=here('outdata/Parameters4.csv'))

tmp <- data.table(matrix(NA, nrow = nrow(age_split), ncol = ncol(PD0)))
names(tmp) <- names(PD0)
tmp <- tmp[,NAME:=age_split$variable]
tmp <- tmp[,DISTRIBUTION:=age_split$value]
age_split <- copy(tmp)
age_split <- age_split[,NAME:=gsub('enrolled.', '',NAME)]

# commented CODE below not working

DENR <- data.table(read_excel(here("indata/Baseline information.xlsx"), sheet = 'Sheet1', range = 'O3:R19'))
DENR <- melt(DENR, variable.name = 'isoz')

# declared/enrolled plot
DENRP <- DENR[grepl('per_index_case', metric), .(isoz, model,variable=metric, value)]
DENRP <- DENRP[,variable:=gsub('_per_index_case', '',variable)]


DENR <- DENR[,model:=toupper(model)]
DENR <- DENR[model=='INT', metric:=paste('int.', metric, sep = '')]
DENR <- DENR[model=='SOC', metric:=paste('soc.', metric, sep = '')]
DENR <- DENR[grepl('enrolled.per.index', metric), .(isoz, variable=metric, care_model=model, value)]
DENR <- DENR[,variable:=gsub('olled', '',variable)]

DENR <- dcast(DENR, isoz~variable)

TBPREV <- fread(here("indata/tb_epi.csv"))
TBPREV[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]
TBPREV <- TBPREV[,.(isoz, tb_prev)]
names(TBPREV) <- c('isoz', 'int.tbprev.symptomatic')

if(SA=='hicoprev'){       
  TBPREV[,int.tbprev.symptomatic:=0.1]              # higher co-prevalence of tuberculosis disease among household child contacts (10%?)
                                                    # 10% (5.0–18.9%) in children < 5 years & 8.4% (2.8–22.6%) for children 5-14 years
}

# pull together 
PP <- merge(PP, DENR, by='isoz') # age disaggregated
PP <- merge(PP, TBPREV, by='isoz')
PPA <- merge(PPA, DENR, by='isoz') # aggregated
PPA <- merge(PPA, TBPREV, by='isoz')
PPA <- rbind(PPA,PPA)
PPA[,age:=rep(agelevels,each=nrow(PPA)/2)]

# adding study data on: frac u5, HIV prevalence and ART coverage
# currently using CMR data

names(INTE) <- names(PD0)
INTE <- INTE[INTE$DISTRIBUTION!="",]
PD0 <- rbind(PD0, age_split, INTE)

PD0 <- PD0[!is.na(PD0$DISTRIBUTION),]
P1 <- parse.parmtable(PD0)             #convert into parameter object
# P2 <- parse.parmtable(P2)             #convert into parameter object
# PZ <- c(P1,P2)
names(P1)
P <- P1
## make base PSA dataset
set.seed(1234) #random number seed
# D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))
D <- makePSA(nreps,P,dbls = list(c('hivartOR:mn','hivartOR:sg')))

## ## NOTE temporary introduction of noise:
## for(nm in setdiff(PD1$NAME,'d.OR.dh.if.TB')){ #loop over probs
##   D[[nm]] <- ilogit(logit(D[[nm]]) + rnorm(nreps)/5)
## }
## D[['d.OR.dh.if.TB']] <- exp(log(D[['d.OR.dh.if.TB']]) + rnorm(nreps)/5)

## use these parameters to construct input data by attribute
D <- makeAttributes(D)
D[,sum(value),by=.(isoz,id)] #CHECK
# D[tb!='noTB',sum(value),by=id] #CHECK
D[,sum(value),by=.(isoz,id,age)] #CHECK

load(file=here('outdata/CDR.Rdata')) #CDR

## make PSA for country CDRs
CDRs <- CDR[qty=='cdr']
CDRs <- CDRs[rep(1:nrow(CDRs),each=max(D$id))]
CDRs[,id:=rep(1:max(D$id),nrow(CDRs)/max(D$id))]
CDRs[,sz:=value*(1-value)/cdr.v-1]
CDRs[,c('a','b'):=.(sz*value,sz*(1-value))]
CDRs[,cdr0:=rbeta(nrow(CDRs),a,b)]
CDRs[,cdr:=cdr0] #basecase
# CDRs[,mn:=value*(1+runif(nrow(CDRs)))]  
# CDRs[,mn:=pmin(mn,1)]
# CDRs[,cdri:=rbeta(nrow(CDRs),mn*cdr.v,(1-mn)*cdr.v)] # TB detection if household visited
CDRs[,cdri:=aCDR(value,cdr.v)] # TB detection if household visited
if(SA=='cdr'){   #sensitivity analysis
  CDRs[,cdr:=runif(nrow(CDRs))]
  CDRs[,cdr:=cdr0*(1-cdr) + cdr] #interpolate between cdr0 & 1
}
CDRs[,summary(cdr)]
CDRs[iso3=='CMR',summary(cdr)]
CDRs[iso3=='UGA',summary(cdr)]

CDRs[,summary(cdr)]*1.5
CDRs[iso3=='CMR',summary(cdr)]*1.5
CDRs[iso3=='UGA',summary(cdr)]*1.5

## merge in CDRs
D <- merge(D,CDRs[,.(isoz=iso3,age,id,cdr,cdri)],by=c('isoz','age','id'))
# D <- rbind(D,D)
# D[,isoz:=rep(c('CMR','UGA'),each=nrow(D)/2)]
# D <- merge(D,PP,by=c('isoz', 'age'),all.x = TRUE)
D <- merge(D,PPA,by=c('isoz', 'age'),all.x = TRUE)


# TODO: figure out how to add country and model of care specific data

if(SA=='ctryeff'){
  D[,tpt.resultOR:=ifelse(isoz=='CMR', cmr.tpt.resultOR, uga.tpt.resultOR)]
  D[,tpt.initiationOR:=ifelse(isoz=='CMR', cmr.tpt.initiationOR, uga.tpt.initiationOR)]
  D[,tpt.completionOR:=ifelse(isoz=='CMR', cmr.tpt.completionOR, uga.tpt.completionOR)]
  D[,tb.resultOR:=ifelse(isoz=='CMR', cmr.tb.resultOR, uga.tb.resultOR)]
  D[,tb.diagnosisOR:=ifelse(isoz=='CMR', cmr.tb.diagnosisOR, uga.tb.diagnosisOR)]
}



## read and make cost data
# csts <- fread(here('indata/testcosts.csv'))         #read cost data
# rcsts <- fread(gh('indata/model_mean_total_costs_age.csv'))    #read cost data
# names(rcsts)
# 
# rcsts <- rcsts[rando=='INT', cascade:=paste('c.int.', cascade, sep = '')]
# rcsts <- rcsts[rando=='SOC', cascade:=paste('c.soc.', cascade, sep = '')]
# 
# rcsts[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
# rcsts[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]
# 
# keeps <- c('isoz', 'age', 'rando', 'cascade','cost.m', 'cost.sd')
# rcsts <- rcsts[,..keeps]
# 
# rcsts <- add_row(rcsts, isoz='CMR', age='5-14', cascade='c.soc.tpt', cost.m=0, cost.sd=0)
# 
# ## turn cost data into PSA
# rcsts[cost.sd==0,cost.sd:=cost.m/40]        #SD such that 95% UI ~ 10% of mean
# allcosts <- rcsts[,.(iso3=isoz, age=age, cost=cascade, cost.m, cost.sd)]
# allcosts[is.na(allcosts)] <- 0 #some quick fix >> setting NA to 0
# C <- MakeCostData(allcosts[iso3=='CMR'],nreps)               #make cost PSA NOTE using CMR cost data

# rcsts <- fread(gh('indata/model_mean_total_costs.csv'))    #read cost data
rcsts <- fread(gh('indata/model_mean_total_costs.csv'))    #read revised cost data

if(SA=='tptru'){
  rcsts <- fread(gh('indata/model_mean_total_costs_red.csv'))    #read cost data
}

names(rcsts)

if(SA=='ugaattcsts'){
rcsts <- rcsts[!(country=='Uganda' & cascade %in% cascade[grepl('att', cascade)]),]
tmp <- rcsts[country=='Cameroon' & cascade %in% cascade[grepl('att', cascade)],]
tmp <- tmp[, country:='Uganda']
rcsts <- rbind(rcsts, tmp)
}
# rcsts[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
rcsts[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]

rcsts <- rcsts[rando=='INT', cascade:=paste('c.int.', cascade, sep = '')]
rcsts <- rcsts[rando=='SOC', cascade:=paste('c.soc.', cascade, sep = '')]

keeps <- c('isoz', 'rando', 'cascade','cost.m', 'cost.sd')
rcsts <- rcsts[,..keeps]

## turn cost data into PSA
rcsts[is.na(rcsts)] <- 0 #some quick fix >> setting NA to 0
rcsts[cost.sd==0,cost.sd:=cost.m/40]        #SD such that 95% UI ~ 10% of mean
# rcsts[cost.sd>100,cost.sd:=cost.m/40]
# rcsts[cost.m>100,cost.m:=cost.m/20]
# rcsts[,cost.m:=cost.m/10]
# rcsts[cascade %in% cascade[grepl('tpt', cascade)],cost.m:=cost.m/10]
allcosts <- rcsts[,.(iso3=isoz, cost=cascade, cost.m, cost.sd)]

keep <- vrz[grepl('c.soc.|c.int.', vrz)]

# check
setdiff(unique(allcosts$cost[allcosts$iso3=='UGA']),
        unique(keep))

setdiff(unique(keep),
        unique(allcosts$cost[allcosts$iso3=='UGA']))

# allcosts[24,cost:='c.soc.inv_inc']
C <- MakeCostData(allcosts[iso3=='CMR'],nreps)               #make cost PSA NOTE using CMR cost data

## compute other parameters (adds by side-effect)
MakeTreeParms(D,P)
## soc and int frac.screened

## checks
D[,sum(value),by=.(isoz,id)] #CHECK
D[,sum(value),by=.(isoz,id,age)] #CHECK
D[,.(isoz,age,F.u5,hivprev.u5)]
D[,.(isoz,age,soc.frac.screened,tb.screeningOR,int.frac.screened)]
## check for leaks
head(SOC.F$checkfun(D)) #SOC arm
head(INT.F$checkfun(D)) #INT arm


## === REACH LOGIC
## Logic for coverage etc:
## index: 558, 341 (int/soc)
## declared: 1889, 1005
## enrolled: 1835, 498
## declared/index = c(1889/558, 1005/241) #3.4, 4.2
## enrolled/index = c(1835/558, 498/241) #3.4, 2.0

## 1) multiply value by 3.4 (or logical equivalent)
# D[,value:=value * 3.4]
D[,value:=value * int.enr_per_index_case]  # applying coverage by country
DENR  # child contacts enrolled per index by country. Huge difference in SOC
DENR[,.(soc.screened=soc.enr_per_index_case/int.enr_per_index_case), by=.(isoz)]     # pooled data close to UGA specific hence no change to UGA result  
## 2) set int.frac.screened = 1, and soc.frac.screened = 2.0/3.4
# D[,c('soc.frac.screened','int.frac.screened'):=.(2.0/3.4,1.0)]
D[,c('soc.frac.screened','int.frac.screened'):=.(soc.enr_per_index_case/int.enr_per_index_case,1.0)]


## TODO approximate for now

names(SOC.F)

## add cost data
D <- merge(D,C,by=c('id'),all.x=TRUE)       #merge into PSA

## === RUN MODEL
arms <- c('SOC','INT')
D <- runallfuns(D,arm=arms)                      #appends anwers

summary(D)



## --- run over different countries
cnmz <- names(C) # C=cost PSA data
cnmz <- cnmz[cnmz!=c('id')]
toget <- c('id',
           'cost.soc','cost.int',
           'cost.screen.soc',	'cost.screen.int',
           'cost.tpt.soc', 'cost.tpt.int',
           'cost.prev.att.soc','cost.prev.att.int',
           'cost.inc.att.soc','cost.inc.att.int',
           'tpt.soc','tpt.int',
           'att.soc','att.int',
           'deaths.soc','deaths.int',
           'LYS','LYS0','value'
           )
toget2 <- c(toget,
            'prevtb.soc','prevtb.int',
            'inctb.soc','inctb.int',
            'incdeaths.soc','incdeaths.int')
notwt <- c('id','LYS','LYS0','value') #variables not to weight against value
lyarm <- c('LYL.soc','LYL.int')
lyarm <- c(lyarm,gsub('\\.','0\\.',lyarm)) #include undiscounted
tosum <- c(setdiff(toget,notwt),lyarm)
tosum2 <- c(tosum,
            'prevtb.soc','prevtb.int',
            'inctb.soc','inctb.int',
            'incdeaths.soc','incdeaths.int')
## heuristic to scale top value for thresholds:
heur <- c('id','value','deaths.int','deaths.soc')
out <- D[,..heur]
out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=c('deaths.int','deaths.soc'),by=id] #sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.int)]
topl <- 5e3
lz <- seq(from = 0,to=topl,length.out = 1000) #threshold vector for CEACs

## containers & loop
allout <- allpout <- list() #tabular outputs
allout2 <- allpout2 <- list() #tabular outputs
ceacl <- NMB <- list()             #CEAC outputs etc
psaout <- psapout <- list()
parmsout <- parmsout <- list()
## NOTE I think there was an additional problem here -
## because countries are contained as separate rows we were summing over both in the below computations
## ie we were looping but then operating over data for both countries
## cn <- isoz[1]

for(cn in isoz){
  dc <- D[isoz==cn]
  cat('running model for:',cn, '0-14 years\n')
  ## --- costs
  ## drop previous costs
  # D <- D[age=='0-4',]
  dc[,c(cnmz):=NULL]
  ## add cost data
  C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
  dc <- merge(dc,C,c('id'),all.x=TRUE)        #merge into PSA
  ## --- DALYs
  ## drop any that are there
  if('LYS' %in% names(D)) dc[,c('LYS','LYS0'):=NULL]
  dc <- merge(dc,LYKc[iso3==cn,.(age,LYS,LYS0)],by='age',all.x = TRUE)        #merge into PSA
  ## --- run model (quietly)
  invisible(capture.output(dc <- runallfuns(dc,arm=arms)))
  ## --- grather outcomes
  out <- dc[,..toget]
  out[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
                   LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
  out <- out[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum,by=id] #sum against popn
  ## non-incremental cost per ATT
  out[,costperATT.soc:=cost.soc/att.soc];
  out[,costperATT.int:=cost.int/att.int]; 
  out[,costperTPT.soc:=cost.soc/tpt.soc];
  out[,costperTPT.int:=cost.int/tpt.int]; 
  ## increments wrt SOC (per child screened at either facility/household)
  out[,Dcost.int:=cost.int-cost.soc];  #inc costs
  out[,Datt.int:=att.int-att.soc];  #inc atts
  out[,Dtpt.int:=tpt.int-tpt.soc];  #inc atts
  out[,Ddeaths.int:=deaths.int-deaths.soc];  #inc deaths
  out[,DLYL0.int:=LYL0.int-LYL0.soc]; #inc LYLs w/o discount
  out[,DLYL.int:=LYL.int-LYL.soc]; #inc LYLs
  ## per whatever
  out[,DcostperATT.int:=Dcost.int/Datt.int];
  out[,DcostperTPT.int:=Dcost.int/Dtpt.int];
  out[,Dcostperdeaths.int:=-Dcost.int/Ddeaths.int];
  out[,DcostperLYS0.int:=-Dcost.int/DLYL0.int];
  out[,DcostperLYS.int:=-Dcost.int/DLYL.int];
  ## summarize
  smy <- outsummary(out)
  outs <- smy$outs; pouts <- smy$pouts;
  outs[,iso3:=cn]; pouts[,iso3:=cn]
  ## capture tabular
  allout[[cn]] <- outs; allpout[[cn]] <- pouts
  ## capture data for NMB
  NMB[[cn]] <- out[,.(iso3=cn,DLYL.int,Dcost.int)]
  ## ceac data
  ceacl[[cn]] <- data.table(iso3=cn,
                            int=make.ceac(out[,.(Q=-DLYL.int,P=Dcost.int)],lz),
                            threshold=lz)
  ## --- grather outcomes for Table 2
  out2 <- dc[,..toget2]
  out2[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
                    LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
  cntcts <- out2[,sum(value),by=id]
  out2 <- out2[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum2,by=id] #sum against popn
  out2[,contacts.int:=cntcts$V1]
  # out2[,contacts.soc:=2.0] #TODO needs changing
  out2[,contacts.soc:=ifelse(cn=='CMR', 1.202703, 1.941176)] # by country TODO needs checking
  out2[,c('LYL0.soc','LYL0.int'):=NULL] #drop
  out2[,c('prevdeaths.soc','prevdeaths.int'):=.(deaths.soc-incdeaths.soc,deaths.int-incdeaths.int)]
  ## increments
  out2[,Dtpt:=tpt.int-tpt.soc]
  out2[,Datt:=att.int-att.soc]
  out2[,DLYL:=LYL.int-LYL.soc]
  out2[,Dprevtb:=prevtb.int-prevtb.soc]
  out2[,Dinctb:=inctb.int-inctb.soc]
  out2[,Dincdeaths:=incdeaths.int-incdeaths.soc]
  out2[,Dprevdeaths:=prevdeaths.int-prevdeaths.soc]
  out2[,Ddeaths:=deaths.int-deaths.soc]
  out2[,Dcontacts:=contacts.int-contacts.soc]
  # out2[,cost.att.soc:=cost.soc/att.soc];
  # out2[,cost.att.int:=cost.int/att.int]; 
  # out2[,cost.tpt.soc:=cost.soc/tpt.soc];
  # out2[,cost.tpt.int:=cost.int/tpt.int];
  out2[,Dcost.screen:=cost.screen.int-cost.screen.soc]
  out2[,Dcost.tpt:=cost.tpt.int-cost.tpt.soc]
  out2[,Dcost.prev.att:=cost.prev.att.int-cost.prev.att.soc]
  out2[,Dcost.inc.att:=cost.inc.att.int-cost.inc.att.soc]
  out2[,Dcost:=cost.int-cost.soc]

  ## summarize
  smy2 <- Table2(out2) #NOTE set per 1000 index cases or HHs - adjust fac in contact_functions.R
  outs2 <- smy2$outs; pouts2 <- smy2$pouts;
  outs2[,iso3:=cn]; pouts2[,iso3:=cn]
  ## capture tabular
  allout2[[cn]] <- outs2; allpout2[[cn]] <- pouts2
  ## capture data for NMB
  # NMB2[[cn]] <- out2[,.(iso3=cn,DLYL,Dcost)]
  # psaout1[[cn]] <- out2[,.(iso3=cn,cost.soc,cost.int,LYL.soc,LYL.int,DLYL,Dcost)]
  # ## ceac data
  # ceacl2[[cn]] <- data.table(iso3=cn,
  #                           int=make.ceac(out2[,.(Q=-DLYL,P=Dcost)],lz),
  #                           threshold=lz)
  psaout[[cn]] <- dc[,.(iso3=cn,
                        cost.soc=sum(cost.soc*value),
                        cost.int=sum(cost.int*value),
                        lyl.soc=sum(deaths.soc*value*LYS),
                        lyl.int=sum(deaths.int*value*LYS)),
                     by=.(id, isoz, age)] #PSA summary
  
  
  parms_epicalc <- c('cfr.notx','cfr.tx', 'p.tbdx.1yr', 'tptRR', 'CDR', 'CDRi')
  parms_inteff <- c('tpt.resultOR','tpt.initiationOR','tpt.completionOR','tb.resultOR','tb.diagnosisOR')
  parms_cstssoc <- vrz[grepl('c.soc', vrz)]
  parms_cstsint <- vrz[grepl('c.int', vrz)]
  
  vars <- c(parms_epicalc, parms_inteff, parms_cstssoc, parms_cstsint) 
  parms <- dc[,lapply(.SD,function(x) mean(x, na.rm = T)),.SDcols=vars,by=.(id, isoz, age)]
  # parms <- cbind(dc[,.(id, iso3=isoz, age)],
  #                dc[,..vars])
  parms[,iso3:=cn]
  parms[,isoz:=NULL]
  parmsout[[cn]] <- parms
}

# summary(D[,..toget])
# D[,.(isoz,age,F.u5,value)] # why are some values 0? TODO BUG ?
# D[,summary(tpt.int/tpt.soc)]
# 
# check1 <- names(labz)[4:14]
# check1 <- c(paste0(check1,'.soc'), paste0(check1,'.int'))
# check2 <- names(labz)[14:21]
# summary(D[,..check1])
# odds <- names(D)[grepl('OR',names(D))]
# summary(D[,..odds])

allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
allout2 <- rbindlist(allout2)
allpout2 <- rbindlist(allpout2)
psaout <- rbindlist(psaout)
parmsout <- rbindlist(parmsout)

ceacl <- rbindlist(ceacl)
NMB <- rbindlist(NMB)

fwrite(allout,file=gh('outdata/allout') + SA + '.csv')
fwrite(allpout,file=gh('outdata/allpout') + SA + '.csv')
fwrite(allout2,file=gh('outdata/allout2') + SA + '.csv')
fwrite(allpout2,file=gh('outdata/allpout2') + SA + '.csv')
save(ceacl,file=gh('outdata/ceacl') + SA + '.Rdata')
save(NMB,file=gh('outdata/NMB') + SA + '.Rdata')

ICERS <- allpout[,.(iso3,ICER=ICER.int)]
fwrite(ICERS,file=gh('outdata/ICERS') + SA + '.csv')
ICERS
## CEAC plot
cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73",
               "#F0E442", "#0072B2","#D55E00", "#CC79A7")
ceaclm <- melt(ceacl,id=c('iso3','threshold'))
ceaclm[,Intervention:=ifelse(variable=='int','Intervention','SOC')]
## name key
ckey <- data.table(iso3=c('CMR','UGA'),
                   country=c('Cameroon','Uganda'))

ceaclm <- merge(ceaclm,ckey,by='iso3',all.x=TRUE)

## plot
GP <- ggplot(ceaclm,aes(threshold,value,
                        col=country,lty=Intervention)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette)
GP

ggsave(GP,file=gh('plots/CEAC') + SA + '.png',w=7,h=5)



## plot
GP <- ggplot(ceaclm[variable=='int'],aes(threshold,value,
                        col=country)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette) ## + xlim(x=c(0,1500))
GP

ggsave(GP,file=gh('plots/CEAC1') + SA + '.png',w=7,h=5)



## ## generate some CEA outputs in graphs/ & outdata/
## ## NOTE these folders need to be created
## ## NOTE need ggpubr, BCEA installed
# MakeCEAoutputs(D, #PSA dataset
#                LYK, #discounted expected life-years by age
#                file.id='test', #string to identify output files
#                Kmax=5e3,wtp=5e3)


## --- run over different countries & age groupd

# ## containers & loop
# allout <- allpout <- list() #tabular outputs
# allout2 <- allpout2 <- list() #tabular outputs
# ceacl <- NMB <- list()             #CEAC outputs etc
# psaout <- psapout <- list()
# ## NOTE I think there was an additional problem here -
# ## because countries are contained as separate rows we were summing over both in the below computations
# ## ie we were looping but then operating over data for both countries
# ## cn <- isoz[1]
#
# for(cn in isoz){
  # dc <- D[isoz==cn & age=='0-4']
#   cat('running model for:',cn,'0-4 years\n')
#   ## --- costs
#   ## drop previous costs
#   # D <- D[age=='0-4',]
#   dc[,c(cnmz):=NULL]
#   ## add cost data
#   C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
#   dc <- merge(dc,C,c('id'),all.x=TRUE)        #merge into PSA
#   ## --- DALYs
#   ## drop any that are there
#   if('LYS' %in% names(D)) dc[,c('LYS','LYS0'):=NULL]
#   dc <- merge(dc,LYKc[iso3==cn,.(age,LYS,LYS0)],by='age',all.x = TRUE)        #merge into PSA
#   ## --- run model (quietly)
#   invisible(capture.output(dc <- runallfuns(dc,arm=arms)))
#   ## --- grather outcomes
#   out <- dc[,..toget]
#   out[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
#                    LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
#   ## out[,sum(value),by=id]                                       #CHECK
#   ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
#   out <- out[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum,by=id] #sum against popn
#   ## non-incremental cost per ATT
#   out[,costperATT.soc:=cost.soc/att.soc];
#   out[,costperATT.int:=cost.int/att.int];
#   out[,costperTPT.soc:=cost.soc/tpt.soc];
#   out[,costperTPT.int:=cost.int/tpt.int];
#   ## increments wrt SOC (per child screened at either facility/household)
#   out[,Dcost.int:=cost.int-cost.soc];  #inc costs
#   out[,Datt.int:=att.int-att.soc];  #inc atts
#   out[,Dtpt.int:=tpt.int-tpt.soc];  #inc atts
#   out[,Ddeaths.int:=deaths.int-deaths.soc];  #inc deaths
#   out[,DLYL0.int:=LYL0.int-LYL0.soc]; #inc LYLs w/o discount
#   out[,DLYL.int:=LYL.int-LYL.soc]; #inc LYLs
#   ## per whatever
#   out[,DcostperATT.int:=Dcost.int/Datt.int];
#   out[,DcostperTPT.int:=Dcost.int/Dtpt.int];
#   out[,Dcostperdeaths.int:=-Dcost.int/Ddeaths.int];
#   out[,DcostperLYS0.int:=-Dcost.int/DLYL0.int];
#   out[,DcostperLYS.int:=-Dcost.int/DLYL.int];
#   ## summarize
#   smy <- outsummary(out)
#   outs <- smy$outs; pouts <- smy$pouts;
#   outs[,iso3:=cn]; pouts[,iso3:=cn]
#   ## capture tabular
#   allout[[cn]] <- outs; allpout[[cn]] <- pouts
#   ## capture data for NMB
#   NMB[[cn]] <- out[,.(iso3=cn,DLYL.int,Dcost.int)]
#   ## ceac data
#   ceacl[[cn]] <- data.table(iso3=cn,
#                             int=make.ceac(out[,.(Q=-DLYL.int,P=Dcost.int)],lz),
#                             threshold=lz)
#   ## --- grather outcomes for Table 2
#   out2 <- dc[,..toget2]
#   out2[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
#                     LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
#   ## out[,sum(value),by=id]                                       #CHECK
#   ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
#   cntcts <- out2[,sum(value),by=id]
#   out2 <- out2[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum2,by=id] #sum against popn
#   out2[,contacts.int:=cntcts$V1]
#   # out2[,contacts.soc:=2.0] #TODO needs changing
#   out2[,contacts.soc:=ifelse(cn=='CMR', 1.202703, 1.941176)] # by country TODO needs checking
#   out2[,c('LYL0.soc','LYL0.int'):=NULL] #drop
#   out2[,c('prevdeaths.soc','prevdeaths.int'):=.(deaths.soc-incdeaths.soc,deaths.int-incdeaths.int)]
#   ## increments
#   out2[,Dtpt:=tpt.int-tpt.soc]
#   out2[,Datt:=att.int-att.soc]
#   out2[,DLYL:=LYL.int-LYL.soc]
#   out2[,Dprevtb:=prevtb.int-prevtb.soc]
#   out2[,Dinctb:=inctb.int-inctb.soc]
#   out2[,Dincdeaths:=incdeaths.int-incdeaths.soc]
#   out2[,Dprevdeaths:=prevdeaths.int-prevdeaths.soc]
#   out2[,Ddeaths:=deaths.int-deaths.soc]
#   out2[,Dcontacts:=contacts.int-contacts.soc]
#   # out2[,cost.att.soc:=cost.soc/att.soc];
#   # out2[,cost.att.int:=cost.int/att.int];
#   # out2[,cost.tpt.soc:=cost.soc/tpt.soc];
#   # out2[,cost.tpt.int:=cost.int/tpt.int];
#   out2[,Dcost.screen:=cost.screen.int-cost.screen.soc]
#   out2[,Dcost.tpt:=cost.tpt.int-cost.tpt.soc]
#   out2[,Dcost.prev.att:=cost.prev.att.int-cost.prev.att.soc]
#   out2[,Dcost.inc.att:=cost.inc.att.int-cost.inc.att.soc]
#   out2[,Dcost:=cost.int-cost.soc]
#
#   ## summarize
#   smy2 <- Table2(out2) #NOTE set per 1000 index cases or HHs - adjust fac in contact_functions.R
#   outs2 <- smy2$outs; pouts2 <- smy2$pouts;
#   outs2[,iso3:=cn]; pouts2[,iso3:=cn]
#   ## capture tabular
#   allout2[[cn]] <- outs2; allpout2[[cn]] <- pouts2
#
#   psaout[[cn]] <- dc[,.(iso3=cn,
#                         cost.SOC=sum(cost.soc*value),
#                         cost.INT=sum(cost.int*value),
#                         lyl.SOC=sum(deaths.soc*value*LYS),
#                         lyl.INT=sum(deaths.int*value*LYS)),
#                      by=id] #PSA summary
#
# }
#
#
# # summary(D[,..toget])
# # D[,.(isoz,age,F.u5,value)] # why are some values 0? TODO BUG ?
# # D[,summary(tpt.int/tpt.soc)]
# #
# # check1 <- names(labz)[4:14]
# # check1 <- c(paste0(check1,'.soc'), paste0(check1,'.int'))
# # check2 <- names(labz)[14:21]
# # summary(D[,..check1])
# # odds <- names(D)[grepl('OR',names(D))]
# # summary(D[,..odds])
#
# allout <- rbindlist(allout)
# allpout <- rbindlist(allpout)
# allout2 <- rbindlist(allout2)
# allpout2 <- rbindlist(allpout2)
# psaout <- rbindlist(psaout)
#
# ceacl <- rbindlist(ceacl)
# NMB <- rbindlist(NMB)
#
# fwrite(allout,file=gh('outdata/alloutY') + SA + '.csv')
# fwrite(allpout,file=gh('outdata/allpoutY') + SA + '.csv')
# fwrite(allout2,file=gh('outdata/allout2Y') + SA + '.csv')
# fwrite(allpout2,file=gh('outdata/allpout2Y') + SA + '.csv')
# save(ceacl,file=gh('outdata/ceaclY') + SA + '.Rdata')
# save(NMB,file=gh('outdata/NMBY') + SA + '.Rdata')
#
# ICERSY <- allpout[,.(iso3,ICER=ICER.int)]
# fwrite(ICERSY,file=gh('outdata/ICERSY') + SA + '.csv')
#
# ## CEAC plot
# cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73",
#                         "#F0E442", "#0072B2","#D55E00", "#CC79A7")
#
# ceaclm <- melt(ceacl,id=c('iso3','threshold'))
# ceaclm[,Intervention:=ifelse(variable=='int','Intervention','SOC')]
#
# ## name key
# ckey <- data.table(iso3=c('CMR','UGA'),
#                    country=c('Cameroon','Uganda'))
#
# ceaclm <- merge(ceaclm,ckey,by='iso3',all.x=TRUE)
#
# ## plot
# GP <- ggplot(ceaclm,aes(threshold,value,
#                         col=country,lty=Intervention)) +
#   geom_line() +
#   theme_classic() +
#   theme(legend.position = 'top',legend.title = element_blank())+
#   ggpubr::grids()+
#   ylab('Probability cost-effective')+
#   xlab('Cost-effectiveness threshold (USD/DALY)')+
#   scale_colour_manual(values=cbPalette)
# GP
#
# ggsave(GP,file=gh('plots/CEACY') + SA + '.png',w=7,h=5)
#
#
# ## plot
# GP <- ggplot(ceaclm[variable=='int'],aes(threshold,value,
#                                          col=country)) +
#   geom_line() +
#   theme_classic() +
#   theme(legend.position = 'top',legend.title = element_blank())+
#   ggpubr::grids()+
#   ylab('Probability cost-effective')+
#   xlab('Cost-effectiveness threshold (USD/DALY)')+
#   scale_colour_manual(values=cbPalette) ## + xlim(x=c(0,1500))
# GP
#
# ggsave(GP,file=gh('plots/CEAC1Y') + SA + '.png',w=7,h=5)
#
# ## --- 5-14 YEARS
#
# ## containers & loop
# allout <- allpout <- list() #tabular outputs
# allout2 <- allpout2 <- list() #tabular outputs
# ceacl <- NMB <- list()             #CEAC outputs etc
# psaout <- psapout <- list()
# ## NOTE I think there was an additional problem here -
# ## because countries are contained as separate rows we were summing over both in the below computations
# ## ie we were looping but then operating over data for both countries
# ## cn <- isoz[1]
#
# for(cn in isoz){
#   dc <- D[isoz==cn & age=='5-14']
#   cat('running model for:',cn,'5-14 years\n')
#   # cat('running model for:',cn, dc$age[1], 'years \n')
#   ## --- costs
#   ## drop previous costs
#   # D <- D[age=='0-4',]
#   dc[,c(cnmz):=NULL]
#   ## add cost data
#   C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
#   dc <- merge(dc,C,c('id'),all.x=TRUE)        #merge into PSA
#   ## --- DALYs
#   ## drop any that are there
#   if('LYS' %in% names(D)) dc[,c('LYS','LYS0'):=NULL]
#   dc <- merge(dc,LYKc[iso3==cn,.(age,LYS,LYS0)],by='age',all.x = TRUE)        #merge into PSA
#   ## --- run model (quietly)
#   invisible(capture.output(dc <- runallfuns(dc,arm=arms)))
#   ## --- grather outcomes
#   out <- dc[,..toget]
#   out[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
#                    LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
#   ## out[,sum(value),by=id]                                       #CHECK
#   ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
#   out <- out[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum,by=id] #sum against popn
#   ## non-incremental cost per ATT
#   out[,costperATT.soc:=cost.soc/att.soc];
#   out[,costperATT.int:=cost.int/att.int];
#   out[,costperTPT.soc:=cost.soc/tpt.soc];
#   out[,costperTPT.int:=cost.int/tpt.int];
#   ## increments wrt SOC (per child screened at either facility/household)
#   out[,Dcost.int:=cost.int-cost.soc];  #inc costs
#   out[,Datt.int:=att.int-att.soc];  #inc atts
#   out[,Dtpt.int:=tpt.int-tpt.soc];  #inc atts
#   out[,Ddeaths.int:=deaths.int-deaths.soc];  #inc deaths
#   out[,DLYL0.int:=LYL0.int-LYL0.soc]; #inc LYLs w/o discount
#   out[,DLYL.int:=LYL.int-LYL.soc]; #inc LYLs
#   ## per whatever
#   out[,DcostperATT.int:=Dcost.int/Datt.int];
#   out[,DcostperTPT.int:=Dcost.int/Dtpt.int];
#   out[,Dcostperdeaths.int:=-Dcost.int/Ddeaths.int];
#   out[,DcostperLYS0.int:=-Dcost.int/DLYL0.int];
#   out[,DcostperLYS.int:=-Dcost.int/DLYL.int];
#   ## summarize
#   smy <- outsummary(out)
#   outs <- smy$outs; pouts <- smy$pouts;
#   outs[,iso3:=cn]; pouts[,iso3:=cn]
#   ## capture tabular
#   allout[[cn]] <- outs; allpout[[cn]] <- pouts
#   ## capture data for NMB
#   NMB[[cn]] <- out[,.(iso3=cn,DLYL.int,Dcost.int)]
#   ## ceac data
#   ceacl[[cn]] <- data.table(iso3=cn,
#                             int=make.ceac(out[,.(Q=-DLYL.int,P=Dcost.int)],lz),
#                             threshold=lz)
#   ## --- grather outcomes for Table 2
#   out2 <- dc[,..toget2]
#   out2[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
#                     LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
#   ## out[,sum(value),by=id]                                       #CHECK
#   ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
#   cntcts <- out2[,sum(value),by=id]
#   out2 <- out2[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum2,by=id] #sum against popn
#   out2[,contacts.int:=cntcts$V1]
#   # out2[,contacts.soc:=2.0] #TODO needs changing
#   out2[,contacts.soc:=ifelse(cn=='CMR', 1.202703, 1.941176)] # by country TODO needs checking
#   out2[,c('LYL0.soc','LYL0.int'):=NULL] #drop
#   out2[,c('prevdeaths.soc','prevdeaths.int'):=.(deaths.soc-incdeaths.soc,deaths.int-incdeaths.int)]
#   ## increments
#   out2[,Dtpt:=tpt.int-tpt.soc]
#   out2[,Datt:=att.int-att.soc]
#   out2[,DLYL:=LYL.int-LYL.soc]
#   out2[,Dprevtb:=prevtb.int-prevtb.soc]
#   out2[,Dinctb:=inctb.int-inctb.soc]
#   out2[,Dincdeaths:=incdeaths.int-incdeaths.soc]
#   out2[,Dprevdeaths:=prevdeaths.int-prevdeaths.soc]
#   out2[,Ddeaths:=deaths.int-deaths.soc]
#   out2[,Dcontacts:=contacts.int-contacts.soc]
#   # out2[,cost.att.soc:=cost.soc/att.soc];
#   # out2[,cost.att.int:=cost.int/att.int];
#   # out2[,cost.tpt.soc:=cost.soc/tpt.soc];
#   # out2[,cost.tpt.int:=cost.int/tpt.int];
#   out2[,Dcost.screen:=cost.screen.int-cost.screen.soc]
#   out2[,Dcost.tpt:=cost.tpt.int-cost.tpt.soc]
#   out2[,Dcost.prev.att:=cost.prev.att.int-cost.prev.att.soc]
#   out2[,Dcost.inc.att:=cost.inc.att.int-cost.inc.att.soc]
#   out2[,Dcost:=cost.int-cost.soc]
#
#   ## summarize
#   smy2 <- Table2(out2) #NOTE set per 1000 index cases or HHs - adjust fac in contact_functions.R
#   outs2 <- smy2$outs; pouts2 <- smy2$pouts;
#   outs2[,iso3:=cn]; pouts2[,iso3:=cn]
#   ## capture tabular
#   allout2[[cn]] <- outs2; allpout2[[cn]] <- pouts2
#
#   psaout[[cn]] <- dc[,.(iso3=cn,
#                         cost.SOC=sum(cost.soc*value),
#                         cost.INT=sum(cost.int*value),
#                         lyl.SOC=sum(deaths.soc*value*LYS),
#                         lyl.INT=sum(deaths.int*value*LYS)),
#                      by=id] #PSA summary
#
# }
#
#
# # summary(D[,..toget])
# # D[,.(isoz,age,F.u5,value)] # why are some values 0? TODO BUG ?
# # D[,summary(tpt.int/tpt.soc)]
# #
# # check1 <- names(labz)[4:14]
# # check1 <- c(paste0(check1,'.soc'), paste0(check1,'.int'))
# # check2 <- names(labz)[14:21]
# # summary(D[,..check1])
# # odds <- names(D)[grepl('OR',names(D))]
# # summary(D[,..odds])
#
# allout <- rbindlist(allout)
# allpout <- rbindlist(allpout)
# allout2 <- rbindlist(allout2)
# allpout2 <- rbindlist(allpout2)
# psaout <- rbindlist(psaout)
#
# ceacl <- rbindlist(ceacl)
# NMB <- rbindlist(NMB)
#
# fwrite(allout,file=gh('outdata/alloutO') + SA + '.csv')
# fwrite(allpout,file=gh('outdata/allpoutO') + SA + '.csv')
# fwrite(allout2,file=gh('outdata/allout2O') + SA + '.csv')
# fwrite(allpout2,file=gh('outdata/allpout2O') + SA + '.csv')
# save(ceacl,file=gh('outdata/ceaclY') + SA + '.Rdata')
# save(NMB,file=gh('outdata/NMBY') + SA + '.Rdata')
#
#
# ICERSO <- allpout[,.(iso3,ICER=ICER.int)]
# fwrite(ICERSO,file=gh('outdata/ICERSO') + SA + '.csv')
#
# ## CEAC plot
# cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73",
#                         "#F0E442", "#0072B2","#D55E00", "#CC79A7")
#
# ceaclm <- melt(ceacl,id=c('iso3','threshold'))
# ceaclm[,Intervention:=ifelse(variable=='int','Intervention','SOC')]
#
# ## name key
# ckey <- data.table(iso3=c('CMR','UGA'),
#                    country=c('Cameroon','Uganda'))
#
# ceaclm <- merge(ceaclm,ckey,by='iso3',all.x=TRUE)
#
# ## plot
# GP <- ggplot(ceaclm,aes(threshold,value,
#                         col=country,lty=Intervention)) +
#   geom_line() +
#   theme_classic() +
#   theme(legend.position = 'top',legend.title = element_blank())+
#   ggpubr::grids()+
#   ylab('Probability cost-effective')+
#   xlab('Cost-effectiveness threshold (USD/DALY)')+
#   scale_colour_manual(values=cbPalette)
# GP
#
# ggsave(GP,file=gh('plots/CEACO') + SA + '.png',w=7,h=5)
#
#
# ## plot
# GP <- ggplot(ceaclm[variable=='int'],aes(threshold,value,
#                                          col=country)) +
#   geom_line() +
#   theme_classic() +
#   theme(legend.position = 'top',legend.title = element_blank())+
#   ggpubr::grids()+
#   ylab('Probability cost-effective')+
#   xlab('Cost-effectiveness threshold (USD/DALY)')+
#   scale_colour_manual(values=cbPalette) ## + xlim(x=c(0,1500))
# GP
#
# ggsave(GP,file=gh('plots/CEAC1O') + SA + '.png',w=7,h=5)
ICERS
