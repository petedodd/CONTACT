## top level handler for CONTACT
rm(list=ls())
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
  LYKc <- merge(LYKc,LYKc0[,.(iso3,age,LYS0=LYS)],by=c('iso3','age'))
  LYK <- LYKc[,.(LYS=mean(LYS),LYS0=mean(LYS0)),by=.(age)] #averaged life-years 4 generic tests
  save(LYKc,file=here('indata/LYKc.Rdata'))
  save(LYK,file=here('indata/LYK.Rdata'))
} else {
  load(file=here('indata/LYKc.Rdata'))
  load(file=here('indata/LYK.Rdata'))
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

tmp1 <- getLNparms(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2,med=FALSE)
tmp[,DISTRIBUTION:=paste0("LN(",tmp1$mu,",",tmp1$sig,")")] #LN distributions

INTE[,DISTRIBUTION:=tmp$DISTRIBUTION]
INTE[DISTRIBUTION=='LN(NA,NA)', DISTRIBUTION:=NA]

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

PP <- dcast(PP, isoz+age~variable)
# 

# load enrollment data (number of contacts per index case & fraction under 5 years)
# also available but not currently being used: declared per household/enrolled, enrolled per household, frac declared u5, frac declared/enrolled hiv
library(readxl)
base_popn <- data.table(read_excel("~/Dropbox/CONTACT/indata/Baseline information .xlsx", sheet = 'Sheet1', range = 'Z13:AB17'))
base_popn <- melt(base_popn, variable.name = 'isoz')
age_split <- base_popn[grepl('enrolled.', metric), .(isoz, variable=metric, value)]
age_split <- age_split[isoz=='CMR', variable:=paste('cmr.', variable, sep = '')]
age_split <- age_split[isoz=='UGA', variable:=paste('uga.', variable, sep = '')]

tmp <- data.table(matrix(NA, nrow = nrow(age_split), ncol = ncol(PD0)))
names(tmp) <- names(PD0)
tmp <- tmp[,NAME:=age_split$variable]
tmp <- tmp[,DISTRIBUTION:=age_split$value]
age_split <- copy(tmp)
age_split <- age_split[,NAME:=gsub('enrolled.', '',NAME)]

# commented CODE below not working
enrolled <- data.table(read_excel("~/Dropbox/CONTACT/indata/Baseline information .xlsx", sheet = 'Sheet1', range = 'O3:R19'))

enrolled <- melt(enrolled, variable.name = 'isoz')

enrolled <- enrolled[,model:=toupper(model)]
enrolled <- enrolled[model=='INT', metric:=paste('int.', metric, sep = '')]
enrolled <- enrolled[model=='SOC', metric:=paste('soc.', metric, sep = '')]
enrolled <- enrolled[grepl('enrolled.per.index', metric), .(isoz, variable=metric, care_model=model, value)]

enrolled <- dcast(enrolled, isoz~variable)
PP <- merge(PP, enrolled, by='isoz') # age disaggregated
PPA <- merge(PPA, enrolled, by='isoz') # aggregated
PPA <- rbind(PPA,PPA)
PPA[,age:=rep(agelevels,each=nrow(PPA)/2)]

# age_split <- base_popn[grepl('enrolled.u5|enrolled.hiv', metric), .(isoz, variable=metric, care_model=model, value)]
# age_split <- age_split[isoz=='CMR', variable:=paste('cmr.', variable, sep = '')]
# age_split <- age_split[isoz=='UGA', variable:=paste('uga.', variable, sep = '')]
# tmp <- data.table(matrix(NA, nrow = nrow(age_split), ncol = ncol(PD0)))
# names(tmp) <- names(PD0)
# tmp <- tmp[,NAME:=age_split$variable]
# tmp <- tmp[,DISTRIBUTION:=age_split$value]
# age_split <- copy(tmp)
# age_split <- age_split[,NAME:=gsub('enrolled.', '',NAME)]


# # load cohort data (fraction under 5 and HIV)
# load(here('indata/cohortdata.Rdata'))
# 
# # fraction under 5 years
# popn <- setDT(popn_age_country)
# popn <- popn[country=='Cameroun', country:='Cameroon']
# popn <- popn[age_cat=='u5', variable:='frac.u5']
# popn <- popn[age_cat=='o5', variable:='frac.o5']
# popn[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
# popn[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]
# popn <- popn[variable=='frac.u5',.(isoz, age, variable, prop)]
# # popn <- popn[rando=='ITV', rando:='INT']
# popn <- dcast(popn, isoz~variable)
# 
# # HIV prevalence
# popn_hiv <- setDT(popn_hiv_age_country)
# popn_hiv <- popn_hiv[country=='Cameroun', country:='Cameroon']
# popn_hiv <- popn_hiv[hiv_status=='Positive', ]
# popn_hiv <- popn_hiv[age_cat=='u5', variable:='frac.hiv.u5']
# popn_hiv <- popn_hiv[age_cat=='o5', variable:='frac.hiv.o5']
# popn_hiv[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
# popn_hiv[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]
# popn_hiv <- popn_hiv[,.(isoz, age, variable, prop)]
# 
# tmp <- popn_hiv[isoz=='CMR',]
# tmp <- tmp[,age:='0-4']
# tmp <- tmp[,variable:='frac.hiv.u5']
# tmp <- tmp[,prop:=0]
# 
# popn_hiv <- rbind(popn_hiv, tmp)
# # popn <- popn[rando=='ITV', rando:='INT']
# popn_hiv <- dcast(popn_hiv, isoz~variable)
# 
# # ART coverage
# popn_art <- setDT(popn_art_age_country)
# popn_art <- popn_art[country=='Cameroun', country:='Cameroon']
# popn_art <- popn_art[art_status=='Yes', ]
# popn_art <- popn_art[age_cat=='u5', variable:='artcov.u5']
# popn_art <- popn_art[age_cat=='o5', variable:='artcov.o5']
# popn_art[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
# popn_art[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]
# popn_art <- popn_art[,.(isoz, age, variable, prop)]
# 
# # TODO: check under 5 art coverage - looks too low. Setting it to O5 coverage for now
# # popn_art <- popn_art[age=='0-4',prop:=prop[age=='5-14']] 
# popn_art$prop[popn_art$age=='0-4'] <- popn_art$prop[popn_art$age=='5-14']
# 
# 
# # popn <- popn[rando=='ITV', rando:='INT']
# popn_art <- dcast(popn_art, isoz~variable)
# 
# 
# PP <- merge(PP, popn, by='isoz')
# PP <- merge(PP, popn_hiv, by='isoz')
# PP <- merge(PP, popn_art, by='isoz')
# P2 <- melt(PP)


## combine different parameter types
# drop <- c(names(popn)[2], names(popn_hiv)[2:3])

# adding study data on: frac u5, HIV prevalence and ART coverage
# currently using CMR data
# TODO: figure out how to add country and model of care specific data
names(INTE) <- names(PD0)
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

# D <- rbind(D,D)
# D[,isoz:=rep(c('CMR','UGA'),each=nrow(D)/2)]
# D <- merge(D,PP,by=c('isoz', 'age'),all.x = TRUE)
D <- merge(D,PPA,by=c('isoz', 'age'),all.x = TRUE)
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

rcsts <- fread(gh('indata/model_mean_total_costs.csv'))    #read cost data
names(rcsts)

rcsts <- rcsts[rando=='INT', cascade:=paste('c.int.', cascade, sep = '')]
rcsts <- rcsts[rando=='SOC', cascade:=paste('c.soc.', cascade, sep = '')]

# rcsts[, age:=ifelse(age_cat=='u5', '0-4', '5-14')]
rcsts[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]

keeps <- c('isoz', 'rando', 'cascade','cost.m', 'cost.sd')
rcsts <- rcsts[,..keeps]

## turn cost data into PSA
rcsts[is.na(rcsts)] <- 0 #some quick fix >> setting NA to 0
rcsts[cost.sd==0,cost.sd:=cost.m/40]        #SD such that 95% UI ~ 10% of mean
# rcsts[cost.sd>100,cost.sd:=cost.m/40]
# rcsts[cost.m>100,cost.m:=cost.m/20]
allcosts <- rcsts[,.(iso3=isoz, cost=cascade, cost.m, cost.sd)]

keep <- vrz[grepl('c.soc.|c.int.', vrz)]

# check
setdiff(unique(allcosts$cost[allcosts$iso3=='UGA']),
        unique(keep))

setdiff(unique(keep),
        unique(allcosts$cost[allcosts$iso3=='UGA']))

allcosts[24,cost:='c.soc.inv_inc']
C <- MakeCostData(allcosts[iso3=='CMR'],nreps)               #make cost PSA NOTE using CMR cost data
## compute other parameters (adds by side-effect)
MakeTreeParms(D,P)

D[,sum(value),by=.(isoz,id)] #CHECK
D[,sum(value),by=.(isoz,id,age)] #CHECK
D[,.(isoz,age,F.u5,hivprev.u5)]
D[,.(isoz,age,soc.frac.screened,tb.screeningOR,int.frac.screened)]
## check for leaks
head(SOC.F$checkfun(D)) #SOC arm
head(INT.F$checkfun(D)) #INT arm

names(SOC.F)

## add cost data
D <- merge(D,C,by=c('id'),all.x=TRUE)       #merge into PSA

## run model
arms <- c('SOC','INT')
D <- runallfuns(D,arm=arms)                      #appends anwers

## --- run over different countries
cnmz <- names(C) # C=cost PSA data
cnmz <- cnmz[cnmz!=c('id')]
toget <- c('id',
           'cost.soc','cost.int',
           'tpt.soc','tpt.int',
           'att.soc','att.int',
           'deaths.soc','deaths.int',
           'LYS','LYS0','value'
           )
notwt <- c('id','LYS','LYS0','value') #variables not to weight against value
lyarm <- c('LYL.soc','LYL.int')
lyarm <- c(lyarm,gsub('\\.','0\\.',lyarm)) #include undiscounted
tosum <- c(setdiff(toget,notwt),lyarm)
## heuristic to scale top value for thresholds:
heur <- c('id','value','deaths.int','deaths.soc')
out <- D[,..heur]
out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=c('deaths.int','deaths.soc'),by=id] #sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.int)]
topl <- 10000*50
lz <- seq(from = 0,to=topl,length.out = 1000) #threshold vector for CEACs

## containers & loop
allout <- allpout <- list() #tabular outputs
ceacl <- NMB <- list()             #CEAC outputs etc
## cn <- isoz[1]
for(cn in isoz){
  cat('running model for:',cn,'\n')
  ## --- costs
  ## drop previous costs
  # D <- D[age=='0-4',]
  D[,c(cnmz):=NULL]
  ## add cost data
  C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
  D <- merge(D,C,c('id'),all.x=TRUE)        #merge into PSA
  ## --- DALYs
  ## drop any that are there
  if('LYS' %in% names(D)) D[,c('LYS','LYS0'):=NULL]
  D <- merge(D,LYKc[iso3==cn,.(age,LYS,LYS0)],by='age',all.x = TRUE)        #merge into PSA
  ## --- run model (quietly)
  invisible(capture.output(D <- runallfuns(D,arm=arms)))
  ## --- grather outcomes
  out <- D[,..toget]
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
}

D[,..toget]
D[,.(isoz,age,F.u5,value)] # why are some values 0?
D[,summary(tpt.int/tpt.soc)]

check1 <- names(labz)[4:14]
check1 <- c(paste0(check1,'.soc'), paste0(check1,'.int'))
check2 <- names(labz)[14:21]
summary(D[,..check1])
odds <- names(D)[grepl('OR',names(D))]
summary(D[,..odds])

allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
ceacl <- rbindlist(ceacl)
NMB <- rbindlist(NMB)

fwrite(allout,file=gh('outdata/allout.csv'))
fwrite(allpout,file=gh('outdata/allpout.csv'))
save(ceacl,file=gh('outdata/ceacl.Rdata'))
save(NMB,file=gh('outdata/NMB.Rdata'))


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

ggsave(GP,file=gh('plots/CEAC.png'),w=7,h=5)



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

ggsave(GP,file=gh('plot/CEAC1.png'),w=7,h=5)



## ## generate some CEA outputs in graphs/ & outdata/
## ## NOTE these folders need to be created
## ## NOTE need ggpubr, BCEA installed
# MakeCEAoutputs(D, #PSA dataset
#                LYK, #discounted expected life-years by age
#                file.id='test', #string to identify output files
#                Kmax=5e3,wtp=5e3)
