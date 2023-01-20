## Generate manuscript tables/results for CONTACT
# rm(list=ls())
library(here)
library(tidyverse)

## load other scripts
# source(here('R/contact_run.R'))           # call if there is need to generate results otherwise they're already available

cost_summary <- fread(gh('indata/mean_cost_summary.csv')) 
# cost_summary <- fread(here('indata/mean_cost_summary.csv'))

cost_summary <- cbind(cost_summary[country=='Cameroon',.('Cascade step'=cascade_cost, Control=SOC, Intervention=INT)], 
                      cost_summary[country=='Uganda', .(Control=SOC, Intervention=INT)])

cost_summary$`Cascade step` <- c('TB symptom screening', 'Prevalent TB investigations', 'Prevalent TB treatment', 
                                 'TB preventive therapy', 'Incident TB investigations','Incident TB treatment')

fwrite(cost_summary,file=here('outdata/table1.csv'))

# Table 2: main results
table2 <- fread(gh('outdata/allpout2.csv')) 
table2 <- melt(table2, id.vars = 'iso3')
table2$variable[grepl('.soc', table2$variable)] 

soc <- table2[grepl('.soc', variable),]
soc <- soc[,variable:=gsub('.soc', '', variable)]
names(soc)[3] <- 'Control'

int <- table2[grepl('.int', variable),]
int <- int[!grepl('ICER', variable),]
int <- int[,variable:=gsub('.int', '', variable)]
names(int)[3] <- 'Intervention'

d <- table2[grepl('D|ICER', variable),]
d <- d[,variable:=gsub('D', '', variable)]
names(d)[3] <- 'Increment'

merged <- merge(soc, int, by=c('iso3', 'variable'), all.x = TRUE, all.y = TRUE)
merged <- merge(merged, d,by=c('iso3', 'variable'), all.x = TRUE, all.y = TRUE)

vars <- c('contacts','tpt','att',
          'prevtb','inctb','prevdeaths','incdeaths','deaths','LYL',
          'cost.screen','cost.tpt','cost.prev.att','cost.inc.att','cost', 'ICER.int')
var_labs <- c('Household contacts screened','TPT courses','ATT courses',
              'Prevalent TB','Incident TB','Prevalent TB deaths','Incident TB deaths','Total deaths','Discounted LYL',
              'Screening cost','TPT cost','Prevalent ATT cost','Incident TPT cost','Total cost','ICER')
merged$variable <- factor(merged$variable, 
                          levels = vars)   
setorder(merged, iso3, variable)

table2_format <- cbind(merged[iso3=='CMR',.(variable, Control, Intervention, Increment)], merged[iso3=='UGA', .(Control, Intervention, Increment)])
table2_format$variable <- var_labs
table2_format
fwrite(table2_format,file=here('outdata/table2.csv'))

# Supplementary Table: Outcomes parameters
parameters1 <- setDT(PD0[2:24,])
parameters1 <- parameters1[!NAME %in% c('hivprev.u5', 'hivprev.o5'),]
parameters1 <- parameters1[,.(NAME, DISTRIBUTION, 'MEAN (IQR)'=`MEDIAN..IQR.`, DESCRIPTION, SOURCE)]

fwrite(parameters1,file=here('outdata/tableS1.csv'))

# Supplementary Table: Cascade of care parameters (aggregated)
# Need ParmsTab2, 
# ParmsTab2 <- fread(gh('outdata/Parameters2.csv')) 
ParmsTab2 <- fread(here('indata/proportions.csv'))

ParmsTab2 <- ParmsTab2[,variable:=gsub('int.|soc.', '',variable)]
# ParmsTab2 <- dcast(ParmsTab2, country+isoz+variable~care_model)
unique(ParmsTab2$variable)
ParmsTab2 <- ParmsTab2[!variable %in% c("frac.bac.noclin.dx" ,"frac.bac.7d.ass","frac.bac.7d.noclin.dx","frac.bac.7d.ltfu","frac.noclin.dx","frac.clin.7d.ass","frac.noclin.7d.dx","frac.clin.7d.ltfu","frac.tpt.not.assessed","frac.tpt.not.reassessed","frac.tpt.not.eligible","frac.tpt.not.reeligible",  "frac.pre.tpt.ltfu" ,"frac.pre.tpt.reltfu","frac.tpt.notcompleted")]
unique(ParmsTab2$variable)
length(unique(ParmsTab2$variable))
vars <- c('frac.screened','frac.asymp','frac.symp',
          'frac.need.rescr','frac.rescreened','frac.rescr.asymp','frac.rescr.symp',
          'frac.tpt.assessed','frac.tpt.eligible','frac.tpt.initiated','frac.tpt.reassessed','frac.tpt.reeligible','frac.tpt.reinitiated','frac.tpt.completed',
          'frac.symp.attending', "frac.tb.evaluation",
          'frac.bac.assess','frac.bac.dx','frac.bac.clin.dx','frac.bac.7d.clin.dx','frac.clin.assess','frac.clin.dx',
          'frac.clin.7d.dx')

var_labs <- c('Screened for TB symptoms','Negative screening','Positive screening',
              'Need re-screening','Re-screened for TB symptoms','Negative re-screening','Positive re-screening',
              'Assessed for TPT','Eligible for TPT','Initiated on TPT',
              'Re-assessed for TPT','Re-eligible for TPT','Re-initiated on TPT','TPT completed',
              'Attending facility referral',"Evaluated for TB",
              'Bacteriological assessment','Bacteriological diagnosis','Clinical diagnosis after bacteriological assessment','Clinical diagnosis after bacteriological assessment & 7 days',
              'Clinical diagnosis without bacteriological assessment',
              'Clinical diagnosis after bacteriological assessment & 7 days','Clinical diagnosis without bacteriological assessment & 7 days')

ParmsTab2$variable <- factor(ParmsTab2$variable, 
                          levels = vars)   
setorder(ParmsTab2, country, variable)

ParmsTab2 <- cbind(ParmsTab2[country=='Cameroon',.(variable, Control=SOC, Intervention=INT)], ParmsTab2[country=='Uganda', .(Control=SOC, Intervention=INT)])
ParmsTab2$variable <- var_labs
ParmsTab2
fwrite(ParmsTab2,file=here('outdata/tableS2.csv'))

# Supplementary Table: Cascade of care parameters (disaggregated)
# Need ParmsTab2, 
ParmsTab2a <- fread(gh('outdata/Parameters2a.csv')) 
ParmsTab2a <- ParmsTab2a[,variable:=gsub('int.|soc.', '',variable)]
ParmsTab2a <- dcast(ParmsTab2a, country+isoz+variable+age~care_model)
unique(ParmsTab2a$variable)
ParmsTab2a <- ParmsTab2a[!variable %in% c('frac.bac.7d.noclin.dx','frac.bac.noclin.dx','frac.clin.7d.noclin.dx','frac.noclin.dx')]
vars <- c('frac.screened','frac.rescreened',
          'frac.asymp','frac.tpt.eligible','frac.tpt.initiated','frac.tpt.completed',
          'frac.symp','frac.rescr.symp','frac.symp.attending',
          'frac.bac.assess','frac.bac.dx','frac.bac.clin.dx','frac.clin.dx',
          'frac.bac.7d.clin.dx','frac.clin.7d.clin.dx')

var_labs <- c('Screened for TB symptoms','Re-screened for TB symptoms',
              'Negative for TB symptoms','frac.tpt.eligible','frac.tpt.initiated','frac.tpt.completed',
              'Positive for TB symptoms','Positive for TB symptoms (@re-screening)','Attending facility referral',
              'Bacteriological assessment','Bacteriological diagnosis','Clinical diagnosis after bacteriological assessment',
              'Clinical diagnosis without bacteriological assessment',
              'Clinical diagnosis after bacteriological assessment & 7 days','Clinical diagnosis without bacteriological assessment & 7 days')

ParmsTab2a$variable <- factor(ParmsTab2a$variable, 
                             levels = vars)   
setorder(ParmsTab2a, country, isoz, variable)

ParmsTab2a <- cbind(ParmsTab2a[isoz=='CMR',.(age, variable, Control=SOC, Intervention=INT)], ParmsTab2a[isoz=='UGA', .(Control=SOC, Intervention=INT)])
ParmsTab2a$variable <- rep(var_labs, each=2)
ParmsTab2a
fwrite(ParmsTab2a,file=here('outdata/tableS2a.csv'))

# Supplementary Table: Cascade of care parameters (disaggregated)
# Need DENR, 
DENR <- data.table(read_excel(here("indata/Baseline information.xlsx"), sheet = 'Sheet1', range = 'O3:R19'))
DENR <- melt(DENR, variable.name = 'isoz')
DENR <- DENR[,country:=ifelse(isoz=='CMR', 'Cameroon', 'Uganda')]

TBPREV <- fread(here("indata/tb_epi.csv"))
TBPREV[, isoz:=ifelse(country=='Uganda', 'UGA', 'CMR')]
TBPREV <- TBPREV[,model:='int']
TBPREV <- melt(TBPREV, id.vars = c('country','isoz','model'), variable.name = 'metric')

ParmsTab3 <- rbind(DENR, TBPREV)
ParmsTab3 <- ParmsTab3[, model:=ifelse(model=='soc', 'Control', 'Intervention')]
ParmsTab3 <- dcast(ParmsTab3, country+isoz+metric~model)
ParmsTab3 <- ParmsTab3[!metric %in% c('declared_per_household', 'enrolled_per_household', 'frac.declared.u5', 
                                      'frac.declared.hiv')]
unique(ParmsTab3$metric)
vars <- c('declared_per_index_case','enrolled_per_index_case','frac.enrolled.u5','frac.enrolled.hiv','tb_prev','tb_inc')
var_labs <- c('Child contacts declared per index case','Child contacts enrolled per index case','Child contacts enrolled under 5 years','Child contacts enrolled HIV+','TB co-prevalence intervention','TB incidence intervention')

ParmsTab3$metric <- factor(ParmsTab3$metric, 
                              levels = vars)   
setorder(ParmsTab3, country, isoz, metric)

ParmsTab3 <- cbind(ParmsTab3[isoz=='CMR',.(country, metric, Control, Intervention)], 
                   ParmsTab3[isoz=='UGA', .(Control, Intervention)])
ParmsTab3$metric <- var_labs
ParmsTab3
fwrite(ParmsTab3,file=here('outdata/tableS3.csv'))

# Supplementary Table: intervention effects 
# Need Parameters3, 
ParmsTab4 <- fread(gh('outdata/Parameters3.csv')) 
ParmsTab4 <- ParmsTab4[NAME %in% c('tpt.resultOR','tpt.initiationOR','tpt.completionOR','tb.resultOR','tb.diagnosisOR'),]
ParmsTab4 <- ParmsTab4[,.(NAME, DISTRIBUTION, `MEDIAN (IQR)`, DESCRIPTION, SOURCE)]
ParmsTab4$NAME <- c("negative.screening", "tpt.initiation", "tpt.completion", "positive.screening", "tb.diagnosis"  )
ParmsTab4
fwrite(ParmsTab4,file=here('outdata/tableS4.csv'))

# Supplementary Results Table: 0-4 years results
# Table 2: supplementary results
table2Y <- fread(gh('outdata/allpout2Y.csv')) 
table2Y <- melt(table2Y, id.vars = 'iso3')
table2Y$variable[grepl('.soc', table2Y$variable)] 

soc <- table2Y[grepl('.soc', variable),]
soc <- soc[,variable:=gsub('.soc', '', variable)]
names(soc)[3] <- 'Control'

int <- table2Y[grepl('.int', variable),]
int <- int[!grepl('ICER', variable),]
int <- int[,variable:=gsub('.int', '', variable)]
names(int)[3] <- 'Intervention'

d <- table2Y[grepl('D|ICER', variable),]
d <- d[,variable:=gsub('D', '', variable)]
names(d)[3] <- 'Increment'

merged <- merge(soc, int, by=c('iso3', 'variable'), all.x = TRUE, all.y = TRUE)
merged <- merge(merged, d,by=c('iso3', 'variable'), all.x = TRUE, all.y = TRUE)

vars <- c('contacts','tpt','att',
          'prevtb','inctb','prevdeaths','incdeaths','deaths','LYL',
          'cost.screen','cost.tpt','cost.prev.att','cost.inc.att','cost', 'ICER.int')
var_labs <- c('Household contacts screened','TPT courses','ATT courses',
              'Prevalent TB','Incident TB','Prevalent TB deaths','Incident TB deaths','Total deaths','Discounted LYL',
              'Screening cost','TPT cost','Prevalent ATT cost','Incident TPT cost','Total cost','ICER')
merged$variable <- factor(merged$variable, 
                          levels = vars)   
setorder(merged, iso3, variable)

table2Y_format <- cbind(merged[iso3=='CMR',.(variable, Control, Intervention, Increment)], merged[iso3=='UGA', .(Control, Intervention, Increment)])
table2Y_format$variable <- var_labs
table2Y_format
fwrite(table2Y_format,file=here('outdata/table2Y.csv'))

# Supplementary Table: 5-14 years results
# Table 2: main results
table2O <- fread(gh('outdata/allpout2O.csv')) 
table2O <- melt(table2O, id.vars = 'iso3')
table2O$variable[grepl('.soc', table2O$variable)] 

soc <- table2O[grepl('.soc', variable),]
soc <- soc[,variable:=gsub('.soc', '', variable)]
names(soc)[3] <- 'Control'

int <- table2O[grepl('.int', variable),]
int <- int[!grepl('ICER', variable),]
int <- int[,variable:=gsub('.int', '', variable)]
names(int)[3] <- 'Intervention'

d <- table2O[grepl('D|ICER', variable),]
d <- d[,variable:=gsub('D', '', variable)]
names(d)[3] <- 'Increment'

merged <- merge(soc, int, by=c('iso3', 'variable'), all.x = TRUE, all.y = TRUE)
merged <- merge(merged, d,by=c('iso3', 'variable'), all.x = TRUE, all.y = TRUE)

vars <- c('contacts','tpt','att',
          'prevtb','inctb','prevdeaths','incdeaths','deaths','LYL',
          'cost.screen','cost.tpt','cost.prev.att','cost.inc.att','cost', 'ICER.int')
var_labs <- c('Household contacts screened','TPT courses','ATT courses',
              'Prevalent TB','Incident TB','Prevalent TB deaths','Incident TB deaths','Total deaths','Discounted LYL',
              'Screening cost','TPT cost','Prevalent ATT cost','Incident TPT cost','Total cost','ICER')
merged$variable <- factor(merged$variable, 
                          levels = vars)   
setorder(merged, iso3, variable)

table2O_format <- cbind(merged[iso3=='CMR',.(variable, Control, Intervention, Increment)], merged[iso3=='UGA', .(Control, Intervention, Increment)])
table2O_format$variable <- var_labs
table2O_format
fwrite(table2O_format,file=here('outdata/table2O.csv'))

#  All below is work in progress
# Looking at CE plane & other metrics
## ## scraps for development of below fn
## data <- dc
## Kmax <- 1e3
## file.id <- 'test'
## wtp <- 500

## NOTE this is more illustrative for now
## NOTE needs a folder called graphs/ creating (which is currently excluded from the repo)
## some automatic CEA outputs
file.id='KEN';Kmax=5e3;wtp=5e3;
MakeCEAoutputs <- function(data,LY,
                           file.id='',Kmax=5e3,wtp=5e3,
                           arms=c('SOC','INT')){
  data <- merge(D,LYK,by='age') #add age
  DS <- d[isoz=='CMR',.(cost.SOC=sum(cost.soc*value),
                        cost.INT=sum(cost.int*value),
                        lyl.SOC=sum(deaths.soc*value*LYS),
                        lyl.INT=sum(deaths.int*value*LYS)),
          by=id] #PSA summary
  
  ## prep for BCEA
  DS <- data.frame(psaout[iso3=='CMR',])
  LYS <- CST <- matrix(nrow=nreps,ncol=2)
  LYS[,1] <- 1-DS$lyl.SOC #NOTE this is life years lost
  LYS[,2] <- 1-DS$lyl.INT
  CST[,1] <- DS$cost.SOC
  CST[,2] <- DS$cost.INT
  ## BCEA outputs
  M <- bcea(e=LYS,c=CST,interventions = arms,Kmax=Kmax)
  print(summary(M))
  
  M$k[M$ceac>=0.5][1]
  print(with(M,ceac[which(k==830)]))
  print(with(M,ceac[which(k==1660) ]))
  
  print(with(M,ceac[which(k==430) ]))
  print(with(M,ceac[which(k==860) ]))
  
  fn <- paste0(here('outdata/kstar_'),file.id,'.txt')
  cat(M$kstar,file = fn)
  fn <- paste0(here('outdata/ICER_'),file.id,'.txt')
  cat(M$ICER,file = fn)
  
  ## NOTE may need more configuration
  ceac.plot(M,graph='ggplot2') +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/CEAC_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)
  
  ceplane.plot(M,graph='ggplot2',wtp=wtp,title = "")+
    scale_x_continuous(label=comma) +
    theme_classic() +
    theme(legend.position = 'top') + ggpubr::grids()
  fn <- paste0(here('graphs/CE_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)
  
  eib.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EIB_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)
  
  evi.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EVI_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)
  
}

psaout[[cn]] <- dc[,.(iso3=cn,
                      cost.SOC=sum(cost.soc*value),
                      cost.INT=sum(cost.int*value),
                      lyl.SOC=sum(deaths.soc*value*LYS),
                      lyl.INT=sum(deaths.int*value*LYS)),
                   by=id] #PSA summary
