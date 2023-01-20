library(here)
library(data.table)
library(glue)
library(googlesheets4)


## NOTE new sheet
## create an ID to access the googlesheets results sheet
yourl <- "https://docs.google.com/spreadsheets/d/1OultOF_rkxnaXM_Jnu81yuRSvVB2rg3R8dNi1uKJsWg/edit#gid=203919644"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(basename,filename,sheetid
                             ){
  filename <- gsub("\\.csv$","",filename) #safety in case csv included at and
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file=fn)
  write_sheet(tmp,sheetid,sheet=filename)
}

# ## upload relevant table data
# upload.to.sheets(here('model/outdata/'),'cascadetab',shid) #first will need to re-authenticate
# 
# ## rest can be run as block
# flz1 <- c(
#  "ACFcascade.1.csv",   "ACFcascadecdr.1.csv","ACFcascadehi.1.csv",
#  "ACFcascadelo.1.csv", "ACFcascadetxd.1.csv","cascadetab.csv",
#  "CEAC50.csv",         "CEAC50cdr.csv",      "CEAC50hi.csv",
#  "CEAC50lo.csv",       "CEAC50pt.1.csv",     "CEAC50ptcdr.1.csv",
#  "CEAC50pthi.1.csv",   "CEAC50ptlo.1.csv",   "CEAC50pttxd.1.csv",
#  "CEAC50txd.csv",      "ICERagept.1.csv",    "ICERageptcdr.1.csv",
#  "ICERagepthi.1.csv",  "ICERageptlo.1.csv",  "ICERagepttxd.1.csv",
#  "ICERall.1.csv",      "ICERallcdr.1.csv",   "ICERallhi.1.csv",
#  "ICERalllo.1.csv",    "ICERalltxd.1.csv")
# for( fn in flz1)
#   upload.to.sheets(here('model/outdata/'),fn,shid)
# 
# Sys.sleep(120) #wait a bit so as not to annoy google
# 
# flz2 <- c("ICERatt.csv",
#  "ICERattcdr.csv",     "ICERatthi.csv",      "ICERattlo.csv",
#  "ICERatttxd.csv",     "ICERpt.1.csv",       "ICERptcdr.1.csv",
#  "ICERpthi.1.csv",     "ICERptlo.1.csv",     "ICERpttxd.1.csv",
#  "ICERSatt.csv",       "ICERSattcdr.csv",    "ICERSatthi.csv",
#  "ICERSattlo.csv",     "ICERSatttxd.csv",    "PTC.csv",
#  "ptsuccess.csv",      "txsuccess.csv"
# )
# for( fn in flz2)
#   upload.to.sheets(here('model/outdata/'),fn,shid)


## need article tables
yurl <- "https://docs.google.com/spreadsheets/d/1vRO20Qjuek5nwlJTUf2E5aIVyWV8Hwj-YzX2JW5kNwc/edit#gid=8659339"
shidneat <- as.character(as_sheets_id(yurl))

## ---- Table 1 -------
## build 1st
# load(here('model/data/Table1ATT.Rdata'))
# load(here('model/data/Table1PT.Rdata'))
# load(here('model/data/Table1PTcost.Rdata'))
# 
# setcolorder(Table1PTcost,names(Table1PT))
# 
Table1 <- fread(gh('outdata/table1.csv')) 

write_sheet(Table1,shidneat,sheet="Tab1RAW")

# upload.to.sheets(here('outdata//'),'table1',shidneat) #first will need to re-authenticate

## ---- Table 2 -------

## 0-14 years 
Table2 <- fread(gh('outdata/table2.csv')) 

write_sheet(Table2,shidneat,sheet="Tab2RAW")

## 0-4 years 
Table2Y <- fread(gh('outdata/table2Y.csv')) 

write_sheet(Table2Y,shidneat,sheet="Tab2YRAW")

## 5-14 years 
Table2O <- fread(gh('outdata/table2O.csv')) 

write_sheet(Table2O,shidneat,sheet="Tab2ORAW")

## --- SA table ---
sa.base <- fread(here('outdata/ICERS.csv'))
sa.hi <- fread(here('outdata/ICERShi.csv'))
sa.lo <- fread(here('outdata/ICERSlo.csv'))
sa.tptru <- fread(here('outdata/ICERStptru.csv'))
sa.hicoprev <- fread(here('outdata/ICERShicoprev.csv'))

names(sa.base)[2] <- 'Base case'
names(sa.hi)[2] <- '5% discount rate'
names(sa.lo)[2] <- '0% discount rate'
names(sa.tptru)[2] <- 'Baseline resource use (TPT visits)'
names(sa.hicoprev)[2] <- 'Higher co-prevalence of tuberculosis disease'

SAll <- Reduce(merge,list(sa.base,sa.hi,sa.lo,sa.tptru,sa.hicoprev))

write_sheet(SAll,shidneat,sheet="SAll")

## --- Main parameters table ---
ParmsTable1 <- fread(gh('outdata/parameters1.csv')) 

write_sheet(ParmsTable1,shidneat,sheet="ParmsTab1RAW")

# Other parameter tables
flz1 <- c(
 "tableS1.csv",   "tableS2.csv","tableS2a.csv",
 "tableS3.csv", "tableS4.csv")
for( fn in flz1)
  upload.to.sheets(here('outdata/'),fn,shidneat)






