library(here)
library(data.table)
library(ggplot2)

load(here('outdata/out.Rdata'))

effects <- out[,.(tpt.ratio=tpt.int/tpt.soc,
                  att.ratio=att.int/att.soc,
                  deaths.ratio=deaths.int/deaths.soc,
                  Ddeaths=deaths.int-deaths.soc,
                  DLYL=LYL.int-LYL.soc,
                  Dcost=cost.int-cost.soc)]

summary(effects)

## text:
## 35% increase in tpt
## 15% reduction in ATT
## 22% reduction in deaths
## diff: 1% of LYL

## not much death
