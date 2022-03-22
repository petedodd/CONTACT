## looking at tree 7
rm(list=ls())
library(here)
library(HEdtree)
library(data.tree)
library(data.table)
library(glue)

# set_here('C:/Users/cm1nm/Dropbox/CONTACT')

## === outcomes subtree ===
tb <- txt2tree(here('indata/DSTB.outcomes.txt')) # tb dx
notb <- txt2tree(here('indata/noTB.outcomes.txt')) # no tb
tpt <- txt2tree(here('indata/TPT.outcomes.txt')) # tpt

## default prob/cost:
tb$Set(p=1)
notb$Set(p=1)
tpt$Set(p=1)
tb$Set(cost=0)
notb$Set(cost=0)
tpt$Set(cost=0)

## set probabilities
## NOTE these namings actually get overwritten by CSV read-ins below
## -- tb outcomes:
tb$`no DS-TB tx`$p <- 'p.ptltfu'
tb$`no DS-TB tx`$dies$p <- 'p.cfr.notx'
tb$`no DS-TB tx`$survives$p <- '1-p.cfr.notx'
tb$`DS-TB tx`$p <- '1-p.ptltfu'
tb$`DS-TB tx`$dies$p <- 'p.cfr.tx'
tb$`DS-TB tx`$survives$p <- '1-p.cfr.tx'

## -- tpt outcomes:
tpt$`TPT tx`$p <- '1-p.ptptltfu'
tpt$`TPT tx`$`TB <1yr`$p <- 'p.tbdx.<1yr'
tpt$`TPT tx`$`no TB <1yr`$p <- '1-p.tbdx.<1yr'
tpt$`no TPT tx`$p <- 'p.ptptltfu'
tpt$`no TPT tx`$`TB <1yr`$p <- 'p.tbdx.<1yr'
tpt$`no TPT tx`$`no TB <1yr`$p <- '1-p.tbdx.<1yr'

print(tptx,'p')

# no tb outcomes
notb$dies$p <- 'p.cfr.notb'
notb$survives$p <- '1-p.cfr.notb'
print(notb,'p')

## ====== function to add outcomes & counters
AddOutcomes <- function(D){
  ## === cost and probs (defaults)
  D$Set(p=1)
  D$Set(cost=0)
  
  ## === merge to create final tree ===
  MergeByName(D,asymp,'Asymptomatic child contact management')
  MergeByName(D,symp,'Symptomatic child contact management')
  MergeByName(D,tpt,'TPT outcomes')
  MergeByName(D,tb,'TB outcomes')
  MergeByName(D,notb,'no TB outcomes')
  
  ## ===========  other counters
  ## check
  D$Set(check=1)
  D$Set(check=0,filterFun=function(x) length(x$children)>0)
  
  ## deaths
  D$Set(deaths=0)
  D$Set(deaths=1,filterFun=function(x) (x$name=='dies'))
  
  ## lives
  D$Set(lives=0)
  D$Set(lives=1,filterFun=function(x) (x$name=='survives'))
  
  ## dx tb
  D$Set(tbdx=0)
  D$Set(tbdx=1,filterFun=function(x) x$name=='DS-TB outcomes')
  
  # ## dx bac
  # D$Set(dxb=0)
  # D$Set(dxb=1,
  #       filterFun=function(x)x$name=='TB diagnosed (bacteriological)')
  ## ATT
  D$Set(att=0)
  D$Set(att=1,
        filterFun=function(x)x$name=='DS-TB tx')
        
  ## tpt tx
  D$Set(tpt=0)
  D$Set(tpt=1,filterFun=function(x) x$name=='TPT outcomes')
        
  ## TPT
  D$Set(tpttx=0)
  D$Set(tpttx=1,
        filterFun=function(x)x$name=='TPT tx')
        
  ## dx tb  tpt
  D$Set(tbdxtpt=0)
  D$Set(tbdxtpt=1,
        filterFun=function(x)x$name=='TB <1yr')        
  
  # ## referrals
  # D$Set(refers=0)
  # D$Set(refers=1,filterFun=function(x) grepl('Refer',x$name))
  
  return(D)
}

# main tree structures
contact <- txt2tree(here('indata/Child.hh.contact.txt'))
asymp <- txt2tree(here('indata/Asymptomatic.child.contact.txt'))
symp <- txt2tree(here('indata/Symptomatic.child.contact.txt'))

## merge in extras, make model of care branches, write out
SOC <- AddOutcomes(contact)

INT <- Clone(SOC)
SOC$name <- 'Facility-based model'
INT$name <- 'Community-based model'

tree2file(SOC,filename = here('indata/CSV/SOC.csv'),
          'p','cost','deaths','lives','tbdx','att','tpt','tpttx','check')

tree2file(INT,filename = here('indata/CSV/INT.csv'),
          'p','cost','deaths','lives','tbdx','att','tpt','tpttx','check')
