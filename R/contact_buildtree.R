## looking at tree 7
rm(list=ls())
library(here)
library(HEdtree)
library(data.tree)
library(DiagrammeR)
library(data.table)
library(glue)

# set_here('C:/Users/cm1nm/Dropbox/CONTACT')

## === outcomes subtree ===
tb <- txt2tree(here('indata/TB.outcomes.txt')) # tb dx
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
tb$`no TB tx`$p <- 'p.ptltfu'
tb$`no TB tx`$dies$p <- 'p.cfr.notx'
tb$`no TB tx`$survives$p <- '1-p.cfr.notx'
tb$`TB tx`$p <- '1-p.ptltfu'
tb$`TB tx`$dies$p <- 'p.cfr.tx'
tb$`TB tx`$survives$p <- '1-p.cfr.tx'

## -- tpt outcomes:
tpt$`Eligible for TPT`$TPT$`TB disease <1yr`$p <- 'p.tbdx.1yr*tptRR'
tpt$`Eligible for TPT`$`no TPT`$`TB disease <1yr`$p <- 'p.tbdx.1yr'
tpt$`Eligible for TPT`$TPT$`no TB disease <1yr`$p <- '1-p.tbdx.1yr*tptRR'
tpt$`Eligible for TPT`$`no TPT`$`no TB disease <1yr`$p <- '1-p.tbdx.1yr'
tpt$`Not eligible for TPT`$`TB disease <1yr`$p <- 'p.tbdx.1yr'
tpt$`Not eligible for TPT`$`no TB disease <1yr`$p <- '1-p.tbdx.1yr'
print(tpt,'p')

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
  MergeByName(D,asymp,'Asymptomatic child contact',leavesonly = TRUE)
  MergeByName(D,symp,'Symptomatic child contact')
  MergeByName(D,tpt,'TPT outcomes')
  # MergeByName(D,tpt,'no TPT outcomes')
  MergeByName(D,tb,'TB outcomes')
  MergeByName(D,notb,'no TB outcomes')
  # MergeByName(D,notb,'Child HH contact not screened')
  # MergeByName(D,notb,'No clinical re(evaluation) with/without  CXR')
  
  ## ===========  other counters
  ## check
  D$Set(check=1)
  D$Set(check=0,filterFun=function(x) length(x$children)>0)
  
  ## deaths
  D$Set(deaths=0)
  D$Set(deaths=1,filterFun=function(x) (x$name=='dies'))
  
  ## deaths for incident TB
  D$Set(incdeaths=0)
  D$`Child contact screened`$`Asymptomatic child contact`$`TPT outcomes`$Set(incdeaths=1,filterFun=function(x)x$name=='dies')
  D$`Child contact screened`$`Asymptomatic child contact`$`TPT outcomes`$Set(incdeaths=1,filterFun=function(x)x$name=='dies')
  D$`Child contact screened`$`Asymptomatic child contact`$`TPT outcomes`$Set(incdeaths=1,filterFun=function(x)x$name=='dies')
  
  ## lives
  D$Set(lives=0)
  D$Set(lives=1,filterFun=function(x) (x$name=='survives'))
  
  ## inctb
  D$Set(inctb=0)
  D$Set(inctb=1,filterFun=function(x)x$name=='TB disease <1yr')
  
  ## coprevalent TB
  
  ## dx clinical
  D$Set(tbdxc=0)
  D$Set(tbdxc=1,filterFun=function(x) x$name=='TB diagnosed (clinical)')
  
  ## dx bac
  D$Set(tbdxb=0)
  D$Set(tbdxb=1,filterFun=function(x) x$name=='TB diagnosed (bacteriological)')
  
  ## coprevalent TB total
  D$Set(prevtb=0)
  D$Set(prevtb=1,filterFun=function(x)x$name %in% c('TB diagnosed (clinical)', 'TB diagnosed (bacteriological)'))
  
  ## ATT in incidence
  D$Set(atti=0)
  D$`Child contact screened`$`Asymptomatic child contact`$`TPT outcomes`$Set(atti=1,filterFun=function(x)  x$name=='TB tx')
  
  ## ATT courses
  D$Set(att=0)
  D$Set(att=1,
        filterFun=function(x)x$name=='TB tx')
        
  ## tpt eligible
  D$Set(tpte=0)
  D$Set(tpte=1,filterFun=function(x) x$name=='Eligible for TPT')
        
  ## TPT courses
  D$Set(tpt=0)
  D$Set(tpt=1,
        filterFun=function(x)x$name=='TPT')
  
  return(D)
}

# main tree structures
contact <- txt2tree(here('indata/Child.hh.contact.txt'))
asymp <- txt2tree(here('indata/Asymptomatic.child.contact.txt'))
symp <- txt2tree(here('indata/Symptomatic.child.contact.txt'))

## merge in extras, make model of care branches, write out
tempTree <- AddOutcomes(contact)

## === SOC
SOC <- Node$new('Facility-based model')
SOC$AddChildNode(tempTree)

SOC$name <- 'Facility-based model'

tree2file(SOC,filename = here('indata/CSV/SOC.csv'),
          'p','cost','tpte', 'tpt', 'inctb','tbdxc','tbdxb', 'prevtb', 'att','lives','incdeaths','deaths','check',
          'TPT.screened','TPT.asymptomatic','TPT.eligible','TPT.treated',
          'TB.symptoms','TB.evaluated','TB.diagnosed','TB.treated')

# ## create version with probs/costs
fn <- here('indata/CSV/SOC1.csv')
if(file.exists(fn)){
        ## read
        labz <- fread(fn)
        SOC$Set(p=labz$p)
        SOC$Set(cost=labz$cost)
        ## TPT/TB cascade counters
        SOC$Set(TPT.screened=labz$TPT.screened)
        SOC$Set(TPT.asymptomatic=labz$TPT.asymptomatic)
        SOC$Set(TPT.eligible=labz$TPT.eligible)
        SOC$Set(TPT.treated=labz$TPT.treated)
        SOC$Set(TB.symptoms=labz$TB.symptoms)
        SOC$Set(TB.evaluated=labz$TB.evaluated)
        SOC$Set(TB.diagnosed=labz$TB.diagnosed)
        SOC$Set(TB.treated=labz$TB.treated)
        ## save out
        tree2file(SOC,filename = here('indata/CSV/SOC2.csv'),
                  'p','cost','deaths','lives','refers','dxc','dxb','att',
                  'check',
                  'TPT.screened','TPT.asymptomatic','TPT.eligible','TPT.treated',
                  'TB.symptoms','TB.evaluated','TB.diagnosed','TB.treated')
}

## === INT
INT <- Clone(SOC)
INT$name <- 'Community-based model'

tree2file(INT,filename = here('indata/CSV/INT.csv'),
          'p','cost','tpte', 'tpt', 'inctb','tbdxc','tbdxb', 'prevtb', 'att','lives','incdeaths','deaths','check',
          'TPT.screened','TPT.asymptomatic','TPT.eligible','TPT.treated',
          'TB.symptoms','TB.evaluated','TB.diagnosed','TB.treated')

## create version with probs/costs
fn <- here('indata/CSV/INT1.csv')
if(file.exists(fn)){
        ## read
        labz <- fread(fn)
        INT$Set(p=labz$p)
        INT$Set(cost=labz$cost)
        ## TPT/TB cascade counters
        INT$Set(TPT.screened=labz$TPT.screened)
        INT$Set(TPT.asymptomatic=labz$TPT.asymptomatic)
        INT$Set(TPT.eligible=labz$TPT.eligible)
        INT$Set(TPT.treated=labz$TPT.treated)
        INT$Set(TB.symptoms=labz$TB.symptoms)
        INT$Set(TB.evaluated=labz$TB.evaluated)
        INT$Set(TB.diagnosed=labz$TB.diagnosed)
        INT$Set(TB.treated=labz$TB.treated)
        ## save out
        tree2file(INT,filename = here('indata/CSV/INT2.csv'),
                  'p','cost','deaths','lives','refers','dxc','dxb','att',
                  'check',
                  'TPT.screened','TPT.asymptomatic','TPT.eligible','TPT.treated',
                  'TB.symptoms','TB.evaluated','TB.diagnosed','TB.treated')
}

## make functions
fnmz <- c('check','cost','deaths','att',
          'lives','refers','dxc','dxb',
          'DH.presented','DH.screened','DH.presumed','DH.treated',
          'PHC.presented','PHC.screened','PHC.presumed','PHC.treated')

SOC.F <- makeTfuns(SOC,fnmz)
INT.F <- makeTfuns(INT,fnmz)



## running all function
runallfuns <- function(D,arm='all'){
        done <- FALSE
        if('SOC' %in% arm | arm[1]=='all'){
                cat('Running functions for SOC:\n')
                for(nm in names(SOC.F)){
                        snm <- gsub('fun','',nm)
                        snma <- paste0(snm,'.soc')
                        D[[snma]] <- SOC.F[[nm]](D)
                        cat('...',snm,' run...\n')
                        done <- TRUE
                }
        }
        if('INT' %in% arm | arm[1]=='all'){
                cat('Running functions for INT:\n')
                for(nm in names(INT.F)){
                        snm <- gsub('fun','',nm)
                        snma <- paste0(snm,'.int')
                        D[[snma]] <- INT.F[[nm]](D)
                        cat('...',snm,' run...\n')
                        done <- TRUE
                }
        }


        if(!done)stop('Functions not run! Likely unrecognised arm supplied.')
        return(D)
}



## --- CHECKS
showAllParmz <- function(TREE){
        B <- showParmz(TREE)
        ## get calx
        cx <- B$calcs
        cx <- gsub("\\*|\\+|-|\\/|\\(|\\)"," ",cx)
        cx <- paste(cx,collapse=" ")
        cx <- strsplit(cx,split=" ")[[1]]
        cx <- cx[cx!=""]
        cx <- cx[cx!="1"]
        ## get non calcs
        cx <- c(cx,B$vars)
        unique(cx)
}

makeTestData <- function(ncheck,vnames){
        A <- data.table(vnames,value=runif(length(vnames)))
        A <- A[rep(1:length(vnames),each=ncheck)]
        idz <- rep(1:ncheck,length(vnames))
        A[,id:=idz]
        A[,value:=runif(nrow(A))]
        dcast(A,id~vnames,value.var = 'value')
}


## checking
vrz.soc <- showAllParmz(SOC)
vrz.spd <- showAllParmz(IPD)
vrz.idh <- showAllParmz(IDH)
vrz.iph <- showAllParmz(IPH)
vrz <- c(vrz.soc,
         vrz.spd,
         vrz.idh,
         vrz.iph
)
vrz <- unique(vrz)
A <- makeTestData(50,vrz)


## checks
IPD.F$checkfun(A) #NOTE OK
IPH.F$checkfun(A) #NOTE OK
IDH.F$checkfun(A) #NOTE OK
SOC.F$checkfun(A) #NOTE OK

## full graph out
export_graph(ToDiagrammeRGraph(contact),
             file_name=here('plots/contact.pdf'))
