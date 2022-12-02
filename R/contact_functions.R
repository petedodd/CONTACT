## TODO questions:
## 
## - stool/sputum
## - flag assumption = groups for SA or pending data
## 

## ========= UTILITIES ===============
logit <- function(x) log(odds(x))
ilogit <- function(x) iodds(exp(x))
odds <- function(x) x/(1-x)
iodds <- function(x) x/(1+x)
lo <- function(x) quantile(x,probs = 0.025, na.rm=TRUE)
hi <- function(x) quantile(x,probs = 1-0.025, na.rm=TRUE)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

brkt <- function(M,L,H,ndp=0) paste0(round(M,ndp),' (',
                                     round(L,ndp),' - ',
                                     round(H,ndp),')')
gm <- function(x) exp(mean(log(x))) #geometric mean
gh <- function(x) glue(here(x))

## ========= OUTCOMES ===============
## TODO - remove excess RNG here
CFRtxY <- function(age,hiv=0,art=0,P){#NB optimized for clarity not speed
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$ontx.u5$r(length(age))
  tmp[age=='5-14'] <- P$ontx.o5$r(sum(age=='5-14'))  #NB this could be achieved in the tree model
  ## hivartOR
  Z <- P$hivartOR$r(length(age))
  hor <- rep(1,length(age))
  tmp <- logit(tmp)                     #transform
  tmp[hiv>0] <- tmp[hiv>0]+Z[hiv>0,1]
  tmp[art>0] <- tmp[art>0]+Z[art>0,2]
  tmp <- ilogit(tmp)                    #inverse transform
  tmp
}
## CFRtxY(1:10,P)                            #test
## summary(CFRtxY(1:1e3,P))
## summary(CFRtxY(1:1e3,hiv=1))
## summary(CFRtxY(1:1e3,hiv=1,art=1,P))


## == CFR off tx
CFRtxN <- function(age,hiv=0,art=0,P){
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$notx.u5$r(length(age))          #default a<5 and hiv=art=0
  tmp[age!='5-14' & hiv>0 & art==0] <- P$notxH.u5$r(sum(age!='5-14' & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[age!='5-14' & hiv>0 & art>0] <- P$notxHA.u5$r(sum(age!='5-14' & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[age=='5-14'] <- P$notx.o5$r(sum(age=='5-14'))    #o5, HIV-ve
  tmp[age=='5-14' & hiv>0 & art==0] <- P$notxH.o5$r(sum(age=='5-14' & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[age=='5-14' & hiv>0 & art>0] <- P$notxHA.o5$r(sum(age=='5-14' & hiv>0 & art>0)) #o5,HIV+,ART+
  tmp
}
## CFRtxN(1:10,P)                            #test
## summary(CFRtxN(1:1e3,P))
## summary(CFRtxN(1:1e3,hiv=1,P))
## summary(CFRtxN(1:1e3,hiv=1,art=1,P))




## add CFRs to data by side-effect
AddCFRs <- function(D,P){
  ## d.cfr.notx & d.cfr.tx
  D[,c('cfr.notx','cfr.tx'):=0] #NOTE neglect non-TB mortality
  ## CFR on  ATT
  D[,cfr.tx:=CFRtxY(age,hiv,art,P)]
  ## CFR w/o ATT
  D[,cfr.notx:=CFRtxN(age,hiv,art,P)]
}

## == LTBI infection probability
#NB this is LTBI given not active: it is taken to be max(0,LTBI-coprev)
ltbi.prev <- function(age,coprev,hinco=FALSE){
        if(length(age)>1 & length(hinco)==1) hinco <- rep(hinco,length(age))
        tmp <- P$LTBI04$r(length(age))
        tmp[hinco] <- P$LTBI04hi$r(sum(hinco))
        tmp[age>=5] <- P$LTBI514$r(sum(age>=5))
        tmp[age>=5 & hinco] <- P$LTBI514hi$r(sum(age>=5 & hinco))
        tmp
        ## pmax(0,tmp - coprev) # already taken into account with decision tree
}
# ltbi.prev(1:10,0.1,hinco=TRUE)

## progression probability
progprob <- function(age,hiv=0,art=0, P){
        if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
        if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
        ans <- P$prog04$r(length(age))
        ans[age>=5] <- P$prog514$r(sum(age>=5))
        if(any(hiv>0)){ #treat as IRR for escape
                hr <- P$hivpi$r(sum(hiv>0))
                ans[hiv>0] <- 1-(1-ans[hiv>0])^hr
        }
        if(any(art>0)){ #treat as IRR for escape
                hr <- P$artp$r(sum(art>0))
                ans[art>0] <- 1-(1-ans[art>0])^hr
        }
        ans
}

## add CFRs to data by side-effect
AddProgProb <- function(D, P){
  D[,p.tbdx.1yr:=progprob(age,hiv,art,P) * ltbi.prev(age,0)] #NOTE handling of coprev happens explicitly
}

## progprob(c(rep(3,5),rep(10,5)), P=P)


## === IPT efficacy
TPTrr <- function(age,hiv=0,
                  tst='none'            #a flag for: given to TST+ or not
){
        if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
        if(tst=='none'){
                ans <- P$tptRR$r(length(age))
                ans[hiv>0] <- P$tptRRhivpos$r(sum(hiv>0))  #HIV+
        } else {
                ans <- P$tptRRhivpos$r(length(age))  #NOTE only applied to HIV-ves
        }
        ans
}
# TPTrr(1:10, P)
# summary(TPTrr(runif(1e3),hiv=0)) #0.37
# summary(TPTrr(runif(1e3),hiv=1)) #0.35
# summary(TPTrr(runif(1e3),tst='yes')) #0.09
AddTPTrr <- function(D,P){
        D[,tptRR0:=TPTrr(age,hiv)]  #base efficacy of IPT
        D[,tptRR1:=TPTrr(age,tst="+ve")]   #base efficacy of PT: TST+ve
        D[,tptRR:=tptRR0] 
}

## ======= COMBINED LABELLER ===========

## additional labels from data (may overwrite some initial version currently)
AddDataDrivenLabels <- function(D){
        
        # 
        # D[,int.enrolled_per_index_case:=ifelse(age=='0-4', int.enrolled_per_index_case*int.frac.enrolled.u5, int.enrolled_per_index_case*(1-int.frac.enrolled.u5))]
        # D[,soc.enrolled_per_index_case:=ifelse(age=='0-4', soc.enrolled_per_index_case*soc.frac.enrolled.u5, soc.enrolled_per_index_case*(1-soc.frac.enrolled.u5))]
        #
        D[,int.frac.screened:=soc.frac.screened]
        
        D[,soc.frac.symp.ltfu:=(1-soc.frac.symp.attending)]
        D[,int.frac.symp.ltfu:=(1-int.frac.symp.attending)]
        
        # TPT eligibility and initiation for 5-14 years 
        D[,soc.frac.tpt.eligible:=ifelse(age=='0-4' | age=='5-14' & hiv==1, soc.frac.tpt.eligible, 0)]
        D[,int.frac.tpt.eligible:=ifelse(age=='0-4' | age=='5-14' & hiv==1, int.frac.tpt.eligible, 0)]
        
        
        # fraction under 5
        # D[,F.u5:=ifelse(isoz=='CMR' & arms=='SOC', cmr.soc.frac.u5, 0)]
        # D[,F.u5:=ifelse(isoz=='CMR' & arms=='INT', cmr.int.frac.u5, F.u5)]
        # D[,F.u5:=ifelse(isoz=='UGA' & arms=='SOC', uga.soc.frac.u5, F.u5)]
        # D[,F.u5:=ifelse(isoz=='UGA' & arms=='INT', uga.int.frac.u5, F.u5)]
        
        # 
        # D[,soc.frac.tpt.initiated:=ifelse(age=='5-14' & hiv==1, soc.frac.tpt.initiated, 0)]
        # D[,int.frac.tpt.initiated:=ifelse(age=='5-14' & hiv==1, int.frac.tpt.initiated, 0)]
        
        # frac u5, HIV prevalence and ART coverage
        # D[,hivprev.u5:=ifelse(isoz=='CMR' & arms=='SOC', cmr.soc.frac.hiv, 0)]
        # D[,hivprev.u5:=ifelse(isoz=='CMR' & arms=='INT', cmr.int.frac.hiv, hivprev.u5)]
        # D[,hivprev.u5:=ifelse(isoz=='UGA' & arms=='SOC', uga.soc.frac.hiv, hivprev.u5)]
        # D[,hivprev.u5:=ifelse(isoz=='UGA' & arms=='INT', uga.int.frac.hiv, hivprev.u5)]
        # D[,hivprev.o5:=hivprev.u5]
}
## combined function to add the labels to the tree prior to calculations
MakeTreeParms <- function(D,P){
  ## -- use of other functions
  # AddSampleTests(D) #samples/tests
  AddCFRs(D,P) #outcomes
  AddTPTrr(D,P)
  AddProgProb(D, P)
  ## new labels from data
  AddDataDrivenLabels(D)
}

## ======= EPIDEMIOLOGY ===========

makeAttributes <- function(D){
    nrep <- nrow(D)
    D[,id:=1:nrep]
    fx <- list(age=agelevels,hiv=hivlevels,art=artlevels,isoz=isoz,arms=c('SOC','INT'))
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
    ## --- age
    D[,F.u5:=ifelse(isoz=='CMR' & arms=='SOC', cmr.soc.frac.u5, 0)]
    D[,F.u5:=ifelse(isoz=='CMR' & arms=='INT', cmr.int.frac.u5, F.u5)]
    D[,F.u5:=ifelse(isoz=='UGA' & arms=='SOC', uga.soc.frac.u5, F.u5)]
    D[,F.u5:=ifelse(isoz=='UGA' & arms=='INT', uga.int.frac.u5, F.u5)]

    D[,hivprev.u5:=ifelse(isoz=='CMR' & arms=='SOC', cmr.soc.frac.hiv, 0)]
    D[,hivprev.u5:=ifelse(isoz=='CMR' & arms=='INT', cmr.int.frac.hiv, hivprev.u5)]
    D[,hivprev.u5:=ifelse(isoz=='UGA' & arms=='SOC', uga.soc.frac.hiv, hivprev.u5)]
    D[,hivprev.u5:=ifelse(isoz=='UGA' & arms=='INT', uga.int.frac.hiv, hivprev.u5)]
    D[,hivprev.o5:=hivprev.u5]
              
    D[,value:=ifelse(age=='5-14',1-F.u5,F.u5),
      by=.(id,isoz,arms)] #NOTE value first set here
    ## --- HIV/ART
    D[,h01:=0]
    D[age!='5-14',h10:=hivprev.u5*(1-artcov)]
    D[age=='5-14',h10:=hivprev.o5*(1-artcov)]
    D[age!='5-14',h00:=1-hivprev.u5]
    D[age=='5-14',h00:=1-hivprev.o5]
    D[age!='5-14',h11:=hivprev.u5*artcov]
    D[age=='5-14',h11:=hivprev.o5*artcov]
    D[hiv==0 & art==0,value:=value*h00]
    D[hiv==0 & art==1,value:=value*h01]
    D[hiv==1 & art==0,value:=value*h10]
    D[hiv==1 & art==1,value:=value*h11]
    D[,c('h00','h01','h10','h11'):=NULL]
    ## --- TB
    ## ## (old version) calculate true TB prev among initial care seeking as:
    ## ## tbi = f x tbp + (1-f) x tbd
    ## ## where: f=fraction initially seeking care at PHC; tbp=prev @ phc; tbd=prev @ dh
    ## ## NOTE the 'underlying' TB prev in care-seekers in controlled by soc parms
    ## D[,tbi:= d.soc.pphc * d.phc.tbinprsmptv + (1-d.soc.pphc) * d.dh.tbinprsmptv]
    ## D[tb=='noTB',value:=value*(1-tbi)]
    ## D[tb=='TB-',value:=value*tbi*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    ## D[tb=='TB+',value:=value*tbi*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
    ## D[,tbi:=NULL]                            #remove temporary variable
    # ## (new version) based on data:
    # D[,tbi:=ifelse(age=='5-14',d.TBprev.ICS.o5,d.TBprev.ICS.u5)]
    # D[tb=='noTB',value:=value*(1-tbi)]
    # D[tb=='TB-',value:=value*tbi*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    # D[tb=='TB+',value:=value*tbi*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
    # D[,tbi:=NULL]                            #remove temporary variable
    return(D)
}



## function for generating random sample of costs
MakeCostData <- function(csts,          #base data table of cost data
                         nrep,          #number of replicates being used in PSA
                         anmz=NULL     #attribute names (if any)
                         ){
  if(nrow(csts[cost.sd>0 & cost.m==0])>0) warning(paste0('Some cost input variables have zero mean & SD>0. These will be treated as fixed variables:\n',paste0(csts[cost.sd>0 & cost.m==0,cost],collapse='\n')))
  if(is.null(anmz)& any(csts[,table(cost)]>1)) warning('Some cost names occur >1 times, but no attributes have been specified! This is unlikely to do what you want.')
  csts[cost.m>0,gmsc:=cost.sd^2/cost.m]
  csts[!is.na(gmsc) & gmsc > 0, gmk:=cost.m/gmsc]
  NR <- nrow(csts)
  csts <- csts[rep(1:NR,nrep)]
  csts[,id:=rep(1:nrep,each=NR)]
  csts[,rnd:=!is.na(gmsc) & !is.na(gmk) & gmk>0 & gmsc > 0]
  csts[rnd==TRUE,value:=rgamma(sum(rnd),shape=gmk,scale = gmsc)] #random sample from gamma distribution
  csts[rnd!=TRUE,value:=cost.m]                                  #fixed values
  ## csts[,cnms:=paste0('c_',cost)]
  csts[,cnms:=paste0(cost)]
  F <- 'id '
  if(!is.null(anmz)) F <- paste0(F,'+ ',paste(anmz,collapse='+')) #split out by attributes if included
  F <- paste0(F, ' ~ cnms')
  dcast(csts,as.formula(F),value.var = 'value')      #id ~ cnms
}
## NOTE
## if attributes are included, all costs need to be specified by them even if this means duplicating those without dependence


## making life years
GetLifeYears <- function(isolist,discount.rate,yearfrom){
    ## template:
    LYT <- data.table(age=0:14,
                      age_group=c(rep('0-4',5),rep('5-14',10)),
                      LYS=0.0)
    ## make country/age key
    LYK <- list()
    for(iso in isolist){
        ## iso <- cn
        tmp <- copy(LYT)
        tmp[,iso3:=iso]
        for(ag in tmp$age)
            tmp[age==ag,LYS:=discly::discly(iso3=iso,
                                            age=ag,
                                            yearnow=yearfrom,
                                            sex='Total',
                                            endyear = 2098,
                                            HR=1,
                                            dr=discount.rate,
                                            hiv='both'
                                            )]
        LYK[[iso]] <- tmp
    }
    LYK <- rbindlist(LYK)
    ## assume unweighted & collapse
    LYK <- LYK[,.(LYS=mean(LYS)),by=.(iso3,age=age_group)]
    setkey(LYK,age)
    LYK
}

## ## scraps for development of below fn
## data <- merge(D,LYK,by='age') #add age
## Kmax <- 1e3
## file.id <- 'test'
## wtp <- 500

## NOTE this is more illustrative for now
## NOTE needs a folder called graphs/ creating (which is currently excluded from the repo)
## some automatic CEA outputs
MakeCEAoutputs <- function(data,LY,
                           file.id='',Kmax=5e3,wtp=5e3,
                           arms=c('SOC','INT')){
  data <- merge(data,LY,by='age') #add age
  DS <- data[,.(cost.SOC=sum(cost.soc*value),
                cost.INT=sum(cost.int*value),
                lyl.SOC=sum(deaths.soc*value*LYS),
                lyl.INT=sum(deaths.int*value*LYS)),
             by=id] #PSA summary

  ## prep for BCEA
  LYS <- CST <- matrix(nrow=nreps,ncol=2)
  LYS[,1] <- 1-DS$lyl.SOC #NOTE this is life years lost
  LYS[,2] <- 1-DS$lyl.INT
  CST[,1] <- DS$cost.SOC
  CST[,2] <- DS$cost.INT
  ## BCEA outputs
  M <- bcea(e=LYS,c=CST,ref=1,interventions = arms,Kmax=Kmax)
  print(summary(M))

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

  ceplane.plot(M,graph='ggplot2',wtp=wtp)+
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

## --- for reformatting costs
reformatCosts <- function(rcsts){
  iextra <- outer(isoz,c('.lo','.hi','drop'),paste0)
  iextra <- c(t(iextra)); iextra <- rev(rev(iextra)[-1])
  nnmz <- c('drop','DESCRIPTION','NAME',iextra)
  names(rcsts)[1:length(nnmz)] <- nnmz
  drop <- grep('drop',names(rcsts),value=TRUE)
  rcsts[,c(drop):=NULL]
  rcsts[is.na(rcsts)] <- 0.1 #dummy
  rcsts <- melt(rcsts,id=c('NAME','DESCRIPTION'))
  rcsts[,DESCRIPTION:=NULL]
  rcsts[,c('iso3','hilo'):=tstrsplit(variable,split="\\.")]
  rcsts <- dcast(rcsts,iso3 + NAME ~ hilo,value.var = 'value')
  rcsts[,c('cost.m','cost.sd'):=.((lo+hi)/2,(hi-lo)/3.92)]
  rcsts <- rcsts[,.(iso3,cost=NAME,cost.m,cost.sd)]
  rcsts
}

## calculate mid/lo/hi
MLH <- function(dat){
  nnmz <- names(dat)
  lnmz <- paste0(nnmz,'.lo')
  hnmz <- paste0(nnmz,'.hi')
  mnmz <- paste0(nnmz,'.mid')
  L <- dat[,lapply(.SD,lo),.SDcols=nnmz]
  M <- dat[,lapply(.SD,mean),.SDcols=nnmz]
  H <- dat[,lapply(.SD,hi),.SDcols=nnmz]
  setnames(L,nnmz,lnmz); setnames(M,nnmz,mnmz); setnames(H,nnmz,hnmz);
  list(L=L,M=M,H=H)
}
## MLH(out[,.(DcostperATT.int,DcostperATT.soc)]) #test


## =========== output formatters
outsummary <- function(out){

  ## mid/lo/hi
  outa <- MLH(out[,.(costperATT.soc,costperATT.int,
                     DcostperATT.int,
                     costperTPT.soc,costperTPT.int,
                     DcostperTPT.int,
                     Ddeaths.int,
                     DLYL.int,
                     DLYL0.int,
                     DcostperLYS0.int,
                     DcostperLYS.int,
                     Dcostperdeaths.int,
                     Dcost.int)])

  ## more bespoke statistics
  outi <- out[,.(ICER.int= -mean(Dcost.int) / mean(DLYL.int))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine

  ## pretty version
  pouts <- outs[,.(costperATT.soc = brkt(costperATT.soc.mid,costperATT.soc.lo,costperATT.soc.hi),
                   costperATT.int = brkt(costperATT.int.mid,costperATT.int.lo,costperATT.int.hi),
                   DcostperATT.int = brkt(DcostperATT.int.mid,DcostperATT.int.lo,DcostperATT.int.hi), # TODO::quick gap measure check!!
                   costperTPT.soc = brkt(costperTPT.soc.mid,costperTPT.soc.lo,costperTPT.soc.hi),
                   costperTPT.int = brkt(costperTPT.int.mid,costperTPT.int.lo,costperTPT.int.hi),
                   DcostperTPT.int = brkt(DcostperTPT.int.mid,DcostperTPT.int.lo,DcostperTPT.int.hi),
                   DcostperLYS0.int = brkt(DcostperLYS0.int.mid,
                                           DcostperLYS0.int.lo,DcostperLYS0.int.hi),
                  DcostperLYS.int = brkt(DcostperLYS.int.mid,DcostperLYS.int.lo,DcostperLYS.int.hi),
                  Dcostperdeaths.int = brkt(Dcostperdeaths.int.mid,
                                            Dcostperdeaths.int.lo,Dcostperdeaths.int.hi),
                  DcostPerChildContact.int = brkt(Dcost.int.mid,Dcost.int.lo,Dcost.int.hi),
                  DdeathsPer100kChildContacts.int = brkt(-1e5*Ddeaths.int.mid,
                                               -1e5*Ddeaths.int.hi,-1e5*Ddeaths.int.lo),
                  DLYS0Per100kChildContacts.int = brkt(-1e5*DLYL0.int.mid,
                                               -1e5*DLYL0.int.hi,-1e5*DLYL0.int.lo),
                  DLYSPer100kChildContacts.int = brkt(-1e5*DLYL.int.mid,
                                            -1e5*DLYL.int.hi,-1e5*DLYL.int.lo),
                  ICER.int=round(ICER.int,0))]

  ## return value
  list(outs=outs,pouts=pouts)
}


## ---- utilities for making CEACs
make.ceac <- function(CEA,lamz){
    crv <- lamz
    for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
    crv
}
