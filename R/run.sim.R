# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' @title A Unified Platform for Comprehensive Evaluation of Phase I Clinical Trial Designs
#'
#' @description This is a description
#'
#' @param prob Vector of true DLT rates
#' @param nSim Number of studies to simulate
#' @param dose Dose levels of recruited patients
#' @param ttl Target toxicity level
#' @param n The maximum sample size of the trial
#' @param m The cohort size
#' @param Max_FoldChange Maximum fold change
#' @param epsilon1
#' @param epsilon2
#' @param pod1 Eliminate the current dose and higher if the posterior probability of exceeding the target toxicity level exceeds this value
#' @param alpha Prior shape1 parameter
#' @param beta Prior shape2 parameter
#' @param cohort.size Vector of cohort sizes
#' @param theta Safety threshold for cohort-sequence design
#' @param alpha1 Prior shape1 parameter for cohort-sequence design
#' @param beta1 Prior shape2 parameter for cohort-sequence design
#'
#'
#' @return \code{run.sim( )} returns a list of simulation results
#' @export
#'
#' @examples
#' run.sim(prob=c(0.21, 0.30, 0.37, 0.43, 0.48, 0.52, 0.56, 0.61, 0.64, 0.66))
#'
run.sim <- function(prob,nSim=10,dose=3^{0:(length(prob)-1)},ttl=0.3,n=60,m=3,Max_FoldChange=3,epsilon1 = 0.05,epsilon2 = 0.05,pod1=0.95,
                  alpha=1,beta=1,cohort.size=c(3,5,8,10),theta=0.25,alpha1=1,beta1=4) {
  if(length(prob)!=length(dose)){
    stop('prob and dose should have equal length')
  }
  if(!all(prob == cummax(prob))){
    stop('prob should be an increasing sequence')
  }
  if(!all(dose == cummax(dose))){
    stop('dose should be an increasing sequence')
  }
  if(!all(prob>=0 & prob<=1)){
    stop('Each value in prob should be between 0 and 1')
  }
  if(!all(dose>0)){
    stop('Each value in dose should be greater than 0')
  }
  datS <- GenerateData(ttl=ttl,nSim=nSim,dose0=dose,n1=n,prob=prob)
  bd=get.boundary(target=ttl, ncohort=5, cohortsize=3)#, p.saf = 0.6 * ttl, p.tox = 1.4 * ttl)
  bd1=c(bd$lambda_e, bd$lambda_d)
  out1=out2=out3=out4=out5=out6=out7=NULL
  im=1
  jobid=1
  iSim=1
  for (iSim in 1:nSim){
    print(iSim)
    dat1=datS[[iSim]]
    fiti1=RFun_BOIN(dati=dat1,n,m,dose0=dose,ttl=ttl,bd=bd,pod1=pod1,MinN_atMTD=6,Max_FoldChange=Max_FoldChange,iTer=TRUE)
    fiti2=RFun_33(dati=dat1,n,m,dose0=dose)
    fiti3=RFun_i33(dati=dat1,n,m,dose0=dose,ttl=ttl,epsilon1=epsilon1,epsilon2=epsilon2,pod1=pod1,MinN_atMTD=6,Max_FoldChange=Max_FoldChange)
    fiti4=RFun_mTPI(dati=dat1,n,m,dose0=dose,ttl=ttl,alpha=alpha,beta=beta,epsilon1=epsilon1,epsilon2=epsilon2,pod1=pod1,MinN_atMTD=6,Max_FoldChange=Max_FoldChange)
    fiti5=RFun_mTPI2(dati=dat1,n,m,dose0=dose,ttl=ttl,alpha=alpha,beta=beta,epsilon1=epsilon1,epsilon2=epsilon2,pod1=pod1,MinN_atMTD=6,Max_FoldChange=Max_FoldChange)
    fiti6=RFun_CS(dati=dat1,n,dose0=dose,cohort.size,theta,ttl=ttl,epsilon1=epsilon1,epsilon2=epsilon2,pod1=pod1,MinN_atMTD=6,Max_FoldChange=3,alpha1=1,beta1=4)


    ###
    fitSi=list(); nameSi=NULL; isi=1
    fitSi[[isi]]=fiti1; nameSi=c(nameSi, "BOIN"); isi=isi+1
    fitSi[[isi]]=fiti2; nameSi=c(nameSi, "3+3"); isi=isi+1
    fitSi[[isi]]=fiti3; nameSi=c(nameSi, "i3+3"); isi=isi+1
    fitSi[[isi]]=fiti4; nameSi=c(nameSi, "mTPI"); isi=isi+1
    fitSi[[isi]]=fiti5; nameSi=c(nameSi, "mTPI2"); isi=isi+1
    fitSi[[isi]]=fiti6; nameSi=c(nameSi, "CS"); isi=isi+1

    ###
    outi=NULL
    for (i in 1:length(fitSi)) {
      fitij=fitSi[[i]]
      outi=rbind(outi, data.frame(Scenario=im, jobid=jobid, simid=iSim, method=nameSi[i], N=sum(fitij$outN[2,]), DLT=sum(fitij$outN[1,]), MTD=dose[tail(fitij$path$D1,1)], stop_reason=fitij$stop_reason))
    }
    out1=rbind(out1, outi)

    outi=NULL
    for (i in 1:length(fitSi)) {
      fitij=fitSi[[i]]
      outi=rbind(outi, data.frame(Scenario=im, jobid=jobid, simid=iSim, method=nameSi[i], fitij$outData))
    }
    out2=rbind(out2, outi)
  }
  out=list(out1,out2)
  return(out)
}
