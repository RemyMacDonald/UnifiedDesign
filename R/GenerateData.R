GenerateData<-function(ttl=ttl,nSim=nSim,dose0=dose, n1=n,prob=prob){
  datS=list()
  iSim=1;
  for (iSim in 1:nSim) {
    i=1; tem=NULL
    for (i in 1:length(dose0)) {
      tem=cbind(tem, rbinom(n1, size=1, prob=prob[i]))
    }
    colnames(tem)=paste0("D",1:ncol(tem))
    dat1=data.frame(tem)
    # apply(dat1,2,mean)
    datS[[iSim]]=dat1
  }
  return(datS=datS)
}




