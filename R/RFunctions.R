RFun_BOIN = function(dati, n, m, dose0, ttl=0.3, bd=c(0,1), pod1=0.95, MinN_atMTD=6, Max_FoldChange=3, iTer=TRUE) {
  # dati=dat1; n=n1; m=mi; pod1=0.95; iTer=TRUE; ttl=0.3

  # bd=get.boundary(target=ttl, ncohort=n/m, cohortsize=m)
  # # bd=get.boundary(ttl, n/m, m)
  # bd=bd$full_boundary_tab

  nd=length(dose0)

  ## Safety bound
  ni=1:n; neli=NULL
  for (i in 1:length(ni)) {
    nj=ni[i]; yj=0:nj
    ppj=pbeta(ttl, 1+yj, 1+nj-yj, lower.tail=FALSE)

    ii=which(ppj>pod1) # which DLT gives probability > 0.95
    if (length(ii)>0) {
      neli=c(neli,yj[min(ii)]) # which DLT does the minimum correspond to
    } else {
      neli=c(neli, NA)
    }
  }
  neli[1:2]=NA
  deli=max(dose0)+100
  id1=nd+1
  # which(bd$full_boundary_tab[4,]!=neli)


  id=1; Yi=NULL; dosei=NULL

  iec=1; outD=NULL; outD0=NULL; outP=NULL; Yj=NULL
  outDTL=NULL; outNd=NULL
  outN=outN0=numeric(nd)
  outM=outDLT=NULL
  fiti=NULL
  id2=NULL
  repeat {

    temi=dati[,id]; indexi=which(!is.na(temi)) # DLTs at current dose and indices of non-NA values
    dati[indexi[1:m],id]=NA # make the first m rows NA
    Yi=c(Yi,temi[indexi[1:m]]) # DLT's in cohort
    dosei=c(dosei,rep(dose0[id],m)) # The m doses taken

    outD0=c(outD0, id) # Dose taken
    outDLT=c(outDLT, sum(temi[indexi[1:m]])) # Sum DLTs
    outM=c(outM, m)

    outDTLi=numeric(nd) # Vector of n_d zeros
    outDTLi[id]=sum(temi[indexi[1:m]]) # Sum DLTs
    outDTL=rbind(outDTL, outDTLi)

    outNdi=numeric(nd) # Vector of n_d zeros
    outNdi[id]=m # Number of patients
    outNd=rbind(outNd, outNdi)

    outN[id]=outN[id]+sum(temi[indexi[1:m]])
    outN0[id]=outN0[id]+m


    if (outN[id]>=neli[outN0[id]] & !is.na(neli[outN0[id]]) & iTer) {

      deli=c(deli, dose0[id])
      id=id-1
      outD=c(outD,id)

      id1=pmin(id,id1) #highest candidate dose

      if (id<=0) {
        stop_reason="All toxic"
        break
      }


    } else {

      # indexi=which(bd[1,]==outN0[id])
      # if (outN[id]<=bd[2,indexi]) {
      #   id=pmin(pmin(id+1,nd), id1)
      # } else if (outN[id]>=bd[3,indexi]) {
      #   id=pmin(pmax(1,id-1), id1)
      # }

      if ((outN[id]/outN0[id])<=bd[1]) {
        id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
        id=pmin(pmin(id+1,nd), id2)
        id=pmin(id,id1)

      } else if ((outN[id]/outN0[id])>=bd[2]) {
        id=pmin(pmax(1,id-1), id1)
      }

      outD=c(outD,id)
    }


    if (sum(outN0)+m > n) {
      stop_reason="MaxN reached"

      while (outN0[id] < MinN_atMTD && id >= 2) {
        id=id-1
      }

      break
    }
    # if (iec>=nec) break
    # iec=iec+1
  }

  fiti=select.mtd(target=ttl, npts=outN0, ntox=outN)

  outD[length(outD)]=pmin(fiti$MTD,ifelse(id1==0,1000,id1))
  # outD[length(outD)]=fiti$MTD

  colnames(outDTL)=colnames(outNd)=paste0("P",1:nd)
  rownames(outDTL)=rownames(outNd)=1:nrow(outDTL)

  tmpi=outNd
  tmpi[tmpi>0]=paste0("/",tmpi[tmpi>0])
  outND=matrix(paste0(outDTL,tmpi),ncol=nd)
  colnames(outND)=paste0("P",1:nd)
  rownames(outND)=1:nrow(outDTL)
  outND[outND=="00"]=""; outND=as.data.frame(outND)

  outDTL[outDTL==0]=""; outDTL=as.data.frame(outDTL)
  outNd[outNd==0]=""; outNd=as.data.frame(outNd)


  outi2=rbind(outN, outN0)
  colnames(outi2)=paste0("P",1:nd)
  rownames(outi2)=c("DLT","Total")

  outi=data.frame(D0=outD0, D1=outD)


  outi3=cbind(Dose1=outD0, Dose=dose0[outD0], DLT=outDLT, M=outM)

  out=list(Y=Yi, dose=dosei, BOIN=fiti, outData=outi3, outN=outi2, path=outi, outND=outND, stop_reason=stop_reason)
  return(out)
}

RFun_mTPI = function(dati, n, m, dose0, ttl=0.3, alpha=1, beta=1, epsilon1 = 0.05, epsilon2 = 0.05, pod1=0.95, MinN_atMTD=6, Max_FoldChange=3) {

  nd=length(dose0)
  Under=c(0,ttl-epsilon1)
  Equiv=c(ttl-epsilon1,ttl+epsilon2)
  Over=c(ttl+epsilon2,1)
  ## Safety bound
  ni=1:n; neli=NULL
  for (i in 1:length(ni)) {
    nj=ni[i]; yj=0:nj
    ppj=pbeta(ttl, 1+yj, 1+nj-yj, lower.tail=FALSE)

    ii=which(ppj>pod1)
    if (length(ii)>0) {
      neli=c(neli,yj[min(ii)])
    } else {
      neli=c(neli, NA)
    }
  }
  neli[1:2]=NA
  deli=max(dose0)+100
  id1=nd+1

  id=1; Yi=NULL; dosei=NULL

  iec=1; outD=NULL; outD0=NULL; outP=NULL; Yj=NULL
  outDTL=NULL; outNd=NULL
  outN=outN0=numeric(nd)
  outM=outDLT=NULL
  fiti=NULL
  alpha.post=NULL; beta.post=NULL
  repeat {

    temi=dati[,id]; indexi=which(!is.na(temi)) # DLTs at current dose and indices of non-NA values
    dati[indexi[1:m],id]=NA # make the first m rows NA
    Yi=c(Yi,temi[indexi[1:m]]) # DLT's in cohort
    dosei=c(dosei,rep(dose0[id],m)) # The m doses taken

    outD0=c(outD0, id) # Dose taken
    outDLT=c(outDLT, sum(temi[indexi[1:m]])) # Sum DLTs
    outM=c(outM, m)

    outDTLi=numeric(nd) # Vector of n_d zeros
    outDTLi[id]=sum(temi[indexi[1:m]]) # Sum DLTs
    outDTL=rbind(outDTL, outDTLi)

    outNdi=numeric(nd) # Vector of n_d zeros
    outNdi[id]=m # Number of patients
    outNd=rbind(outNd, outNdi)

    outN[id]=outN[id]+sum(temi[indexi[1:m]])
    outN0[id]=outN0[id]+m

    alpha.post[id] <- alpha+outN[id]
    beta.post[id] <- beta+outN0[id]-outN[id]

    UPM.under <- (pbeta(Under[2],alpha.post[id],beta.post[id])-pbeta(Under[1],alpha.post[id],beta.post[id]))/(Under[2]-Under[1])
    UPM.equiv <- (pbeta(Equiv[2],alpha.post[id],beta.post[id])-pbeta(Equiv[1],alpha.post[id],beta.post[id]))/(Equiv[2]-Equiv[1])
    UPM.over <- (pbeta(Over[2],alpha.post[id],beta.post[id])-pbeta(Over[1],alpha.post[id],beta.post[id]))/(Over[2]-Over[1])
    max <- which.max(c(UPM.under,UPM.equiv,UPM.over))

    if (outN[id]>=neli[outN0[id]] & !is.na(neli[outN0[id]])) {

      deli=c(deli, dose0[id])
      id=id-1
      outD=c(outD,id)

      id1=pmin(id,id1)

      if (id<=0) {
        stop_reason="All toxic"
        break
      }


    } else {

      if (max==1) {
        id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
        id=pmin(pmin(id+1,nd), id2)
        id=pmin(id,id1)

      } else if (max==3) {
        id=pmin(pmax(1,id-1), id1)
      }

      outD=c(outD,id)
    }


    if (sum(outN0)+m > n) {
      stop_reason="MaxN reached"

      while (outN0[id] < MinN_atMTD && id >= 2) {
        id=id-1
      }

      break
    }
    # if (iec>=nec) break
    # iec=iec+1
  }

  fiti=select.mtd(target=ttl, npts=outN0, ntox=outN)

  outD[length(outD)]=pmin(fiti$MTD,ifelse(id1==0,1000,id1))
  # outD[length(outD)]=fiti$MTD

  colnames(outDTL)=colnames(outNd)=paste0("P",1:nd)
  rownames(outDTL)=rownames(outNd)=1:nrow(outDTL)

  tmpi=outNd
  tmpi[tmpi>0]=paste0("/",tmpi[tmpi>0])
  outND=matrix(paste0(outDTL,tmpi),ncol=nd)
  colnames(outND)=paste0("P",1:nd)
  rownames(outND)=1:nrow(outDTL)
  outND[outND=="00"]=""; outND=as.data.frame(outND)

  outDTL[outDTL==0]=""; outDTL=as.data.frame(outDTL)
  outNd[outNd==0]=""; outNd=as.data.frame(outNd)


  outi2=rbind(outN, outN0)
  colnames(outi2)=paste0("P",1:nd)
  rownames(outi2)=c("DLT","Total")

  outi=data.frame(D0=outD0, D1=outD)


  outi3=cbind(Dose1=outD0, Dose=dose0[outD0], DLT=outDLT, M=outM)

  out=list(Y=Yi, dose=dosei, BOIN=fiti, outData=outi3, outN=outi2, path=outi, outND=outND, stop_reason=stop_reason)
  return(out)
}

RFun_mTPI2 = function(dati, n, m, dose0, ttl=0.3, alpha=1, beta=1, epsilon1 = 0.05, epsilon2 = 0.05, pod1=0.95, MinN_atMTD=6, Max_FoldChange=3) {

  nd=length(dose0)
  Equiv=c(ttl-epsilon1,ttl+epsilon2)
  Under<-Over<-list()
  temp <- ttl-epsilon1-(epsilon1+epsilon2)
  i=1
  while(temp>0){
    Under[[i]] <- c(ttl-epsilon1-i*(epsilon1+epsilon2),ttl-epsilon1-(i-1)*(epsilon1+epsilon2))
    temp <- ttl-epsilon1-(i+1)*(epsilon1+epsilon2)
    i=i+1
  }
  if(Under[[length(Under)]][1]>0){
    Under[[length(Under)+1]] <- c(0,Under[[length(Under)]][1])
  }

  temp <- ttl+epsilon2+(epsilon1+epsilon2)
  i=1
  while(temp<1){
    Over[[i]] <- c(ttl+epsilon2+(i-1)*(epsilon1+epsilon2),ttl+epsilon2+i*(epsilon1+epsilon2))
    temp <- ttl+epsilon2+(i+1)*(epsilon1+epsilon2)
    i=i+1
  }
  if(Over[[length(Over)]][2]<1){
    Over[[length(Over)+1]] <- c(Over[[length(Over)]][2],1)
  }

  ## Safety bound
  ni=1:n; neli=NULL
  for (i in 1:length(ni)) {
    nj=ni[i]; yj=0:nj
    ppj=pbeta(ttl, 1+yj, 1+nj-yj, lower.tail=FALSE)

    ii=which(ppj>pod1)
    if (length(ii)>0) {
      neli=c(neli,yj[min(ii)])
    } else {
      neli=c(neli, NA)
    }
  }
  neli[1:2]=NA
  deli=max(dose0)+100
  id1=nd+1

  id=1; Yi=NULL; dosei=NULL

  iec=1; outD=NULL; outD0=NULL; outP=NULL; Yj=NULL
  outDTL=NULL; outNd=NULL
  outN=outN0=numeric(nd)
  outM=outDLT=NULL
  fiti=NULL
  alpha.post=NULL; beta.post=NULL
  repeat {

    temi=dati[,id]; indexi=which(!is.na(temi)) # DLTs at current dose and indices of non-NA values
    dati[indexi[1:m],id]=NA # make the first m rows NA
    Yi=c(Yi,temi[indexi[1:m]]) # DLT's in cohort
    dosei=c(dosei,rep(dose0[id],m)) # The m doses taken

    outD0=c(outD0, id) # Dose taken
    outDLT=c(outDLT, sum(temi[indexi[1:m]])) # Sum DLTs
    outM=c(outM, m)

    outDTLi=numeric(nd) # Vector of n_d zeros
    outDTLi[id]=sum(temi[indexi[1:m]]) # Sum DLTs
    outDTL=rbind(outDTL, outDTLi)

    outNdi=numeric(nd) # Vector of n_d zeros
    outNdi[id]=m # Number of patients
    outNd=rbind(outNd, outNdi)

    outN[id]=outN[id]+sum(temi[indexi[1:m]])
    outN0[id]=outN0[id]+m

    alpha.post[id] <- alpha+outN[id]
    beta.post[id] <- beta+outN0[id]-outN[id]

    max=(pbeta(Equiv[2],alpha.post[id],beta.post[id])-pbeta(Equiv[1],alpha.post[id],beta.post[id]))/(Equiv[2]-Equiv[1])
    argmax <- "Equiv"

    for(i in 1:length(Under)){
      UMP <- (pbeta(Under[[i]][2],alpha.post[id],beta.post[id])-pbeta(Under[[i]][1],alpha.post[id],beta.post[id]))/(Under[[i]][2]-Under[[i]][1])
      if(UMP>max){
        max=UMP
        argmax <- "Under"
      }
    }

    for(i in 1:length(Over)){
      UMP <- (pbeta(Over[[i]][2],alpha.post[id],beta.post[id])-pbeta(Over[[i]][1],alpha.post[id],beta.post[id]))/(Over[[i]][2]-Over[[i]][1])
      if(UMP>max){
        max=UMP
        argmax <- "Over"
      }
    }

    if (outN[id]>=neli[outN0[id]] & !is.na(neli[outN0[id]])) {

      deli=c(deli, dose0[id])
      id=id-1
      outD=c(outD,id)

      id1=pmin(id,id1)

      if (id<=0) {
        stop_reason="All toxic"
        break
      }


    } else {

      if (argmax=="Under") {
        id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
        id=pmin(pmin(id+1,nd), id2)
        id=pmin(id,id1)

      } else if (argmax=="Over") {
        id=pmin(pmax(1,id-1), id1)
      }

      outD=c(outD,id)
    }


    if (sum(outN0)+m > n) {
      stop_reason="MaxN reached"

      while (outN0[id] < MinN_atMTD && id >= 2) {
        id=id-1
      }

      break
    }
    # if (iec>=nec) break
    # iec=iec+1
  }

  fiti=select.mtd(target=ttl, npts=outN0, ntox=outN)

  outD[length(outD)]=pmin(fiti$MTD,ifelse(id1==0,1000,id1))
  # outD[length(outD)]=fiti$MTD

  colnames(outDTL)=colnames(outNd)=paste0("P",1:nd)
  rownames(outDTL)=rownames(outNd)=1:nrow(outDTL)

  tmpi=outNd
  tmpi[tmpi>0]=paste0("/",tmpi[tmpi>0])
  outND=matrix(paste0(outDTL,tmpi),ncol=nd)
  colnames(outND)=paste0("P",1:nd)
  rownames(outND)=1:nrow(outDTL)
  outND[outND=="00"]=""; outND=as.data.frame(outND)

  outDTL[outDTL==0]=""; outDTL=as.data.frame(outDTL)
  outNd[outNd==0]=""; outNd=as.data.frame(outNd)


  outi2=rbind(outN, outN0)
  colnames(outi2)=paste0("P",1:nd)
  rownames(outi2)=c("DLT","Total")

  outi=data.frame(D0=outD0, D1=outD)


  outi3=cbind(Dose1=outD0, Dose=dose0[outD0], DLT=outDLT, M=outM)

  out=list(Y=Yi, dose=dosei, BOIN=fiti, outData=outi3, outN=outi2, path=outi, outND=outND, stop_reason=stop_reason)
  return(out)
}

RFun_33 = function(dati,n,m,dose0){
  nd=length(dose0)
  id=1; Yi=NULL; dosei=NULL
  iec=1; outD=NULL; outD0=NULL; outP=NULL; Yj=NULL
  outDTL=NULL; outNd=NULL
  outN=outN0=numeric(nd)
  outM=outDLT=NULL
  fiti=NULL
  stay=0

  repeat{
    temi=dati[,id]; indexi=which(!is.na(temi))
    dati[indexi[1:m],id]=NA
    Yi=c(Yi,temi[indexi[1:m]])
    dosei=c(dosei,rep(dose0[id],m))

    outD0=c(outD0, id)
    outDLT=c(outDLT, sum(temi[indexi[1:m]]))
    outM=c(outM, m)

    outDTLi=numeric(nd)
    outDTLi[id]=sum(temi[indexi[1:m]])
    outDTL=rbind(outDTL, outDTLi)

    outNdi=numeric(nd)
    outNdi[id]=m
    outNd=rbind(outNd, outNdi)

    outN[id]=outN[id]+sum(temi[indexi[1:m]])
    outN0[id]=outN0[id]+m

    if(stay==1){ # If stay=1 we are using the same dose as before
      if(outDTLi[id]==0){
        id=pmin(id+1,nd)
        stay <- 0
        outD=c(outD,id)
      }
      else if(outDTLi[id]>0){
        id <- pmax(id-1,1)
        outD=c(outD,id)
        MTD <- dose0[id]
        stop_reason <- "> 1 of 6 DLT"
        break
      }
    } else if(outDTLi[id]==0){
      id=pmin(id+1,nd)
      outD=c(outD,id)
    } else if(outDTLi[id]==1){
      id <- id
      stay <- 1
      outD=c(outD,id)
    } else if(outDTLi[id]>=2){
      id <- pmax(id-1,1)
      MTD <- dose0[id]
      outD=c(outD,id)
      stop_reason="2-3 DLT"
      break
    }

    if (sum(outN0)+m > n) {
      stop_reason="MaxN reached"
      break
    }
  }

  colnames(outDTL)=colnames(outNd)=paste0("P",1:nd)
  rownames(outDTL)=rownames(outNd)=1:nrow(outDTL)

  tmpi=outNd
  tmpi[tmpi>0]=paste0("/",tmpi[tmpi>0])
  outND=matrix(paste0(outDTL,tmpi),ncol=nd)
  colnames(outND)=paste0("P",1:nd)
  rownames(outND)=1:nrow(outDTL)
  outND[outND=="00"]=""; outND=as.data.frame(outND)

  outDTL[outDTL==0]=""; outDTL=as.data.frame(outDTL)
  outNd[outNd==0]=""; outNd=as.data.frame(outNd)


  outi2=rbind(outN, outN0)
  colnames(outi2)=paste0("P",1:nd)
  rownames(outi2)=c("DLT","Total")

  outi=data.frame(D0=outD0, D1=outD)


  outi3=cbind(Dose1=outD0, Dose=dose0[outD0], DLT=outDLT, M=outM)

  MTD = dose0[id]

  out=list(Y=Yi, dose=dosei, MTD=MTD, outData=outi3, outN=outi2, path=outi, outND=outND, stop_reason=stop_reason)
}

RFun_i33 = function(dati, n, m, dose0, ttl=0.3, epsilon1 = 0.05, epsilon2 = 0.05, pod1=0.95, MinN_atMTD=6, Max_FoldChange=3) {

  nd=length(dose0)
  Under=c(0,ttl-epsilon1)
  Equiv=c(ttl-epsilon1,ttl+epsilon2)
  Over=c(ttl+epsilon2,1)
  ## Safety bound
  ni=1:n; neli=NULL
  for (i in 1:length(ni)) {
    nj=ni[i]; yj=0:nj
    ppj=pbeta(ttl, 1+yj, 1+nj-yj, lower.tail=FALSE)

    ii=which(ppj>pod1)
    if (length(ii)>0) {
      neli=c(neli,yj[min(ii)])
    } else {
      neli=c(neli, NA)
    }
  }
  neli[1:2]=NA
  deli=max(dose0)+100
  id1=nd+1

  id=1; Yi=NULL; dosei=NULL

  iec=1; outD=NULL; outD0=NULL; outP=NULL; Yj=NULL
  outDTL=NULL; outNd=NULL
  outN=outN0=numeric(nd)
  outM=outDLT=NULL
  fiti=NULL
  alpha.post=NULL; beta.post=NULL
  repeat {

    temi=dati[,id]; indexi=which(!is.na(temi)) # DLTs at current dose and indices of non-NA values
    dati[indexi[1:m],id]=NA # make the first m rows NA
    Yi=c(Yi,temi[indexi[1:m]]) # DLT's in cohort
    dosei=c(dosei,rep(dose0[id],m)) # The m doses taken

    outD0=c(outD0, id) # Dose taken
    outDLT=c(outDLT, sum(temi[indexi[1:m]])) # Sum DLTs
    outM=c(outM, m)

    outDTLi=numeric(nd) # Vector of n_d zeros
    outDTLi[id]=sum(temi[indexi[1:m]])
    outDTL=rbind(outDTL, outDTLi)

    outNdi=numeric(nd)
    outNdi[id]=m
    outNd=rbind(outNd, outNdi)

    outN[id]=outN[id]+sum(temi[indexi[1:m]])
    outN0[id]=outN0[id]+m

    if (outN[id]>=neli[outN0[id]] & !is.na(neli[outN0[id]])) {

      deli=c(deli, dose0[id])
      id=id-1
      outD=c(outD,id)

      id1=pmin(id,id1) #highest candidate dose

      if (id<=0) {
        stop_reason="All toxic"
        break
      }
    } else{
      if(outN[id]/outN0[id]<Equiv[1]){
        id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
        id=min(c(id+1,nd,id2,id1))
      } else if(outN[id]/outN0[id]>=Equiv[1] & outN[id]/outN0[id]<=Equiv[2]){
        id=id
      } else if(outN[id]/outN0[id]>Equiv[2]){
        if((outN[id]-1)/outN0[id]<Equiv[1]){
          id=id
        } else{
          id=pmin(pmax(1,id-1), id1)
        }
      }
      outD=c(outD,id)
    }


    if (sum(outN0)+m > n) {
      stop_reason="MaxN reached"

      while (outN0[id] < MinN_atMTD && id >= 2) {
        id=id-1
      }

      break
    }
    # if (iec>=nec) break
    # iec=iec+1
  }

  fiti=select.mtd(target=ttl, npts=outN0, ntox=outN)

  outD[length(outD)]=pmin(fiti$MTD,ifelse(id1==0,1000,id1))
  # outD[length(outD)]=fiti$MTD

  colnames(outDTL)=colnames(outNd)=paste0("P",1:nd)
  rownames(outDTL)=rownames(outNd)=1:nrow(outDTL)

  tmpi=outNd
  tmpi[tmpi>0]=paste0("/",tmpi[tmpi>0])
  outND=matrix(paste0(outDTL,tmpi),ncol=nd)
  colnames(outND)=paste0("P",1:nd)
  rownames(outND)=1:nrow(outDTL)
  outND[outND=="00"]=""; outND=as.data.frame(outND)

  outDTL[outDTL==0]=""; outDTL=as.data.frame(outDTL)
  outNd[outNd==0]=""; outNd=as.data.frame(outNd)


  outi2=rbind(outN, outN0)
  colnames(outi2)=paste0("P",1:nd)
  rownames(outi2)=c("DLT","Total")

  outi=data.frame(D0=outD0, D1=outD)


  outi3=cbind(Dose1=outD0, Dose=dose0[outD0], DLT=outDLT, M=outM)

  out=list(Y=Yi, dose=dosei, BOIN=fiti, outData=outi3, outN=outi2, path=outi, outND=outND, stop_reason=stop_reason)
  return(out)
}

RFun_CS = function(dati, n, dose0, cohort.size, theta, ttl=0.3, epsilon1 = 0.05, epsilon2 = 0.05, pod1=0.95, MinN_atMTD=6, Max_FoldChange=3,alpha1=1,beta1=4) {

  b=rep(0,length(cohort.size))
  for(i in 1:length(cohort.size)){
    DLT<-0:cohort.size[i]
    post.prob <- 1-pbeta(theta,alpha1+DLT,beta1+cohort.size[i]-DLT)
    b[i]<-which(post.prob==min(post.prob[post.prob>0.1]))
  }

  nd=length(dose0)

  ## Safety bound
  ni=1:n; neli=NULL
  for (i in 1:length(ni)) {
    nj=ni[i]; yj=0:nj
    ppj=pbeta(ttl, 1+yj, 1+nj-yj, lower.tail=FALSE)

    ii=which(ppj>pod1)
    if (length(ii)>0) {
      neli=c(neli,yj[min(ii)])
    } else {
      neli=c(neli, NA)
    }
  }
  neli[1:2]=NA
  deli=max(dose0)+100
  id1=nd

  id=1; Yi=NULL; dosei=NULL

  iec=1; outD=NULL; outD0=NULL; outP=NULL; Yj=NULL
  outDTL=NULL; outNd=NULL
  outN=outN0=numeric(nd)
  outM=outDLT=NULL
  fiti=NULL
  alpha.post=NULL; beta.post=NULL
  curr=0
  cohort.id=1
  m=cohort.size[cohort.id]
  last.action="Start"
  situation=0
  repeat {

    temi=dati[,id]; indexi=which(!is.na(temi)) # DLTs at current dose and indices of non-NA values
    dati[indexi[1:m],id]=NA # make the first m rows NA
    Yi=c(Yi,temi[indexi[1:m]]) # DLT's in cohort
    dosei=c(dosei,rep(dose0[id],m)) # The m doses taken

    outD0=c(outD0, id) # Dose taken
    outDLT=c(outDLT, sum(temi[indexi[1:m]])) # Sum DLTs
    outM=c(outM, m)

    outDTLi=numeric(nd) # Vector of n_d zeros
    outDTLi[id]=sum(temi[indexi[1:m]])
    outDTL=rbind(outDTL, outDTLi)

    outNdi=numeric(nd)
    outNdi[id]=m
    outNd=rbind(outNd, outNdi)

    outN[id]=outN[id]+sum(temi[indexi[1:m]])
    outN0[id]=outN0[id]+m

    if (outN[id]>=neli[outN0[id]] & !is.na(neli[outN0[id]]) & situation!=1 & situation !=2) {

      deli=c(deli, dose0[id])
      id=id-1
      outD=c(outD,id)

      id1=pmin(id,id1)

      cohort.id=length(cohort.size)
      m=cohort.size[cohort.id]
      last.action = "De-escalate"

      if (id<=0) {
        stop_reason="All toxic"
        break
      }


    } else {

      if(situation==1){
        situation=0
        if(DLT.prev+outDTLi[id]<b[cohort.id]){
          if(last.action == "De-escalate"){
            id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
            id=min(c(id+1,nd,id2,id1))
            outD=c(outD,id)
            stop_reason="Safe dose after de-escalation"
            break
          }
          id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
          id=min(c(id+1,nd,id2,id1))
          last.action = "Escalate"
        }
        if(DLT.prev+outDTLi[id]==b[cohort.id]){
          outD=c(outD,id)
          stop_reason="Cannot enroll additional patients"
          break
        }
        if(DLT.prev+outDTLi[id]>b[cohort.id]){
          if(id==2){
            id=id-1
            outD=c(outD,id)
            stop_reason="De-escalate to lowest dose"
            break #de-escalating to the lowest dose
          }
          id=pmin(pmax(1,id-1), id1)
          last.action = "De-escalate"
        }
      } else if(situation==2){
        situation=0
        if(DLT.prev+outDTLi[id]<b[cohort.id+1]){
          if(last.action == "De-escalate"){
            id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
            id=min(c(id+1,nd,id2,id1))
            outD=c(outD,id)
            stop_reason="Safe dose after de-escalation"
            break
          }
          id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
          id=min(c(id+1,nd,id2,id1))
          cohort.id=cohort.id+1
          m=cohort.size[cohort.id]
          last.action = "Escalate"
        } else if(DLT.prev+outDTLi[id]>=b[cohort.id+1]){
          if(id==2){
            id=id-1
            outD=c(outD,id)
            stop_reason="De-escalate to lowest dose"
            break #de-escalating to the lowest dose
          }
          id=pmin(pmax(1,id-1), id1)
          cohort.id=length(cohort.size)
          m=cohort.size[cohort.id]
          last.action = "De-escalate"
        }
      } else if(last.action=="Escalate" & id==id1 & outDTLi[id]<b[cohort.id] & cohort.id<length(cohort.size)){
        cohort.id=length(cohort.size)
        m=cohort.size[cohort.id]-m
        DLT.prev = outDTLi[id]
        situation=1
      } else if(outDTLi[id]<b[cohort.id]){
        if(last.action == "De-escalate"){
          id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
          id=min(c(id+1,nd,id2,id1))
          outD=c(outD,id)
          stop_reason="Safe dose after de-escalation"
          break
        }
        id2=max(which(dose0/dose0[id] < Max_FoldChange + 0.001))
        id=min(c(id+1,nd,id2,id1))
        last.action = "Escalate"
      } else if(outDTLi[id]>b[cohort.id]){
        if(id==2){
          id=id-1
          outD=c(outD,id)
          stop_reason="De-escalate to lowest dose"
          break #de-escalating to the lowest dose
        }
        id=pmin(pmax(1,id-1), id1)
        cohort.id=length(cohort.size)
        m=cohort.size[cohort.id]
        last.action = "De-escalate"
      } else if(outDTLi[id]==b[cohort.id]){
        if(cohort.id==length(cohort.size)){
          outD=c(outD,id)
          stop_reason="Cannot enroll additional patients"
          break
        }
        m=cohort.size[cohort.id+1]-cohort.size[cohort.id]
        DLT.prev <- outDTLi[id]
        situation=2
      }

      outD=c(outD,id)

    }

    if (sum(outN0)+m > n) {
      stop_reason="MaxN reached"

      while (outN0[id] < MinN_atMTD && id >= 2) {
        id=id-1
      }

      break
    }
  }

  fiti=select.mtd(target=ttl, npts=outN0, ntox=outN)

  outD[length(outD)]=pmin(fiti$MTD,ifelse(id1==0,1000,id1))
  # outD[length(outD)]=fiti$MTD

  colnames(outDTL)=colnames(outNd)=paste0("P",1:nd)
  rownames(outDTL)=rownames(outNd)=1:nrow(outDTL)

  tmpi=outNd
  tmpi[tmpi>0]=paste0("/",tmpi[tmpi>0])
  outND=matrix(paste0(outDTL,tmpi),ncol=nd)
  colnames(outND)=paste0("P",1:nd)
  rownames(outND)=1:nrow(outDTL)
  outND[outND=="00"]=""; outND=as.data.frame(outND)

  outDTL[outDTL==0]=""; outDTL=as.data.frame(outDTL)
  outNd[outNd==0]=""; outNd=as.data.frame(outNd)


  outi2=rbind(outN, outN0)
  colnames(outi2)=paste0("P",1:nd)
  rownames(outi2)=c("DLT","Total")

  outi=data.frame(D0=outD0, D1=outD)


  outi3=cbind(Dose1=outD0, Dose=dose0[outD0], DLT=outDLT, M=outM)

  out=list(Y=Yi, dose=dosei, BOIN=fiti, outData=outi3, outN=outi2, path=outi, outND=outND, stop_reason=stop_reason)
  return(out)
}
