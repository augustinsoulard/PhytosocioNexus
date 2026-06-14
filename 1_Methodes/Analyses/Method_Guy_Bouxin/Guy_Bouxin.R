Caf.R

don<-as.matrix(filename)
nbrrel= ncol(don)
nbrsp= nrow(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
nrang= 1+nbrsp+nbrrel
nbrcol=nax*3
nomrang<-list(1:nrang)
nomcol<-list(1:nax*3)
for(j in 1:nax)
{
  nomcol[1+(j-1)*3]<-paste("coord",j)
  nomcol[2+(j-1)*3] <-paste("Cr%",j)
  nomcol[3+(j-1)*3] <-paste("P",j)
}
resulca<-matrix(1:nrang*nbrcol, nrow=nrang,ncol=nbrcol,byrow=TRUE)
colnames(resulca)<- nomcol
#calculation of the transformed and the covariance matrices
donc<- matrix(1 :nbrsp*nbrrel, nrow=nbrsp,ncol=nbrrel)
don1<-matrix(1 :nbrsp*nbrrel, nrow=nbrsp,ncol=nbrrel)
donc<-don
covar1<-matrix(rep(0,nbrrel*nbrrel),nbrrel,nbrrel)
covar1<-covarca(donc,nbrsp,nbrrel)
#distribution histogram of the covariance coefficients with the function hist()
nbrcel= nbrrel*(nbrrel-1)/2 +nbrrel
matcov<-as.matrix(1 :nbrcel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 1:nbrrel)
{
  for (j2 in 1:j1)
  {
    i=i+1
    matcov[i]=covar1[j1,j2]
  }
}
hist(matcov, breaks=40, include.lowest = TRUE, right = TRUE,
     labels=FALSE, border = NULL,
     main = paste("Matrice de covariance CA"), plot=TRUE,
     axes = TRUE,
     warn.unused = TRUE)
#calculation of the trace, of the eigenvectors and eigenvalues
Lam<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lam2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
R<-eigen(covar1, symmetric=TRUE, only.values=FALSE)
V <- R$vectors; Lam <- R$values
trace = sum(Lam)
trace
Lam2<-Lam*100/trace
#calculation of the coordinates and relative contributions of the individuals on axes
hr<-rowSums(don)
hc<-colSums(don)
somme=sum(don)
coordrel<-matrix(1:nbrrel*nax,nrow=nbrrel,ncol=nax)
contrerel<-matrix(1:nbrrel*nax,nrow=nbrrel,ncol=nax)
for(j in 1:nbrrel)
{
  for(b in 1:nax)
  {
    coordrel[j,b]=V[j,b]*sqrt(Lam[b])*sqrt(somme/hc[j])
    contrerel[j,b] = 100 * V[j,b] * V[j,b] 
  }
}
#calculation of the coordinates and relative contributions of the species on axes
coordsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
for (b in 1:nax)
{
  for(i in 1:nbrsp)
  {
    cosp=0
    for( j in 1:nbrrel)
    {
      cosp= cosp + don[i,j]*V[j, b]*sqrt(somme/hc[j])/hr[i]
    }
    coordsp[i,b]=cosp
    contrelsp[i, b]= cosp^2*100*hr[i] /(Lam[b]*somme)
  }
}
#permutations of the matrix and new calculations
prob1<-matrix(1:nax*1,nrow=1,ncol=nax,byrow=TRUE)
prob2<-matrix(1:nax*nbrrel,nrow=nbrrel,ncol=nax,byrow=TRUE)
prob3<-matrix(1:nax*nbrsp,nrow=nbrsp,ncol=nax,byrow=TRUE)
prob1[]<-0
prob2[]<-0
prob3[]<-0
nombreper=0
for(np in 1:nper)
{
  Lamp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  Lamp2<-matrix(rep(0,1*nbrrel),1,nbrrel)
  donc<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  don1<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  hc<-matrix(rep(0,1*nbrrel),1,nbrrel)
  covar1<-matrix(rep(0,nbrrel*nbrrel),nbrrel,nbrrel)
  for(i in 1 :nbrsp)
  {
    donc[i,]<-sample(don[i,],nbrrel,replace=FALSE)
  }
  hc<-colSums(donc)
  if (min(hc)==0){
    next
  }
  nombreper=nombreper+1
  covar1<-covarca(donc,nbrsp,nbrrel)
  #calculation of the eigenvectors, of the eigenvalues and of the trace after permutations
  Rp<-eigen(covar1, symmetric=TRUE, only.values=FALSE)
  Vp <- Rp$vectors; Lamp <- Rp$values
  trace2=sum(Lamp)
  Lamp2<-Lamp*100/trace2
  for(b in 1:nax)
  {
    if (Lamp[b]>=Lam[b])
    {
      prob1[b]=prob1[b]+1
    }
  }
  #calculation of the relative contributions of the individuals on axes after permutations
  contrerelp<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
  for(j in 1:nbrrel)
  {
    for(b in 1:nax)
    {
      contrerelp[j,b] = 100 * Vp[j,b] * Vp[j,b] 
      if(contrerelp[j,b] >=contrerel[j,b]){ 
        prob2[j,b]=prob2[j,b]+1
      }
    }
  }
  #calculation of the relative contributions of the species on axes after permutations
  contrelspp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
  for(b in 1:nax)
  {
    for (i in 1:nbrsp)
    {
      cosp=0
      for( j in 1:nbrrel)
      {
        cosp= cosp + donc[i,j]*Vp[j, b]*sqrt(somme/hc[j])/(hr[i])
      }
      contrelspp[i, b]= cosp^2*100*hr[i]/(Lamp[b]*somme)
      if (contrelspp[i,b]>=contrelsp[i,b])
      {
        prob3[i,b]=prob3[i,b]+1
      }
    }
  }
}
nombreper
prob1=prob1/nombreper*100
prob2=prob2/nombreper*100
prob3=prob3/nombreper*100
for(j in 1:nax)
{
  resulca[1,3*(j-1)+1]=Lam[j]
  resulca[1,3*(j-1)+2]=Lam2[j]
  resulca[1,3*(j-1)+3]=prob1[j]
}
nomrang[1]<- paste("Vp,Vp% et P")
for(i in 1 :nbrsp)
{
  for(j in 1:nax)
  {
    resulca[i+1,3*(j-1)+1]=coordsp[i,j]
    resulca[i+1,3*(j-1)+2]=contrelsp[i,j]
    resulca[i+1,3*(j-1)+3]=prob3[i,j]
  }
  nomrang[i+1]=nomsp[i]
}
for(i in 1 :nbrrel)
{
  for(j in 1:nax)
  {
    resulca[i+1+nbrsp,3*(j-1)+1]=coordrel[i,j]
    resulca[i+1+nbrsp,3*(j-1)+2]=contrerel[i,j]
    resulca[i+1+nbrsp,3*(j-1)+3]=prob2[i,j]
  }
  nomrang[i+1+nbrsp]=nomrel[i]
}
rownames(resulca)<-nomrang



covarca.R

covarca<-function(donc,nbrsp,nbrrel)
{
  somme<-0  
  hr<-rowSums(donc)
  hc<-colSums(donc)
  somme=sum(donc)
  for(i in 1:nbrsp) 
  {     
    for(j in 1 :nbrrel)
    {   
      don1[i,j]=donc[i,j]*(somme/(hc[j]*hr[i]))-1
    } 
  }
  for(j1 in 1:nbrrel)
  {
    for(j2 in 1:nbrrel)
    {
      for (i in 1:nbrsp)
      {
        covar1[j1,j2]= covar1[j1,j2]+don1[i,j1]*don1[i,j2]*hr[i]*sqrt(hc[j1]*hc[j2])
      }
      covar1[j1,j2]= covar1[j1,j2]/somme^2
    }
  }
  covar1
}

covarnsca.R

covarnsca<-function(donc,nbrsp,nbrrel)
{
  somme = 0
  hr<-rowSums(donc)
  hc<-colSums(donc)
  somme=sum(donc)
  for(i in 1:nbrsp)
  { 
    for(j in 1:nbrrel)
    { 
      don1[i,j]= (donc[i,j]/hc[j]-hr[i]/somme)*nbrsp
    }
  }
  for(j1 in 1:nbrrel)
  {
    for(j2 in 1:nbrrel)
    {
      for (i in 1:nbrsp)
      {
        covar1[j1,j2] = covar1[j1,j2] +don1[i,j1]*don1[i,j2]*sqrt(hc[j1]*hc[j2])
      }
      covar1[j1,j2]= covar1[j1,j2]/(nbrsp*somme)
    }
  }
  covar1
}

Cs.R

don<-as.matrix(filename)
ng=length(co)
nbrsp= nrow(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
nbrrel2=0
nbrrel=0
for(i in 1:ng)
{
  nbrrel2=nbrrel2+co[i]+3
  nbrrel=nbrrel+co[i]
}
nbrrel2=nbrrel2+2
sligne<-rowSums(don)
nomcol1<-list(1:nbrrel2)
#Writing the frequencies of the complete table 
res<-matrix(1 :nbrrel2*nbrsp,nrow=nbrsp,ncol=nbrrel2,byrow=TRUE)
for(i in 1:nbrsp)
{
  res[i,1]=sligne[i]
  res[i,2]=sligne[i]/nbrrel
}
nomcol1[1]<-paste("sum")
nomcol1[2]<-paste("fr")
#Frequency calculation in each cluster
k=0 ;l=0
sumg<-matrix(rep(0,nbrsp*ng), nbrsp,ng)
frg<-matrix(rep(0,nbrsp*ng), nbrsp,ng)
for(p in 1:ng)
{
  nrg=co[p]
  k=l+1
  l=l+nrg
  for(i in 1:nbrsp)
  {
    for(j in k: l)
    {
      sumg[i,p]=sumg[i,p]+don[i,j] 
    }
    frg[i,p]=sumg[i,p]/nrg
  }
}
#Line permutation
prg<-matrix(rep(0,nbrsp*ng), nbrsp,ng)
nbrrel
for(np in 1:nper)
{
  donp<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  for(i in 1:nbrsp)
  {
    donp[i,]<-sample(don[i,],nbrrel,replace=FALSE)
  }
  #For each line
  for(i in 1 :nbrsp)
  {
    #For each cluster
    #nrg : releve number in each cluster
    k=0 ; l=0
    for(g in 1 : ng)
    {
      nrg=co[g]
      k=l+1 ;l=l+nrg
      sumgp=0 
      for(j in k :l)
      {
        sumgp=sumgp+donp[i,j]
      }
      if (sumgp>=sumg[i,g])
      {
        prg[i,g]=prg[i,g]+1
      }
    }
  }
}
for(i in 1 :nbrsp)
{
  for(g in 1 :ng)
  {
    prg[i,g]=prg[i,g]/nper
  }
}
#Writing the results 
for(g in 1:ng)
{
  if(g==1)
  {
    nbrrel1=1
    nbrrel2=co[1]
  }
  else
  {
    nbrrel1=nbrrel1+co[g-1]
    nbrrel2=nbrrel2+co[g]
  }
  for(i in 1:nbrsp)
  {
    for(j in nbrrel1: nbrrel2)
    {
      nomcol1[j+2+(g-1)*3]=nomrel[j]
      res[i,j+2+(g-1)*3]=don[i,j]
    }
    res[i,nbrrel2+3+(g-1)*3]=sumg[i,g]
    res[i,nbrrel2+4+(g-1)*3]=frg[i,g]
    res[i,nbrrel2+5+(g-1)*3]=prg[i,g]
  }
  nomcol1[nbrrel2+3+(g-1)*3]<-paste("sum")
  nomcol1[nbrrel2+4+(g-1)*3]<-paste("fr")
  nomcol1[nbrrel2+5+(g-1)*3]<-paste("pr")
}
rownames(res)<-nomsp
colnames(res)<-nomcol1
resul<-data.frame(res)


disj113.R

don<-as.matrix(filename)
nbrrel= ncol(don)
nbrsp= nrow(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
#calcul nombre de lignes du nouveau tableau
nbrl=0
don1<-don
res3<-apply(don,1,max)
for(i in 1:nbrsp)
{
  if(res3[i]>3)
  {
    nbrl=nbrl+1
    for(j in 1:nbrrel)
    {
      if(don[i,j]>3)
      {
        don[i,j]=0
      }
    }
  }
}
res2<-apply(don,1,max)
for(i in 1:nbrsp)
{
  if(res2[i]>1)
  {
    nbrl=nbrl+1
    for(j in 1:nbrrel)
    {
      if(don[i,j]>1)
      {
        don[i,j]=0
      }
    }
  }
}
res1<- apply(don,1,max)
for(i in 1:nbrsp)
{
  if(res1[i]==1)
  {
    nbrl=nbrl+1
  }
}
nbrl

#definition of a new table
nomsp2<-list(1 :nbrl)
don2<-matrix(rep(0,nbrl*nbrrel),nbrl,nbrrel)
nbrl=0
for(i in 1:nbrsp)
{
  if(res3[i]>3)
  {
    nbrl=nbrl+1
    nomsp2[nbrl]<-paste(nomsp[i],">3")
    
    for(j in 1:nbrrel)
    {
      if(don1[i,j]>3)
      {
        don2[nbrl,j]=1
      }
    }
  }
  if(res2[i]>1)
  {
    nbrl=nbrl+1
    nomsp2[nbrl]<-paste(nomsp[i],">1")
    
    for(j in 1:nbrrel)
    {
      if(don1[i,j]>1)
      {
        don2[nbrl,j]=1
      }
    }
  }
  if(res1[i]==1)
  {
    nbrl=nbrl+1
    nomsp2[nbrl]<-paste(nomsp[i])
    for(j in 1:nbrrel)
    {
      if(don1[i,j]>0)
      {
        don2[nbrl,j]=1
      }
    }
  }
}
nbrl
rownames(don2)<-c(nomsp2)
colnames(don2)<-c(nomrel)
resul<-data.frame(don2)

disj12345.R

don<-as.matrix(filename)
nbrrel= ncol(don)
nbrsp= nrow(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
#calcul nombre de lignes du nouveau tableau
nbrl=0
for(i in 1:nbrsp)
{
  for(j in 1:nbrrel)
  {
    if (don[i,j] == 1)
    {
      nbrl=nbrl+1
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] ==2)
    {
      nbrl=nbrl+1
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] ==3)
    {
      nbrl=nbrl+1
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] == 4)
    {
      nbrl=nbrl+1
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] == 5)
    {
      nbrl=nbrl+1
      break
    }
  }
}
nbrl
#definition of a new table
nomsp2<-list(1 :nbrl)
don1<-matrix(rep(0,nbrl*nbrrel),nbrl,nbrrel)
nbrl=0
for(i in 1:nbrsp)
{
  for(j in 1:nbrrel)
  {
    if (don[i,j] == 1)
    {
      nbrl=nbrl+1
      nomsp2[nbrl]<-paste(nomsp[i])
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] == 1)
    {
      don1[nbrl,j] = 1
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] ==2)
    {
      nbrl=nbrl+1
      nomsp2[nbrl]<-paste(nomsp[i],"2")
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] ==2)
    {
      don1[nbrl,j] =1
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==3)
    {
      nbrl=nbrl+1
      nomsp2[nbrl]<-paste(nomsp[i] ,"3")
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==3)
    {
      don1[nbrl,j]=1
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==4)
    {
      nbrl=nbrl+1
      nomsp2[nbrl]<-paste(nomsp[i] ,"4")
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==4)
    {
      don1[nbrl,j]=1
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==5)
    {
      nbrl=nbrl+1
      nomsp2[nbrl]<-paste(nomsp[i] ,"5")
      break
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==5)
    {
      don1[nbrl,j]=1
    }
  }
}
nbrl
rownames(don1)<-c(nomsp2)
colnames(don1)<-c(nomrel)
resul<-data.frame(don1)

mfaca.R

don<-as.matrix(filename)
nomsp<-rownames(don)
nomrel<-colnames(don)
nbrsp= nrow(don)
nbrsp1=nbrsp
nbrrel= ncol(don)
nb=length(co)
nrang=nb+2+nbrsp+nbrrel
nbrcol=nax*2
nomrang<-list(1:nrang)
nomcol<-list(1:nbrcol)
for(j in 1:nax)
{
  nomcol[1+(j-1)*2]<-paste("coord",j)
  nomcol[2+(j-1)*2] <-paste("Cr%",j)
}
resulmfaca<-matrix(1:1*nbrcol*nrang,nrow=nrang,ncol=nbrcol,byrow=TRUE)
colnames(resulmfaca)<-nomcol
#Eigenvalues of the originel table
donc<-matrix(1:nbrsp*nbrrel,nrow=nbrsp,ncol=nbrrel)
donc<-don
hr<-matrix(1:nbrsp)
hc<-matrix(1:nbrrel)
hr<-rowSums(don)
hc<-colSums(don)
somme<-sum(don)
covar1<-matrix(rep(0,nbrrel*nbrrel),nbrrel,nbrrel)
covar1<-covarca(donc,nbrsp,nbrrel)
Lam<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lam2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
R<-eigen(covar1,symmetric=TRUE)
V<- R$vectors
Lam<-R$values
trace=sum(Lam)
Lam2<-Lam*100/trace
for(j in 1 :nax)
{
  resulmfaca[1,2*(j-1)+1]=Lam[j]
  resulmfaca[1,2*(j-1)+2]=Lam2[j]
}
nomrang[1]<-paste("Vp & Vp% of the original table")
#Sub-tables and 1st eigenvalues of the sub-tables
v1st<-matrix(1 :nb)
i1=0
for(k in 1:nb)
{
  nbrsp= co[k]
  donc<-matrix(1 :nbrsp*nbrrel, nrow = nbrsp,ncol = nbrrel)
  for(i in 1 :nbrsp)
  {
    i1=i1+1
    for(j in 1 :nbrrel)
    {
      donc[i,j]=don[i1,j]
    }
  }
  covarst<-matrix(rep(0,nbrrel*nbrrel), nbrrel,nbrrel)
  covarst<-covarca(donc,nbrsp,nbrrel)
  Lamst<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
  Lamst2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
  Rst<-eigen(covarst, symmetric=TRUE, only.values=FALSE)
  Lamst<- Rst$values
  tracest = sum(Lamst)
  Lamst2=Lamst*100/tracest
  v1st[k]=Lamst[1]
  for(j in 1 :nax)
  {
    resulmfaca[1+k,2*(j-1)+1]=Lamst[j]
    resulmfaca[1+k,2*(j-1)+2]=Lamst2[j]
  }
  nomrang[1+k]<-paste("Vp & Vp% du sous-tableau",k)
}
#Weighting of the sub-tables by the 1st eigenvalues of the sub-tables
nbrsp=nbrsp1
donc<-matrix(1:nbrsp*nbrrel,nrow=nbrsp,ncol=nbrrel)
v1<-matrix(1 :nbrsp)
i1=0
for(k in 1:nb)
{
  nbrspp= co[k]
  for(i in 1:nbrspp)
  {
    i1=i1+1
    v1[i1]=v1st[k]
  }
}
for(i in 1:nbrsp)
{
  for(j in 1:nbrrel)
  {
    donc[i,j]=don[i,j]/sqrt(v1[i]) 
  }
}
covar1p<-matrix(1 :nbrrel*nbrrel,nrow = nbrrel,ncol = nbrrel)
covar1p<-covarca(donc,nbrsp,nbrrel)
#Distribution histogram of the covariance coefficients with the function hist()
nbrcel=nbrrel*(nbrrel-1)/2+nbrrel
matcov<-as.matrix(1 :nbrcel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 1:nbrrel)
{
  for (j2 in 1:j1)
  {
    i=i+1
    matcov[i]=covar1p[j1,j2]
  }
}
hist(matcov, breaks=40, include.lowest = TRUE, right = TRUE,
     border = NULL,
     main = paste("Covariance distribution-mfaca"), plot=TRUE,
     axes = TRUE,
     warn.unused = TRUE)
Lamp<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lamp2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Rp<-eigen(covar1p, symmetric=TRUE, only.values=FALSE)
Vp<- Rp$vectors; Lamp <- Rp$values
tracep= sum(Lamp)
Lamp2<-Lamp*100/tracep
for(j in 1 :nax)
{
  resulmfaca[2+nb,2*(j-1)+1]=Lamp[j]
  resulmfaca[2+nb,2*(j-1)+2]=Lamp2[j]
}
nomrang[2+nb]<-paste("Vp & Vp% of the weighted table")
coordsp<-matrix(rep(0,nbrsp*nax),nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
#Calculation of the coordinates and relative contributions of the individuals, on axes
coordrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
contrelrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
hr<-rowSums(donc)
hc<-colSums(donc)
somme=sum(donc)
for(j in 1:nbrrel)
{
  for(b in 1 :nax)
  {
    coordrel[j,b]=Vp[j,b]*sqrt(Lamp[b])/sqrt(somme/ hc[j])
    contrelrel[j,b]=100* Vp[j,b]* Vp[j,b]
  }
}
#Calculation of the coordinates and relative contributions of the species, on axes
coordsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
for(b in 1 :nax)
{
  for(i in 1 : nbrsp)
  {
    cosp=0
    for(j in 1 :nbrrel)
    {
      cosp= cosp + donc[i,j]*Vp[j,b]*sqrt(somme/hc[j])/hr[i]
    }
    coordsp[i,b]=cosp
    contrelsp[i,b]= cosp^2*100*hr[i]/(Lamp[b]*somme)
  }
}
for(i in 1 :nbrsp)
{
  nomrang[i+2+nb]<-nomsp[i]
  for(j in 1:nax)
  {
    resulmfaca[i+2+nb,2*(j-1)+1]=coordsp[i,j]
    resulmfaca[i+2+nb,2*(j-1)+2]=contrelsp[i,j]
  }
}
for(i in 1 :nbrrel)
{
  for(j in 1:nax)
  {
    resulmfaca[i+2+nb+nbrsp,2*(j-1)+1]=coordrel[i,j]
    resulmfaca[i+2+nb+nbrsp,2*(j-1)+2]=contrelrel[i,j]
  }
  nomrang[i+2+nb+nbrsp]<-nomrel[i]
}
rownames(resulmfaca)<-nomrang


mfanscafdensity.R

don<-as.matrix(filename)
nomsp<-rownames(don)
nomrel<-colnames(don)
nbrsp= nrow(don)
nbrrel= ncol(don)
somme=sum(don)
hr<-matrix(1:nbrsp)
hc<-matrix(1:nbrrel)
hr<-rowSums(don)
hc<-colSums(don)
nb=length(co)
nrang1=nb+1
nrang2=1+nbrsp+nbrrel
ncol1=2
ncol2=nax*3
nomrang1<-list(1:nrang1)
nomrang2<-list(1:nrang2)
nomcol1<-list(1:ncol1)
nomcol2<-list(1:ncol2)
nomcol1[1]<-paste("line numbers")
nomcol1[2]<-paste("table densities")
resul1mfansca<-matrix(1:nrang1*ncol1,nrow=nrang1,ncol=ncol1,byrow=TRUE)
resul2mfansca<-matrix(1:nrang2*ncol2,nrow=nrang2,ncol=ncol2,byrow=TRUE)
#Eigenvalues of the originel table
donc<- matrix(1 :nbrsp*nbrrel, nrow=nbrsp,ncol=nbrrel)
don1<-matrix(1 :nbrsp*nbrsp, nrow=nbrsp,ncol=nbrsp)
donc<-don
covar1<-matrix(rep(0,nbrrel*nbrrel),nbrrel,nbrrel)
covar1<-covarnsca(donc,nbrsp,nbrrel)
Lam<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lam2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
density<-matrix(1 :nb+1)
R<-eigen(covar1,symmetric=TRUE)
Lam<-R$values
trace<-sum(Lam)
density[1]<-somme/(nbrsp*nbrrel)
Lam2<-Lam*100/trace
resul1mfansca[1,1]<-nbrsp
resul1mfansca[1,2]<-density[1]
nomrang1[1]<-paste("Complete table")
#densities of sub-tables
v1st<-matrix(1 :nb)
i1=0
for(k in 1:nb)
{
  nbrspp= co[k]
  donst<-matrix(1 :nbrspp*nbrrel, nrow = nbrspp,ncol = nbrrel)
  for(i in 1 :nbrspp)
  {
    i1=i1+1
    for(j in 1 :nbrrel)
    {
      donst[i,j]=don[i1,j]
    }
  }
  sommest=sum(donst)
  density[k+1]<-sommest/(nbrrel*nbrspp)
  v1st[k]=density[k+1]
  resul1mfansca[1+k,1]<-nbrspp
  resul1mfansca[1+k,2]<-density[k+1]
  nomrang1[k+1]<-paste("subtable" ,k)
}
#Weighting of the original tables by the densities of the sub-tables
donc<-matrix(1:nbrsp*nbrrel,nrow=nbrsp,ncol=nbrrel)
v1<-matrix(1 :nbrsp)
i1=0
for(k in 1:nb)
{
  nbrspp= co[k]
  for(i in 1:nbrspp)
  {
    i1=i1+1
    v1[i1]=v1st[k]
  }
}
for(i in 1:nbrsp)
{ 
  for(j in 1:nbrrel)
  { 
    donc[i,j]=don[i,j]/sqrt(v1[i])
  }
}
covar1p<-matrix(1 :nbrrel*nbrrel,nrow = nbrrel,ncol = nbrrel)
covar1p<-covarnsca(donc,nbrsp,nbrrel)
somme = sum(donc)
hr<-rowSums(donc)
hc<-colSums(donc)
#Distribution histogram of the covariance - coefficients with the function hist()
nbrcel=nbrrel*(nbrrel-1)/2 + nbrrel
matcov<-as.matrix(1 :nbrrel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 2:nbrrel)
{
  for (j2 in 2:j1)
  {
    i=i+1
    matcov[i]=covar1p[j1,j2]
  }
}
hist(matcov, breaks=40, include.lowest = TRUE, right = TRUE,
     labels=FALSE, border = NULL,
     main = paste("Covariance distribution-mfansca"), plot=TRUE,
     axes = TRUE,
     warn.unused = TRUE)
Lamp<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lamp2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Rp<-eigen(covar1p, symmetric=TRUE, only.values=FALSE)
Vp<- Rp$vectors; Lamp <- Rp$values
tracep= sum(Lamp)
Lamp2<-Lamp*100/tracep
nomrang2[2+nb]<-paste("Vp & Vp% of weighted tables")
coordsp<-matrix(rep(0,nbrsp*nax),nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
#Calculation of the coordinates and relative contributions of the individuals, on axes
coordrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
contrelrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
for(j in 1:nbrrel)
{
  for(b in 1 :nax)
  {
    coordrel[j,b]=Vp[j,b]/sqrt(hc[j])*sqrt(somme)
    contrelrel[j,b]=100* Vp[j,b]* Vp[j,b]
  }
}
#Calculation of the coordinates and relative contributions of the species, on axes
coordsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
for(b in 1 :nax)
{
  for(i in 1 : nbrsp)
  {
    cosp=0
    const= sqrt(nbrsp/(Lamp[b]*somme))
    const
    for(j in 1 :nbrrel)
    {
      cosp= cosp + donc[i,j]*Vp[j, b]/sqrt(hc[j])*const
    }
    coordsp[i,b]=cosp
    contrelsp[i,b]= cosp^2*100
  }
}
#permutations of the weighted matrix and new calculations
prob1<-matrix(1:nax*1,nrow=1,ncol=nax,byrow=TRUE)
prob2<-matrix(1:nax*nbrrel,nrow=nbrrel,ncol=nax,byrow=TRUE)
prob3<-matrix(1:nax*nbrsp,nrow=nbrsp,ncol=nax,byrow=TRUE)
prob1[]<-0
prob2[]<-0
prob3[]<-0
nombreper<-0
for(np in 1:nper)
{
  Lampp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  Lampp2<-matrix(rep(0,1*nbrrel),1,nbrrel)
  donpp<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  don1pp<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  hcpp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  covar1pp<-matrix(1 :nbrrel*nbrrel,nrow=nbrrel,ncol = nbrrel)
  for(i in 1 :nbrsp)
  {
    donpp[i,]<-sample(donc[i,],nbrrel,replace=FALSE)
  }
  hcpp<-colSums(donpp)
  if(min(hcpp)==0){
    next
  }
  nombreper=nombreper+1
  donc<-donpp
  covar1pp<-covarnsca(donc,nbrsp,nbrrel)
  #calculation  of the eigenvectors, of the eigenvalues and of the trace after permutations
  Rpp<-eigen(covar1pp, symmetric=TRUE, only.values=FALSE)
  Vpp <- Rpp$vectors; Lampp <- Rpp$values
  tracepp<-sum(Lampp)
  Lampp2<-Lampp*100/tracepp
  for(b in 1 :nax)
  {
    if (Lampp[b]>=Lamp[b])
    {
      prob1[b]=prob1[b]+1
    }
  }
  #calculation of the relative contributions of the individuals on axes after permutations
  contrerelpp<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
  for(j in 1:nbrrel)
  {
    for(b in 1 :nax)
    {
      contrerelpp[j,b]=100* Vpp[j,b]* Vpp[j,b]
      if(contrerelpp[j,b] >=contrelrel[j,b]){ 
        prob2[j,b]=prob2[j,b]+1
      }
    }
  }
  #calculation of the relative contributions of the species on axes after permutations
  contrelspp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
  for(b in 1 :nax)
  {
    for(i in 1 : nbrsp)
    {
      cosp=0
      const= sqrt(nbrsp/(Lampp[b]*somme))
      for(j in 1 :nbrrel)
      {
        cosp= cosp + donpp[i,j]*Vpp[j, b]/sqrt(hcpp[j])*const
      }
      contrelspp[i,b]=cosp^2*100
      if (contrelspp[i,b] >=contrelsp[i,b]){
        prob3[i,b]=prob3[i,b]+1
      }
    }
  }
}
nombreper
prob1=prob1*100/nombreper
prob2=prob2*100/ nombreper
prob3=prob3*100/ nombreper
for(j in 1:nax)
{
  nomcol2[1+(j-1)*3]<-paste("coord",j)
  nomcol2[2+(j-1)*3] <-paste("CR",j)
  nomcol2[3+(j-1)*3] <-paste("P",j)
}
nomrang2[1]<-paste("Vp, Vp%, P")
for(j in 1:nax)
{
  resul2mfansca[1,3*(j-1)+1]=Lamp[j]
  resul2mfansca[1,3*(j-1)+2]=Lamp2[j]
  resul2mfansca[1,3*(j-1)+3]=prob1[j]
}
for(i in 1 :nbrsp)
{
  nomrang2[i+1]=nomsp[i]
  for(j in 1:nax)
  {
    resul2mfansca[i+1,3*(j-1)+1]=coordsp[i,j]
    resul2mfansca[i+1,3*(j-1)+2]=contrelsp[i,j]
    resul2mfansca[i+1,3*(j-1)+3]=prob3[i,j]
  }
}
for(i in 1 :nbrrel)
{
  nomrang2[i+1+nbrsp]=nomrel[i]
  for(j in 1:nax)
  {
    resul2mfansca[i+1+nbrsp,3*(j-1)+1]=coordrel[i,j]
    resul2mfansca[i+1+nbrsp,3*(j-1)+2]=contrelrel[i,j]
    resul2mfansca[i+1+nbrsp,3*(j-1)+3]=prob2[i,j]
  }
}
rownames(resul1mfansca)<-nomrang1
colnames(resul1mfansca)<-nomcol1
rownames(resul2mfansca)<-nomrang2
colnames(resul2mfansca)<-nomcol2

mfanscafEV.R

don<-as.matrix(filename)
nomsp<-rownames(don)
nomrel<-colnames(don)
nbrsp= nrow(don)
nbrsp1=nbrsp
nbrrel= ncol(don)
somme=sum(don)
don<-log1p(don)
hr<-matrix(1:nbrsp)
hc<-matrix(1:nbrrel)
hr<-rowSums(don)
hc<-colSums(don)
nb=length(co)
nrang1=nb+1
nrang2=1+nbrsp+nbrrel
ncol1=2
ncol2=nax*3
nomrang1<-list(1:nrang1)
nomrang2<-list(1:nrang2)
nomcol1<-list(1:ncol1)
nomcol2<-list(1:ncol2)
nomcol1[1]<-paste("line numbers")
nomcol1[2]<-paste("eigenvalues")
resul1mfansca<-matrix(1:nrang1*ncol1,nrow=nrang1,ncol=ncol1,byrow=TRUE)
resul2mfansca<-matrix(1:nrang2*ncol2,nrow=nrang2,ncol=ncol2,byrow=TRUE)
#Eigenvalues of the originel table
donc<- matrix(1 :nbrsp*nbrrel, nrow=nbrsp,ncol=nbrrel)
don1<- matrix(1 :nbrsp*nbrrel, nrow=nbrsp,ncol=nbrrel)
donc<-don
somme=0
hr<-rowSums(donc)
hc<-colSums(donc)
somme=sum(donc)
covar1<-covarnsca(donc,nbrsp,nbrrel)
Lam<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lam2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
R<-eigen(covar1,symmetric=TRUE)
Lam<-R$values
trace<-sum(Lam)
Lam2<-Lam*100/trace
resul1mfansca[1,1]<-nbrsp
resul1mfansca[1,2]<-Lam[1]
nomrang1[1]<-paste("Complete table")
#Sub-tables and 1st eigenvalues of the sub-tables
v1st<-matrix(1 :nb)
i1=0[]
for(k in 1:nb)
{
  nbrsp= co[k]
  donc<-matrix(1 :nbrsp*nbrrel, nrow = nbrsp,ncol = nbrrel)
  for(i in 1 :nbrsp)
  {
    i1=i1+1
    for(j in 1 :nbrrel)
    {
      donc[i,j]=don[i1,j]
    }
  }
  covarst<-matrix(rep(0,nbrrel*nbrrel), nbrrel,nbrrel)
  covarst<-covarnsca(donc,nbrsp,nbrrel)
  Lamst<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
  Lamst2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
  Rst<-eigen(covarst, symmetric=TRUE, only.values=FALSE)
  Lamst<- Rst$values
  tracest = sum(Lamst)
  Lamst2=Lamst*100/tracest
  v1st[k]=Lamst[1]
  for(j in 1 :nax)
  {
    resul1mfansca[1+k,1]=nbrsp
    resul1mfansca[1+k,2]=Lamst[1]
  }
  nomrang[1+k]<-paste("Vp & Vp% du sous-tableau",k)
}
v1st
#Weighting of the original tables by the eigenvalues of the sub-tables
nbrsp=nbrsp1
donc<-matrix(1:nbrsp*nbrrel,nrow=nbrsp,ncol=nbrrel)
v1<-matrix(1 :nbrsp)
i1=0
for(k in 1:nb)
{
  nbrspp= co[k]
  for(i in 1:nbrspp)
  {
    i1=i1+1
    v1[i1]=v1st[k]
  }
}
for(i in 1:nbrsp)
{ 
  for(j in 1:nbrrel)
  { 
    donc[i,j]=don[i,j]/sqrt(v1[i])
  }
}
covar1p<-matrix(1 :nbrrel*nbrrel,nrow = nbrrel,ncol = nbrrel)
covar1p<-covarnsca(donc,nbrsp,nbrrel)
#Distribution histogram of the covariance - coefficients with the function hist()
nbrcel=nbrrel*(nbrrel-1)/2 + nbrrel
matcov<-as.matrix(1 :nbrrel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 2:nbrrel)
{
  for (j2 in 2:j1)
  {
    i=i+1
    matcov[i]=covar1p[j1,j2]
  }
}
hist(matcov, breaks=40, include.lowest = TRUE, right = TRUE,
     labels=FALSE, border = NULL,
     main = paste("Covariance distribution-mfansca"), plot=TRUE,
     axes = TRUE,
     warn.unused = TRUE)
Lamp<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lamp2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Rp<-eigen(covar1p, symmetric=TRUE, only.values=FALSE)
Vp<- Rp$vectors; Lamp <- Rp$values
tracep= sum(Lamp)
Lamp2<-Lamp*100/tracep
nomrang2[2+nb]<-paste("Vp & Vp% of weighted tables")
coordsp<-matrix(rep(0,nbrsp*nax),nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
#Calculation of the coordinates and relative contributions of the individuals, on axes
coordrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
contrelrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
for(j in 1:nbrrel)
{
  for(b in 1 :nax)
  {
    coordrel[j,b]=Vp[j,b]/sqrt(hc[j])*sqrt(somme)
    contrelrel[j,b]=100* Vp[j,b]* Vp[j,b]
  }
}
#Calculation of the coordinates and relative contributions of the species, on axes
coordsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
for(b in 1 :nax)
{
  for(i in 1 : nbrsp)
  {
    cosp=0
    const= sqrt(nbrsp/(Lamp[b]*somme))
    const
    for(j in 1 :nbrrel)
    {
      cosp= cosp + donc[i,j]*Vp[j, b]/sqrt(hc[j])*const
    }
    coordsp[i,b]=cosp
    contrelsp[i,b]= cosp^2*100
  }
}
#permutations of the weighted matrix and new calculations
prob1<-matrix(1:nax*1,nrow=1,ncol=nax,byrow=TRUE)
prob2<-matrix(1:nax*nbrrel,nrow=nbrrel,ncol=nax,byrow=TRUE)
prob3<-matrix(1:nax*nbrsp,nrow=nbrsp,ncol=nax,byrow=TRUE)
prob1[]<-0
prob2[]<-0
prob3[]<-0
nombreper<-0
for(np in 1:nper)
{
  Lampp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  Lampp2<-matrix(rep(0,1*nbrrel),1,nbrrel)
  donpp<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  don1pp<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  hcpp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  covar1pp<-matrix(1 :nbrrel*nbrrel,nrow=nbrrel,ncol = nbrrel)
  for(i in 1 :nbrsp)
  {
    donpp[i,]<-sample(don[i,],nbrrel,replace=FALSE)
  }
  hcpp<-colSums(donpp)
  if(min(hcpp)==0){
    next
  }
  nombreper=nombreper+1
  donc<-donpp
  covar1pp<-covarnsca(donc,nbrsp,nbrrel)
  #calculation  of the eigenvectors, of the eigenvalues and of the trace after permutations
  Rpp<-eigen(covar1pp, symmetric=TRUE, only.values=FALSE)
  Vpp <- Rpp$vectors; Lampp <- Rpp$values
  tracepp<-sum(Lampp)
  Lampp2<-Lampp*100/tracepp
  for(b in 1 :nax)
  {
    if (Lampp[b]>=Lamp[b])
    {
      prob1[b]=prob1[b]+1
    }
  }
  #calculation of the relative contributions of the individuals on axes after permutations
  contrerelpp<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
  for(j in 1:nbrrel)
  {
    for(b in 1 :nax)
    {
      contrerelpp[j,b]=100* Vpp[j,b]* Vpp[j,b]
      if(contrerelpp[j,b] >=contrelrel[j,b]){ 
        prob2[j,b]=prob2[j,b]+1
      }
    }
  }
  #calculation of the relative contributions of the species on axes after permutations
  contrelspp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
  for(b in 1 :nax)
  {
    for(i in 1 : nbrsp)
    {
      cosp=0
      const= sqrt(nbrsp/(Lampp[b]*somme))
      for(j in 1 :nbrrel)
      {
        cosp= cosp + donpp[i,j]*Vpp[j, b]/sqrt(hcpp[j])*const
      }
      contrelspp[i,b]=cosp^2*100
      if (contrelspp[i,b] >=contrelsp[i,b]){
        prob3[i,b]=prob3[i,b]+1
      }
    }
  }
}
nombreper
prob1=prob1*100/nombreper
prob2=prob2*100/ nombreper
prob3=prob3*100/ nombreper
for(j in 1:nax)
{
  nomcol2[1+(j-1)*3]<-paste("coord",j)
  nomcol2[2+(j-1)*3] <-paste("CR",j)
  nomcol2[3+(j-1)*3] <-paste("P",j)
}
nomrang2[1]<-paste("Vp, Vp%, P")
for(j in 1:nax)
{
  resul2mfansca[1,3*(j-1)+1]=Lamp[j]
  resul2mfansca[1,3*(j-1)+2]=Lamp2[j]
  resul2mfansca[1,3*(j-1)+3]=prob1[j]
}
nomsp
for(i in 1 :nbrsp)
{
  nomrang2[i+1]=nomsp[i]
  for(j in 1:nax)
  {
    resul2mfansca[i+1,3*(j-1)+1]=coordsp[i,j]
    resul2mfansca[i+1,3*(j-1)+2]=contrelsp[i,j]
    resul2mfansca[i+1,3*(j-1)+3]=prob3[i,j]
  }
}
nomrel
for(i in 1 :nbrrel)
{
  nomrang2[i+1+nbrsp]=nomrel[i]
  for(j in 1:nax)
  {
    resul2mfansca[i+1+nbrsp,3*(j-1)+1]=coordrel[i,j]
    resul2mfansca[i+1+nbrsp,3*(j-1)+2]=contrelrel[i,j]
    resul2mfansca[i+1+nbrsp,3*(j-1)+3]=prob2[i,j]
  }
}
rownames(resul1mfansca)<-nomrang1
colnames(resul1mfansca)<-nomcol1
rownames(resul2mfansca)<-nomrang2
colnames(resul2mfansca)<-nomcol2



mfapca.R

don<-as.matrix(filename)
nbrsp= nrow(don)
nbrrel= ncol(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
nb=length(co)
nrang= nb+2+nbrsp+nbrrel
nbrcol=nax*2
nomrang<-list(1:nrang)
nomcol<-list(1:nbrcol)
for(j in 1:nax)
{
  nomcol[1+(j-1)*2]<-paste("coord",j)
  nomcol[2+(j-1)*2] <-paste("Cr%",j)
}
resulmfapca<-matrix(1:nbrcol*nrang, nrow=nrang,ncol=nbrcol,byrow=TRUE)
colnames(resulmfapca)<-nomcol
correl<-cor(t(don))
R<-eigen(correl,symmetric=TRUE)
V<- R$vectors; Lam <- R$values
Lam2=Lam*100/nbrsp
for(j in 1:nax)
{
  resulmfapca[1,2*(j-1)+1]=Lam[j]
  resulmfapca[1,2*(j-1)+2]=Lam2[j]
}
nomrang[1]<- paste("Vp & Vp% du tableau originel")
#Subtables and 1st eigenvalues for each subtable
vp1<-matrix(1 :nb)
i1=0
for(k in 1:nb)
{
  nbrspp= co[k]
  donst<-matrix(1:nbrspp*nbrrel,nrow=nbrspp,ncol=nbrrel)
  for(i in 1 :nbrspp)
  {
    i1=i1+1
    for(j in 1 :nbrrel)
    {
      donst[i,j]=don[i1,j]
    }
  }
  correlst<-cor(t(donst))
  Rst<-eigen(correlst,symmetric=TRUE)
  Lamst<-Rst$values
  Lamst2<-matrix(1:nax)
  vp1[k]=Lamst[1]
  Lamst2=Lamst*100/nbrspp
  for(j in 1:nax)
  {
    resulmfapca[1+k,2*(j-1)+1]=Lamst[j]
    resulmfapca[1+k,2*(j-1)+2]=Lamst2[j]
  }
  nomrang[1+k]<-paste("Vp & Vp% du sous-tableau", k)
}
#Complete weighted table by the 1rst eigenvalues of sub-tables
v1<-matrix(1 :nbrsp)
correlp<-matrix(1:nbrsp*nbrsp,nrow=nbrsp,ncol=nbrsp)
i1=0
for(k in 1:nb)
{
  nbrspp= co[k]
  for(i in 1 :nbrspp)
  {
    i1=i1+1
    v1[i1]=vp1[k]
  }
}
for(i in 1 :nbrsp)
{
  for(j in 1:nbrsp)
  {
    correlp[i,j]= correl[i,j]/sqrt(v1[i]*v1[j])
  }
}
#Distribution histogram of correlations with function hist()
nbrcel=nbrsp*(nbrsp-1)/2
matcor<-as.matrix(1 :nbrcel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 2:nbrsp)
{
  for (j2 in 2:j1)
  {
    i=i+1
    matcor[i]=correlp[j1,j2]
  }
}
hist(matcor, breaks=40, include.lowest = TRUE, right = TRUE,
     labels=FALSE, border = NULL,
     main = paste("Correl distribution-mfapca"), plot=TRUE,
     axes = TRUE,
     warn.unused = TRUE)
Rp<-eigen(correlp,symmetric=TRUE)
Vp<- Rp$vectors; Lamp <- Rp$values
trace =sum(Lamp)
trace
Lamp2<-matrix(1 :nbrsp)
Lamp2=Lamp*100/trace
for(j in 1:nax)
{
  resulmfapca[2+nb,2*(j-1)+1]=Lamp[j]
  resulmfapca[2+nb,2*(j-1)+2]=Lamp2[j]
}
nomrang[2+nb]<- paste("Vp & Vp% of the weighted table")
#Calculation of species coordinates and relative contributions on axes
coordsp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
for(i in 1:nbrsp)
{
  for(a in 1:nax)
  {
    coordsp[i,a] =Vp[i,a]*sqrt(Lamp[a]*v1[i])
    contrelsp[i,a] =coordsp[i,a]^2
  }
}

# Calculation of releve coordinates and relative contributions on axes
coordrel<-matrix(1:nbrrel*nax,nrow=nbrrel,ncol=nax)
contrelrel<-matrix(1:nbrrel*nax,nrow=nbrrel,ncol=nax)
nrmp<-matrix(rep(0,1*nbrrel),1,nbrrel)
#Calculation of weighted reduced centered variables
don1<-matrix(1:nbrsp*nbrrel, nrow = nbrsp,ncol=nbrrel)
sec<-matrix(1 :nbrsp)
moy<-rowMeans(don)
for(i in 1 :nbrsp)
{
  for(j in 1:nbrrel)
  {
    sec[i]=sec[i] +(don[i,j]-moy[i])^2
  }
  sec[i]=sqrt(sec[i]/nbrrel)
}
for(i in 1 :nbrsp)
{
  for(j in 1:nbrrel)
  {
    don1[i,j]=(don[i,j]-moy[i])/(sec[i]*sqrt(v1[i]))
  }
}
for(i in 1:nbrrel)
{
  for (j in 1:nbrsp)
  {
    nrmp[i] = nrmp[i]+(don1[j,i])^2
  }
}
for(i in 1:nbrrel)
{
  for(a in 1:nax)
  {
    coind=0
    for(j in 1:nbrsp)
    {
      coind= coind+ don1[j,i]*Vp[j,a]
    }
    coordrel[i,a]=coind
    contrelrel[i,a]= coordrel[i,a]* coordrel[i,a]/nrmp[i]
  }
}
for(i in 1 :nbrsp)
{
  nomrang[i+2+nb]<-nomsp[i]
  for(j in 1:nax)
  {
    resulmfapca[i+2+nb,2*(j-1)+1]=coordsp[i,j]
    resulmfapca[i+2+nb,2*(j-1)+2]=contrelsp[i,j]
  }
}
for(i in 1 :nbrrel)
{
  for(j in 1:nax)
  {
    resulmfapca[i+nbrsp+2+nb,2*(j-1)+1]=coordrel[i,j]
    resulmfapca[i+nbrsp+2+nb,2*(j-1)+2]=contrelrel[i,j]
  }
  nomrang[i+2+nb+nbrsp]<-nomrel[i]
}
rownames(resulmfapca)<-nomrang

nscaf.R

don<-as.matrix(filename)
nbrsp= nrow(don)
nbrrel= ncol(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
somme=sum(don)
nrang= 1+nbrsp+nbrrel
nbrcol=nax*3
nomrang<-list(1:nrang)
nomcol<-list(1:1+nax*3)
for(j in 1:nax)
{
  nomcol[1+(j-1)*3]<-paste("coord",j)
  nomcol[2+(j-1)*3] <-paste("Cr%",j)
  nomcol[3+(j-1)*3] <-paste("P",j)
}
resulnsca<-matrix(1:nrang*nbrcol, nrow=nrang,ncol=nbrcol,byrow=TRUE)
colnames(resulnsca)<- nomcol
#calculation of the transformed and covariance matrices
donc<- matrix(1 :nbrsp*nbrrel, nrow=nbrsp,ncol=nbrrel)
donc<-don
covar1<-matrix(rep(0,nbrrel*nbrrel),nbrrel,nbrrel)
covar1<-covarnsca(donc,nbrsp,nbrrel)
#distribution histogram of the covariance coefficients with the function hist()
nbrcel=nbrrel*(nbrrel-1)/2 +nbrrel
matcov<-as.matrix(1 :nbrcel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 1:nbrrel)
{
  for (j2 in 1:j1)
  {
    i=i+1
    matcov[i]=covar1[j1,j2]
  }
}
hist(matcov, breaks=40, include.lowest = TRUE, right = TRUE,
     labels=FALSE, border = NULL,
     main = paste("NSCA"), plot=TRUE,
     axes = TRUE,
     warn.unused = TRUE)
#calculation of the trace, of the eigenvectors and eigenvalues
Lam<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
Lam2<-matrix(1:nbrrel,nrow=1,ncol=nbrrel)
R<-eigen(covar1, symmetric=TRUE, only.values=FALSE)
V <- R$vectors; Lam <- R$values
trace = sum(Lam)
trace
Lam2<-Lam*100/trace
#calculation of the coordinates and relative contributions of the individuals on axes
coordrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
contrerel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
hc<-colSums(don)
for(j in 1:nbrrel)
{
  for(b in 1 :nax)
  {
    coordrel[j,b]=V[j,b]/sqrt(hc[j])*sqrt(somme)
    contrerel[j,b]=100* V[j,b]* V[j,b]
  }
}
coordrel
contrerel
#calculation of the coordinates and relative contributions of the species on axes
coordsp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
for(b in 1 :nax)
{
  for(i in 1 : nbrsp)
  {
    cosp=0
    const= sqrt(nbrsp/(Lam[b]*somme))
    for(j in 1 :nbrrel)
    {
      cosp= cosp+don[i,j]*V[j,b]/sqrt(hc[j])*const
    }
    coordsp[i,b]=cosp
    contrelsp[i,b]= cosp^2*100
  }
}
#permutations of the matrix and new calculations
prob1<-matrix(1:nax*1,nrow=1,ncol=nax,byrow=TRUE)
prob2<-matrix(1:nax*nbrrel,nrow=nbrrel,ncol=nax,byrow=TRUE)
prob3<-matrix(1:nax*nbrsp,nrow=nbrsp,ncol=nax,byrow=TRUE)
prob1[]<-0
prob2[]<-0
prob3[]<-0
nombreper<-0
for(np in 1:nper)
{
  Lamp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  Lamp2<-matrix(rep(0,1*nbrrel),1,nbrrel)
  donc<-matrix(rep(0,nbrsp*nbrrel),nbrsp,nbrrel)
  hc<-matrix(rep(0,1*nbrrel),1,nbrrel)
  covar1<-matrix(rep(0,nbrrel*nbrrel),nbrrel,nbrrel)
  for(i in 1 :nbrsp)
  {
    donc[i,]<-sample(don[i,],nbrrel,replace=FALSE)
  }
  hc<-colSums(donc)
  if(min(hc)==0){
    next
  }
  nombreper=nombreper+1
  covar1<-covarnsca(donc,nbrsp,nbrrel)
  #calculation  of the eigenvectors, of the eigenvalues and of the trace after permutations
  Rp<-eigen(covar1, symmetric=TRUE, only.values=FALSE)
  Vp <- Rp$vectors; Lamp <- Rp$values
  trace2=sum(Lamp)
  Lamp2<-Lamp*100/trace2
  for(b in 1 :nax)
  {
    if (Lamp[b]>=Lam[b])
    {
      prob1[b]=prob1[b]+1
    }
  }
  #calculation of the relative contributions of the individuals on axes after permutations
  contrerelp<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
  for(j in 1:nbrrel)
  {
    for(b in 1 :nax)
    {
      contrerelp[j,b]=100* Vp[j,b]* Vp[j,b]
      if(contrerelp[j,b] >=contrerel[j,b]){ 
        prob2[j,b]=prob2[j,b]+1
      }
    }
  }
  #calculation of the relative contributions of the species on axes after permutations
  contrelspp<-matrix(1 :nbrsp*nax,nrow=nbrsp,ncol=nax)
  for(b in 1 :nax)
  {
    for(i in 1 : nbrsp)
    {
      cosp=0
      const= sqrt(nbrsp/(Lamp[b]*somme))
      for(j in 1 :nbrrel)
      {
        cosp=cosp+donc[i,j]*Vp[j,b]/sqrt(hc[j])*const
      }
      contrelspp[i,b]=cosp^2*100
      if (contrelspp[i,b] >=contrelsp[i,b]){
        prob3[i,b]=prob3[i,b]+1
      }
    }
  }
}
nombreper
prob1=prob1*100/nombreper
prob2=prob2*100/ nombreper
prob3=prob3*100/ nombreper
for(j in 1:nax)
{
  resulnsca[1,3*(j-1)+1]=Lam[j]
  resulnsca[1,3*(j-1)+2]=Lam2[j]
  resulnsca[1,3*(j-1)+3]=prob1[j]
}
nomrang[1]<- paste("Vp,Vp% and P")
for(i in 1 :nbrsp)
{
  for(j in 1:nax)
  {
    resulnsca[i+1,3*(j-1)+1]=coordsp[i,j]
    resulnsca[i+1,3*(j-1)+2]=contrelsp[i,j]
    resulnsca[i+1,3*(j-1)+3]=prob3[i,j]
  }
  nomrang[i+1]=nomsp[i]
}
for(i in 1 :nbrrel)
{
  for(j in 1:nax)
  {
    resulnsca[i+1+nbrsp,3*(j-1)+1]=coordrel[i,j]
    resulnsca[i+1+nbrsp,3*(j-1)+2]=contrerel[i,j]
    resulnsca[i+1+nbrsp,3*(j-1)+3]=prob2[i,j]
  }
  nomrang[i+1+nbrsp]=nomrel[i]
}
rownames(resulnsca)<-nomrang

pca.R

don<-as.matrix(filename)
nbrsp= nrow(don)
nbrrel= ncol(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
nrang= 1+nbrsp+nbrrel
nbrcol=nax*3
nomrang<-list(1:nrang)
nomcol<-list(1:nax*3)
nomcol[1]<-paste("sp./releves")
for(j in 1:nax)
{
  nomcol[1+(j-1)*3]<-paste("coord",j)
  nomcol[2+(j-1)*3] <-paste("Cr%",j)
  nomcol[3+(j-1)*3] <-paste("P",j)
}
resulpca<-matrix(1:nrang*nbrcol, nrow=nrang,ncol=nbrcol,byrow=TRUE)
colnames(resulpca)<- nomcol
#calculation of the centered reduced matrix
don1<-matrix(1 :nbrsp*nbrrel, nrow=nbrsp, ncol=nbrrel)
sce<-matrix(rep(0,1*nbrsp),1,nbrsp)
moyenne<-rowMeans(don)
for(i in 1 :nbrsp)
{
  for(j in 1:nbrrel)
  {
    sce[i]=sce[i]+ (don[i,j]-moyenne[i])* (don[i,j]-moyenne[i])
  }
  sce[i]=sqrt(sce[i]/nbrrel)
}
for(i in 1 :nbrsp)
{
  for(j in 1:nbrrel)
  {
    don1[i,j]=(don[i,j]-moyenne[i])/sce[i]
  }
}
#calculation of the correlation matrix, of eigenvectors and of eigenvalues
correlation<-cor(t(don1))
R<-eigen(correlation, symmetric=TRUE)
V <- R$vectors; Lam <- R$values
trace=sum(Lam)
Lam2<-matrix(1 :nbrsp)
Lam2=Lam*100/trace
#correlation distribution histogram with function hist()
nbrcel=nbrsp*nbrsp
matcov<-as.matrix(1 :nbrcel, nrow=1,ncol=nbrcel)
i=0
for (j1 in 1:nbrsp)
{
  for (j2 in 1:nbrsp)
  {
    i=i+1
    matcov[i]<-correlation[j1,j2]
  }
}
hist(matcov, breaks=40, include.lowest = TRUE, right = TRUE,labels=FALSE, border = NULL,main = paste("PCA"), plot=TRUE,axes = TRUE,warn.unused = TRUE)
#calculation of coordinates and relatives contributions of species on axes
coordsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelsp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
contrelrel<-matrix(1 :nbrrel*nax,nrow=nbrrel,ncol=nax)
for(i in 1:nbrsp)
{
  for(j in 1:nax)
  {
    coordsp[i,j] =V[i,j]*sqrt(Lam[j])
    contrelsp[i,j] <- coordsp[i,j]^2
  }
}
# calculation of coordinates and relatives contributions of releves on axes
coordrel<-matrix(rep(0,nbrrel*nax),nbrrel,nax)
nrm<-matrix(rep(0,1*nbrrel),1,nbrrel)
for(i in 1:nbrrel)
{
  for (j in 1:nbrsp)
  {
    nrm[i] = nrm[i]+don1[j,i]^2
  }
}
for(i in 1:nbrrel)
{
  norme=nrm[i]
  for(a in 1:nax)
  {
    for(j in 1:nbrsp)
    {
      coordrel[i,a]= coordrel[i,a]+don1[j,i]*V[j,a]
    }
    contrelrel[i,a]= coordrel[i,a]* coordrel[i,a]/norme
  }
}
#permutations of the matrix and new calculations
prob1<-matrix(1:nax*1,nrow=1,ncol=nax,byrow=TRUE)
prob2<-matrix(1:nax*nbrsp,nrow=nbrsp,ncol=nax,byrow=TRUE)
prob3<-matrix(1:nax*nbrrel,nrow=nbrrel,ncol=nax,byrow=TRUE)
prob1[]<-0
prob2[]<-0
prob3[]<-0
Lamp<-matrix(1:nbrsp,nrow=1,ncol= nbrsp)
Lamp2<-matrix(1:nbrsp,nrow=1,ncol= nbrsp)
don1p<-matrix(1 :nbrsp*nbrrel, nbrsp,nbrrel)
coordspp<-matrix(1:nbrsp*nax,nrow=nbrsp,ncol=nax)
for(np in 1:nper)
{
  for(i in 1 :nbrsp)
  {
    don1p[i,]<-sample(don1[i,],nbrrel,replace=FALSE)
  }
  #calculation of the correlation matrices and of eigenvalues and eigenvectors after permutations
  correlationp = cor(t(don1p))
  Rp<-eigen(correlationp, symmetric=TRUE)
  Vp <- Rp$vectors; Lamp <- Rp$values
  Lamp2<-Lamp*100/nbrsp
  for(a in 1 :nax)
  {
    if (Lamp[a]>=Lam[a]){
      prob1[a]=prob1[a]+1
    }
  }
  #calculation of the coordinates and relative contributions of species on axes after permutations
  for(i in 1:nbrsp)
  {
    for(a in 1:nax)
    {
      coordspp[i,a] =Vp[i,a]*sqrt(Lamp[a])
    }
  }
  # calculation of releve coordinates on axes after permutations
  for(i in 1:nbrsp)
  {
    for(a in 1:nax)
    {
      contrelspp =coordspp[i,a]*coordspp[i,a]
      if (contrelspp>=contrelsp[i,a]){ 
        prob2[i,a]=prob2[i,a]+1
      }
    }
  }
  # calculation of the coordinates and relative contributions of species on axes after permutations
  nrmp<-matrix(rep(0,1*nbrrel),1,nbrrel)
  for(i in 1:nbrrel)
  {
    for (j in 1:nbrsp)
    {
      nrmp[i] = nrmp[i]+(don1p[j,i]* don1p[j,i])
    }
  }
  for(i in 1:nbrrel)
  {
    for(a in 1:nax)
    {
      coind=0
      for(j in 1:nbrsp)
      {
        coind= coind+don1p[j,i]*Vp[j,a]
      }
      contrelrelp= coind* coind/nrmp[i]
      if (contrelrelp>=contrelrel[i,a]){
        prob3[i,a]=prob3[i,a]+1
      }
    }
  }
}
prob1=prob1*100/nper
prob2=prob2*100/ nper
prob3=prob3*100/ nper
for(j in 1:nax)
{
  resulpca[1,3*(j-1)+1]=Lam[j]
  resulpca[1,3*(j-1)+2]=Lam2[j]
  resulpca[1,3*(j-1)+3]=prob1[j]
}
nomrang[1]<- paste("Vp,Vp% and P")
for(i in 1 :nbrsp)
{
  for(j in 1:nax)
  {
    resulpca[i+1,3*(j-1)+1]=coordsp[i,j]
    resulpca[i+1,3*(j-1)+2]=contrelsp[i,j]
    resulpca[i+1,3*(j-1)+3]=prob2[i,j]
  }
  nomrang[i+1]=nomsp[i]
}
for(i in 1 :nbrrel)
{
  for(j in 1:nax)
  {
    resulpca[i+1+nbrsp,3*(j-1)+1]=coordrel[i,j]
    resulpca[i+1+nbrsp,3*(j-1)+2]=contrelrel[i,j]
    resulpca[i+1+nbrsp,3*(j-1)+3]=prob3[i,j]
  }
  nomrang[i+1+nbrsp]=nomrel[i]
}
rownames(resulpca)<-nomrang

phyto1.R

don<-as.matrix(filename)
nbrrel= ncol(don)
nbrsp= nrow(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
res<-matrix(1 :nbrrel*nbrsp,nrow=nbrsp,ncol=nbrrel,byrow=TRUE)
for(i in 1:nbrsp)
{
  for(j in 1:nbrrel)
  {
    if (don[i,j]!=".")
    {
      res[i,j]=1
    }
    else 
    {
      res[i,j]=0
    }
  }
}
rownames(res)<-c(nomsp)
colnames(res)<-c(nomrel)
resul<-data.frame(res)

phyto12345.R

don<-as.matrix(filename)
nbrrel= ncol(don)
nbrsp= nrow(don)
nomsp<-rownames(don)
nomrel<-colnames(don)
#calcul nombre de lignes du nouveau tableau
#definition of a new table
nomsp2<-list(1 :nbrsp)
don1<-matrix(1 :nbrrel*nbrsp,nrow=nbrsp,ncol=nbrrel,byrow=TRUE)
for(i in 1:nbrsp)
{
  nomsp2[i]<-paste(nomsp[i])
}

for(i in 1:nbrsp)
{
  for(j in 1:nbrrel)
  {
    if (don[i,j] ==".")
    {
      don1[i,j]=0
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j] =="+"| don[i,j] =="r"| don[i,j] =="i"| don[i,j] ==1)
    {
      don1[i,j]=1
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]==1)
    {
      don1[i,j]=1
    }
  }
  for(j in 1:nbrrel)
  {
    if (don[i,j]>1)
    {
      don1[i,j]<-paste(don[i,j])
    }
  }
}
rownames(don1)<-c(nomsp2)
colnames(don1)<-c(nomrel)
resul<-data.frame(don1)

