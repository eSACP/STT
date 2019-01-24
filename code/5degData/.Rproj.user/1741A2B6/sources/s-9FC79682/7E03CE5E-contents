## Preliminaries
rm(list=ls())
library(INLA)
library(maps)
library(fields)
library(excursions)
library(RColorBrewer)
source("plotSeaMask.R")
source("defineColorbar.R")

set.seed(101)
objDir <- "./robj/"
season <- "JJA"
make.plot <- TRUE



##########################################################################################################
##
## M O D E L   F I T T I N G   V I A   I N L A
##
##########################################################################################################

##########################################################################################################
## Load gridded temperatures data object data.inla
##   data.inla is a data frame with S * n.year rows:
##        S = 70 measurement locations
##        n.year = 65 annual measurements
##   and 6 columns:
##        t.vec:      Time (regularized to [-3.25,3.25] to reflect a decadal time unit)
##        y.vec:      Temperature anomalies (in degrees Celsius)
##        lon.vec:    Longitudes of measurement locations (repeatedly)
##        lat.vec:    Latitudes of measurement locations (repeatedly)
##        station.id: Indices for the different measurement locations
##        time:       Indices for the different time points
load(file=paste0(objDir,season,"_data.inla.5deg.RData"))
S <- length(unique(data.inla$station.id))
n.year <- length(unique(data.inla$time))
y.stations <- matrix(data.inla$y.vec,byrow=T,nrow=n.year)


##########################################################################################################
## INLA mesh
lon.vec <- data.inla$lon.vec[1:S]
lat.vec <- data.inla$lat.vec[1:S]
coords <- cbind(lon.vec,lat.vec)
mesh <- inla.mesh.2d(loc=coords,
                     offset=c(7.5,15),
                     max.edge=c(10,10),
                     min.angle=c(21,21))


##########################################################################################################
## SPDE model with Matérn covariance structure
spde <- inla.spde2.matern(mesh, alpha=2, theta.prior.prec=1.5, constr=FALSE)


##########################################################################################################
## Projection from SPDE mesh to an output 1deg lon/lat lattice in measurement space
lon.prd <- c(-12,44)
lat.prd <- c(34,72)
stepsize <- 1
nxy <- round(c(diff(range(lon.prd)),diff(range(lat.prd)))/stepsize)
lattice.prd <- inla.mesh.projector(mesh,xlim=range(lon.prd),ylim=range(lat.prd), dims=nxy+1)
lattice.prd <- submesh.grid(matrix(1,nxy[1]+1,nxy[2]+1),list(loc=lattice.prd$lattice$loc,dims=nxy+1))


##########################################################################################################
## Projector matrices of the two GMRFs of the SPDE model, and link between INLA mesh and regular lattice
##    tau field: AR(1) temporal correlation with no covariance in space
##    xi field:  Spatial covariance structure a la Matérn, no temporal correlation
A.proj.tau <- inla.spde.make.A(mesh,
                               loc=as.matrix(cbind(rep(coords[,1],n.year),rep(coords[,2],n.year))),
                               group=data.inla$time,
                               n.group=n.year)
                          
## Dynamic regression field; purely spatial coefficients of a trend component (additional to a linear term)
A.proj.alpha1s <- inla.spde.make.A(mesh,
                                   loc=as.matrix(cbind(rep(lon.vec,n.year),rep(lat.vec,n.year))),
                                   weights=data.inla$t.vec)   # NB! t.vec scaled to be in [-3.25,3.25]

## Projector matrix from irregular INLA mesh to regular lattice in observation space
A.prd.lattice <- inla.spde.make.A(mesh,loc=lattice.prd$loc)


##########################################################################################################
## Spatial indices of the two field components of the SPDE model
ind.tau <- inla.spde.make.index(name="field.tau",n.spde=spde$n.spde,n.group=n.year)
ind.alpha1s <- inla.spde.make.index(name="field.alpha1s",n.spde=spde$n.spde)


##########################################################################################################
## Problem formulation using a 'stack' data structure in INLA
## For construction details, see eg 'stack.est' on page 23 of https://inla.googlecode.com/hg/r-inla.org/case-studies/Blangiardo-et-al-2012/Report-version-Oct-2012.pdf

## Estimation stack for fitting the model
stack.est.D <- inla.stack(data=list(y=data.inla$y.vec),
                          A=list(A.proj.tau,A.proj.alpha1s,1),
                          effects=list(ind.tau,ind.alpha1s,list(trend=data.inla$t.vec)),
                          tag="est.D")

## Prediction stack for alpha.1 + alpha.1s (global fixed effect + spatial coefficient field)
#  Prediction onto mesh in INLA space, mesh$n locations
stack.prd.mesh.D <- inla.stack(data=list(y=NA),
                               A=list(1,1),
                               effects=list(ind.alpha1s,list(trend=rep(1,mesh$n))),
                               tag="prd.mesh.D")
#  Prediction onto lattice in measurement space, 57 x 39 1deg grid cells
stack.prd.lattice.D <- inla.stack(data=list(y=NA),
                                  A=list(A.prd.lattice,1),
                                  effects=list(ind.alpha1s,list(trend=rep(1,lattice.prd$n))),
                                  tag="prd.lattice.D")

## Full stack for estimation and prediction
stack.D <- inla.stack(stack.est.D,stack.prd.mesh.D,stack.prd.lattice.D)


##########################################################################################################
## Fit model via INLA
##

## INLA formula, model D
##    g_s(t) = (alpha_1 + alpha_1s) * t  +  tau_ts   (alpha_1s Matern in space, no temporal correlation)
##             tau_ts = a * tau_(t-1)s + xi_ts       (tau_ts AR(1) in time)
##             xi_ts ~ N(0,Sigma_xi^2)               (xi_ts Matern in space)
formula.D <- y ~ -1 + trend + f(field.alpha1s,model=spde) + 
                              f(field.tau,model=spde,group=field.tau.group,control.group=list(model="ar1"))

## Run INLA
result.D <- inla(formula.D,
                 data=inla.stack.data(stack.D,spde=spde),
                 family="gaussian",
                 verbose=TRUE,
                 control.family=list(hyper=list(prec=list(initial=10,fixed=TRUE))),
                 control.compute=list(config=TRUE),
                 control.predictor=list(A=inla.stack.A(stack.D),compute=TRUE),
                 control.inla=list(int.strategy="eb"))



##########################################################################################################
##
## T R E N D   S I G N I F I C A N C E   A S S E S S M E N T
##
##########################################################################################################

##########################################################################################################
## Marginal credible intervals
##########################################################################################################

## Marginal posteriors of trend field (alpha1+alpha1.s) on lattice
index.lattice.prd <- inla.stack.index(stack.D,tag="prd.lattice.D")$data

## Level 1% and 5%
marg.post.005quant.lattice <- NULL
marg.post.025quant.lattice <- NULL
marg.post.975quant.lattice <- NULL
marg.post.995quant.lattice <- NULL
for (j in 1:length(index.lattice.prd)) {
  marg.post.005quant.lattice <- c(marg.post.005quant.lattice,inla.qmarginal(0.005,result.D$marginals.linear.predictor[[index.lattice.prd[j]]]))
  marg.post.025quant.lattice <- c(marg.post.025quant.lattice,inla.qmarginal(0.025,result.D$marginals.linear.predictor[[index.lattice.prd[j]]]))
  marg.post.975quant.lattice <- c(marg.post.975quant.lattice,inla.qmarginal(0.975,result.D$marginals.linear.predictor[[index.lattice.prd[j]]]))
  marg.post.995quant.lattice <- c(marg.post.995quant.lattice,inla.qmarginal(0.995,result.D$marginals.linear.predictor[[index.lattice.prd[j]]]))
}
marg.significant005 <- as.integer(!(marg.post.025quant.lattice<0 & marg.post.975quant.lattice>0))
marg.significant001 <- as.integer(!(marg.post.005quant.lattice<0 & marg.post.995quant.lattice>0))


##########################################################################################################
## Avoidance excursion set calculations for level u=0
##########################################################################################################
u.curr <- 0

# INLA space (mesh)
exc.mesh <- excursions.inla(result.D,stack.D,tag="prd.mesh.D",u=u.curr,u.link=TRUE,type="!=",F.lim=0.6,method='QC')

# Observation space (lattice.prd)
exc.lattice <- excursions.inla(result.D,stack.D,tag="prd.lattice.D",u=u.curr,u.link=TRUE,type="!=",F.lim=0.6,method='QC')



##########################################################################################################
##
## P L O T T I N G
##
##########################################################################################################
if (make.plot) {

  ########################################################################################################
  ## Preliminaries
  ##   -Plotting area
  ##   -Image matrices
  ##   -Legends and colour palettes
  
  ## Fixed parameters
  num.res <- 5
  resolution <- paste0(num.res,"deg")
  plt.axes <- TRUE
  probs <- c(0.05,0.01)
  lon.prd <- c(-12,41.6)
  lat.prd <- c(35.2,72)
  lon.plt <- c(-18,48)
  lat.plt <- c(25,80)
  
  ## Harmonize mesh and lattice extent for plotting
  idx.lat <- (lat.prd[1]<=mesh$loc[,2]&mesh$loc[,2]<lat.prd[2])
  idx.lon <- (lon.prd[1]<=mesh$loc[,1]&mesh$loc[,1]<lon.prd[2])
  idx.mesh <- idx.lat*idx.lon
  
  ## Number of lon/lat grid cells
  n.lon.data <- ceiling(diff(range(mesh$loc[idx.mesh==TRUE,1]))/num.res)
  n.lat.data <- ceiling(diff(range(mesh$loc[idx.mesh==TRUE,2]))/num.res)
  n.lon.lattice <- ceiling(diff(range(mesh$loc[idx.mesh==TRUE,1]))/1)
  n.lat.lattice <- ceiling(diff(range(mesh$loc[idx.mesh==TRUE,2]))/1)
  
  ## Which mesh points correspond to data points
  id.plot <- mesh$idx$loc
  lon.offset <- 2
  lat.offset <- 2
  lon.max <- max(mesh$loc[id.plot,1])+lon.offset
  lon.min <- min(mesh$loc[id.plot,1])-lon.offset
  lat.max <- max(mesh$loc[id.plot,2])+lat.offset
  lat.min <- min(mesh$loc[id.plot,2])-lat.offset
  ixy.mesh <- which(mesh$loc[,1]>lon.min & 
                      mesh$loc[,1]<lon.max &
                      mesh$loc[,2]>lat.min & 
                      mesh$loc[,2]<lat.max)
  ixy.lattice <- which(lattice.prd$loc[,1]>lon.min & 
                         lattice.prd$loc[,1]<lon.max &
                         lattice.prd$loc[,2]>lat.min & 
                         lattice.prd$loc[,2]<lat.max)
  
  ## Image matrices
  postmean.mesh <- as.image(exc.mesh$mean[ixy.mesh],x=mesh$loc[ixy.mesh,1:2],nx=n.lon.data,ny=n.lat.data)   #counting the number of mesh locations in either direction
  poststdev.mesh <- as.image(sqrt(exc.mesh$vars[ixy.mesh]),x=mesh$loc[ixy.mesh,1:2],nx=n.lon.data,ny=n.lat.data)   #counting the number of mesh locations in either direction
  postmean.lattice <- as.image(exc.lattice$mean[ixy.lattice],x=lattice.prd$loc[ixy.lattice,1:2],nx=n.lon.lattice,ny=n.lat.lattice) 
  poststdev.lattice <- as.image(sqrt(exc.lattice$vars[ixy.lattice]),x=lattice.prd$loc[ixy.lattice,1:2],nx=n.lon.lattice,ny=n.lat.lattice) 
  margsignif005.lattice <- as.image(marg.significant005[ixy.lattice],x=lattice.prd$loc[ixy.lattice,1:2],nx=n.lon.lattice,ny=n.lat.lattice) 
  margsignif005.lattice$z[margsignif005.lattice$z>0] <- 1
  margsignif001.lattice <- as.image(marg.significant001[ixy.lattice],x=lattice.prd$loc[ixy.lattice,1:2],nx=n.lon.lattice,ny=n.lat.lattice) 
  margsignif001.lattice$z[margsignif001.lattice$z>0] <- 1
  
  ## Load legend range values
  load(file=paste0(objDir,"zzRange.RData"))
  rr <- c(-max(abs(zz)),max(abs(zz)))
  
  ## Trend posterior mean palettes
  cmap <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
  cmap.F <- colorRampPalette(brewer.pal(9, "Greens"))(100)
  cmap.RdBu <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
  cmap.rwb <- defineColorbar(z=zz,rr=rr,set_white=0,col_scheme="bwr")
  cmap.rwb.mean <- cmap.rwb$color
  brk <- cmap.rwb$brk
  brk.at <- cmap.rwb$brk.at
  brk.lab <- cmap.rwb$brk.lab
  
  ColorLevels <- seq(rr[1],rr[2],length=length(cmap.rwb.mean))
  ColorLevels.std <- seq(zz.std[1],zz.std[2],length=length(cmap.F))
  
  
  ########################################################################################################
  ## Generate various plots included in the paper based on the 5deg E-OBS data set
  ########################################################################################################
  
  ########################################################################################################
  ## Plot mesh
  ########################################################################################################
  par(mfrow=c(1,1),oma=c(0,0,0,0))
  plot(mesh$loc[,1:2], col="gray30",xlab="lon",ylab="lat")
  map("world",col="gray45",fill=TRUE,add=TRUE)
  plot(mesh, col="gray30", add=TRUE)
  points(coords, col=2, pch=16, cex=0.6)
  points(mesh$loc[mesh$idx$loc,1:2], col=3, pch=16, cex=0.95)
  
  
  ########################################################################################################
  ## Posterior trend distribution (mean and stdev) for mesh and lattice
  ########################################################################################################
  
  par(mfrow=c(1,1),oma=c(0,0,1,5),cex.main=1)
  
  ########################################################################################################
  ## Mean
  
  ## INLA space (mesh)
  ## Grid scale in degrees and choose kernel bandwidth to be 1.5 degrees. 
  look <- image.smooth(postmean.mesh,theta=1.5)
  to.plot <- look$z
  col.plot.min <- which.min(abs(ColorLevels-range(to.plot,na.rm=T)[1]))
  col.plot.max <- which.min(abs(ColorLevels-range(to.plot,na.rm=T)[2]))
  par(oma=c(0,0,0,7))
  image(look,col=cmap.rwb.mean[col.plot.min:col.plot.max],axes=plt.axes,xlim=lon.prd,ylim=lat.prd) 
  map("world",xlim=lon.prd, ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
  plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
  par(oma=c(0,0,1,1))
  image.plot(legend.only=TRUE,breaks=brk,at=brk.at,label=brk.lab,col=cmap.rwb.mean,legend.width=.8,legend.shrink=.90,
             axis.args=list(cex.axis=0.7),legend.args=list(text="Trend (degC/decade)",side=4,font=2,line=2,cex=.9)) 
  
  ## Observation space (lattice.prd)
  ## Grid scale in degrees and choose kernel bandwidth to be .5 degrees. 
  look <- image.smooth(postmean.lattice,theta=.5)
  to.plot <- look$z
  col.plot.min <- which.min(abs(ColorLevels-range(to.plot,na.rm=T)[1]))
  col.plot.max <- which.min(abs(ColorLevels-range(to.plot,na.rm=T)[2]))
  par(oma=c(0,0,0,7))
  image(look,col=cmap.rwb.mean[col.plot.min:col.plot.max],axes=plt.axes,xlim=lon.prd,ylim=lat.prd) 
  map("world",xlim=lon.prd, ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
  plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
  par(oma=c(0,0,1,1))
  image.plot(legend.only=TRUE,breaks=brk,at=brk.at,label=brk.lab,col=cmap.rwb.mean,legend.width=.8,legend.shrink=.90,
             axis.args=list(cex.axis=0.7),legend.args=list(text="Trend (degC/decade)",side=4,font=2,line=2,cex=.9)) 
  
  ########################################################################################################
  ## StDev
  
  ## INLA space (mesh)
  ## Grid scale in degrees and choose kernel bandwidth to be 1.5 degrees. 
  look <- image.smooth(poststdev.mesh,theta=1.5)
  to.plot <- look$z
  col.plot.min <- which.min(abs(ColorLevels.std-range(to.plot,na.rm=T)[1]))
  col.plot.max <- which.min(abs(ColorLevels.std-range(to.plot,na.rm=T)[2]))
  par(oma=c(0,0,0,7))
  image(look,col=cmap.F[col.plot.min:col.plot.max],axes=plt.axes,xlim=lon.prd,ylim=lat.prd) 
  map("world",xlim=lon.prd, ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
  plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
  par(oma=c(0,0,1,1))
  image.plot(legend.only=TRUE,zlim=zz.std,col=cmap.F,legend.width=.8,legend.shrink=.90,
             axis.args=list(cex.axis=0.7),legend.args=list(text="Stdev (degC/decade)",side=4,font=2,line=2,cex=.9)) 
  
  ## Observation space (lattice.prd)
  ## Grid scale in degrees and choose kernel bandwidth to be .5 degrees. 
  look <- image.smooth(poststdev.lattice,theta=.5)
  to.plot <- look$z
  col.plot.min <- which.min(abs(ColorLevels.std-range(to.plot,na.rm=T)[1]))
  col.plot.max <- which.min(abs(ColorLevels.std-range(to.plot,na.rm=T)[2]))
  par(oma=c(0,0,0,7))
  image(look,col=cmap.F[col.plot.min:col.plot.max],axes=plt.axes,xlim=lon.prd,ylim=lat.prd) 
  map("world",xlim=lon.prd, ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
  plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
  par(oma=c(0,0,1,1))
  image.plot(legend.only=TRUE,zlim=zz.std,col=cmap.F,legend.width=.8,legend.shrink=.90,
             axis.args=list(cex.axis=0.7),legend.args=list(text="Stdev (degC/decade)",side=4,font=2,line=2,cex=.9)) 
  
  
  ########################################################################################################
  ## Marginal posterior significances in observation space (lattice.prd)
  ########################################################################################################
  
  ## Level 5% 
  par(oma=c(0,0,0,7))
  image(margsignif005.lattice,col=c("white","red"),axes=plt.axes,xlim=lon.prd,ylim=lat.prd) 
  map("world",xlim=lon.prd, ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
  plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
  par(oma=c(0,0,1,1))
  mtext(text=paste0("1deg lattice (",resolution,"): Posterior marginal trend significances 5% (alpha.1 + alpha.1s)"),side=3,outer=TRUE,cex=1)
  
  ## Level 1% 
  par(oma=c(0,0,0,7))
  image(margsignif001.lattice,col=c("white","red"),axes=plt.axes,xlim=lon.prd,ylim=lat.prd) 
  map("world",xlim=lon.prd, ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
  plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
  par(oma=c(0,0,1,1))
  mtext(text=paste0("1deg lattice (",resolution,"): Posterior marginal trend significances 1% (alpha.1 + alpha.1s)"),side=3,outer=TRUE,cex=1)
  
  
  ########################################################################################################
  ## Avoidance excursion sets
  ########################################################################################################
  
  ########################################################################################################
  ## Loop over different error probabilities
  for (err.prob in probs) {
    
    alpha.str <- paste0(substr(err.prob,1,1),substr(err.prob,3,4))
    
    ## Excursion sets: Level u avoidance sets (only code for lattice included)
    sets.mesh <- continuous(exc.mesh,mesh,alpha=err.prob)
    sets.lattice <- continuous(exc.lattice,lattice.prd,alpha=err.prob)
    
    ## lattice
    par(oma=c(0,0,0,1))
    plot(sets.lattice$M["1"], col = "red", main=paste0("1deg lattice (",resolution,"): Avoidance excursion set, u=",u.curr),xpd=FALSE,xlim=lon.prd,ylim=lat.prd,axes=plt.axes,xaxs="i",yaxs="i")
    map("world",xlim=lon.prd,  ylim=lat.prd,col="gray45",fill=FALSE,add=TRUE)
    plotSeaMask(x.lim=lon.plt,y.lim=lat.plt,l.force='e',plotornot=FALSE)
    
  }
  
}

