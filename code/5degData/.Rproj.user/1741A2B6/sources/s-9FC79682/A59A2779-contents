library(grid)
plotSeaMask <- function(x.lim,y.lim,l.force,plotornot=FALSE) {

  outline <- map("world",fill=TRUE,col='gray70',xlim=x.lim,ylim=y.lim,lforce=l.force,plot=plotornot)
  xrange <- range(outline$x, na.rm=TRUE) # get bounding box
  yrange <- range(outline$y, na.rm=TRUE)
  xbox <- xrange
  ybox <- yrange
  polypath(c(outline$x, NA, c(xbox, rev(xbox))),
           c(outline$y, NA, rep(ybox, each=2)),
           col="light blue", rule="evenodd")
}
