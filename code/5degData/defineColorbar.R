defineColorbar <- function(z,rr=NULL,set_white=NULL,col_scheme="bwr")
{
  if(is.null(rr)) rr = range(z, na.rm = TRUE)
  brk = seq(rr[1],rr[2], length = 500)
  brk.ind = round(seq(1, length(brk), length = 10))
  brk.lab = round(brk[brk.ind], 2)
  brk.at = brk[brk.ind]
  if (is.null(set_white)) set_white = 0.0
  zero.ind = min(which(brk > set_white))/length(brk)
  color <- designer.colors(n = length(brk) - 1, 
                           col = c("darkblue","white", "darkred"), 
                           x = c(0, zero.ind, 1))

  return(list(color=color,brk=brk,brk.at=brk.at,brk.lab=brk.lab))

}