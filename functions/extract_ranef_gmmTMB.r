#' @title random effects data frame for glmmTMB objects
#' @description extracts random effects from glmmTMB objects into a data frame
#' @param m object resulting from model fitted using \code{glmmTMB\{glmmTMB\}}.
#' @return a data frame
#' @details
#' the function largely serves to make it convenient to get the results for the random effects nicely formatted.
#' @seealso
#' \code{\link{glmmTMB}}
#' @examples
#' \dontrun{
#' library(glmmTMB)
#' m <- glmmTMB(count ~ spp + mined + (1|site), zi=~spp + mined + (1|site), family=nbinom2, data=Salamanders)
#' extract.ranef.glmmTMB(m=m)
#' }
#' @export
extract.ranef.glmmTMB<-function(m){
  to.do=summary(m)$varcor
  xx=lapply(to.do, function(x){
    lapply(x, attr, "stddev")
  })
  res=suppressWarnings(data.frame(
    part=rep(names(to.do), unlist(lapply(xx, function(x){sum(unlist(lapply(x, length)))}))),
    grp=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){unlist(lapply(x, length))}))),
    var1=unlist(lapply(xx, function(x){unlist(lapply(x, names))})),
    var2=NA,
    sdcor=unlist(xx)
  ))
  xx=lapply(to.do, function(x){
    lapply(x, attr, "correlation")
  })
  xx=xx[unlist(lapply(xx, length))>0]
  xx=suppressWarnings(data.frame(
    part=rep(names(xx), unlist(sum(unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))))),
    grp=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))),
    var1=unlist(lapply(xx, function(x){lapply(x, function(x){rep(rownames(x), times=ncol(x))})})),
    var2=unlist(lapply(xx, function(x){lapply(x, function(x){rep(colnames(x), each=nrow(x))})})),
    sdcor=unlist(xx)
  ))
  xx=xx[as.character(xx[, "var1"])<as.character(xx[, "var2"]), ]
  res=rbind(res, xx)
  res=res[order(as.character(res$part), as.character(res$grp), as.character(res$var1)), ]
  rownames(res)=NULL
  return(res)
}
