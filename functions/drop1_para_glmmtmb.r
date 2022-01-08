drop1p<-function(model.res, para=F, data=NULL, contr=NULL, n.cores=c("all-1", "all"), to.del=NULL, return.model.results=F){
	#written Nov. 2018
	##last changed Dec. 2019 (added dropping from the zi- and dispformula; but not thoroughly tested yet)
	if(is.null(contr)){contr=glmmTMBControl()}
	##determine terms that can be dropped:
	xx=as.character(model.res$call)
	names(xx)=names(model.res$call)
  model.fe.re=xx["formula"]
  model.zi=as.formula(xx["ziformula"])
	model.disp=as.formula(xx["dispformula"])
	xfamily=xx["family"]
	if(grepl(xfamily, pattern="(", fixed=T)){
		xfamily=unlist(strsplit(xfamily, split="(", fixed=T))
		xfamily[2]=gsub(x=xfamily[2], pattern=")", replacement="", fixed=T)
		xfamily[2]=gsub(x=xfamily[2], pattern="link = ", replacement="", fixed=T)
		xfamily[2]=gsub(x=xfamily[2], pattern="\"", replacement="", fixed=T)
		if(xfamily[2]=="log"){
			xfamily=get(xfamily[1])(make.link(xfamily[2]))
		}else{
			xfamily=get(xfamily[1])(xfamily[2])
		}
	}
	##need to address how weights are recognized
	if(any(names(xx)=="weights")){
		data$XXXweights=data[, xx["weights"]]
	}else{
		data$XXXweights=1
	}
	
	to.drop.fun<-function(mf){
		#xcall=gsub(x=xcall, pattern="\n", replacement=T)
		model.terms=attr(terms(as.formula(mf)), "term.labels")
		model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
		model.terms=strsplit(model.terms, split=":", fixed=T)
		to.del=unlist(lapply(1:length(model.terms), function(m1){
			xx=unlist(lapply((1:length(model.terms))[-m1], function(m2){
				length(intersect(model.terms[[m1]], model.terms[[m2]]))==length(model.terms[[m1]])
			}))
			sum(xx)>0
		}))
		model.terms=model.terms[!to.del]
		model.terms=unlist(lapply(model.terms, paste, collapse=":"))
		return(model.terms)
	}
	xcall=as.character(model.fe.re)
	if(is.null(to.del)){
		model.terms=to.drop.fun(mf=model.fe.re)
	}else{
		model.terms=to.del
	}
	drop.from.zi=NULL
	drop.from.disp=NULL
	if(!all(as.character(model.zi)%in%c("~", "0", "1"))){
		drop.from.zi=to.drop.fun(mf=as.character(model.zi))
	}
	if(!all(as.character(model.disp)%in%c("~", "0", "1"))){
		drop.from.disp=to.drop.fun(mf=as.character(model.disp))
	}
	
	#if data are handed over
	if(length(data)>0){
		#check whether all columns needed are in them
		#figure out names of the variables involved:as.character(model.zi)
		model.terms2=c(attr(terms(as.formula(xcall)), "term.labels"), attr(terms(as.formula(as.character(model.zi))), "term.labels"),
			attr(terms(as.formula(as.character(model.disp))), "term.labels"))
		model.terms2=gsub(x=model.terms2, pattern="||", replacement="|", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="|", replacement="+", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern=":", replacement="+", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="*", replacement="+", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern=" ", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="I(", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="^2)", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="(", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern=")", replacement="", fixed=T)
		model.terms2=unique(unlist(strsplit(model.terms2, split="+",fixed=T)))
		model.terms2=model.terms2[!model.terms2%in%c("0", "1")]
		#deal with cbind response:
		#guess not needed anymore
# 		if(substr(x=as.character(xcall), start=1, stop=6)=="cbind("){
# 			xx=unlist(strsplit(as.character(xcall), split="~", fixed=T))[1]
# 			xx=gsub(xx, pattern="cbind(", replacement="", fixed=T)
# 			xx=gsub(xx, pattern=")", replacement="", fixed=T)
# 			xx=gsub(xx, pattern=" ", replacement="", fixed=T)
# 			model.terms2=c(model.terms2, unlist(strsplit(xx, split=",", fixed=T)))
# 		}
		#deal with sine/cosine included into model:
		#guess also not needed anymore:
# 		model.terms2=model.terms2[!substr(model.terms2, start=1, stop=4)=="sin("]
# 		model.terms2=model.terms2[!substr(model.terms2, start=1, stop=4)=="cos("]
		
		#still missing: offset term... (weights and squared terms should be fine)
		if(any(!model.terms2%in%names(data))){
			xx=model.terms2[!model.terms2%in%names(data)]
			stop(paste(c(paste(xx, collapse=", "), ifelse(length(xx)==1, "is ", "are "), "missing in data"), collape=""))
		}
	}else{
		stop("Error: no data frame handed over to argument 'data'")
	}
	model.terms2=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
	##prepare data:
	#model=paste(model.terms, collapse="+")##get model wrt fixed effects
	ii.data=data
	#all.models=c(xcall, paste(xcall, model.terms, sep="-"))
	all.models=cbind(formula=paste(as.character(model.fe.re), model.terms, sep="-"), 
		ziformula=paste(as.character(model.zi), collapse=""), dispformula=paste(as.character(model.disp), collapse=""))
	if(!is.null(drop.from.zi)){
		all.models=rbind(all.models, cbind(formula=as.character(model.fe.re), 
			ziformula=paste(paste(as.character(model.zi), collapse=""), drop.from.zi, sep="-"), 
			dispformula=paste(as.character(model.disp), collapse="")))
	}
	if(!is.null(drop.from.disp)){
		all.models=rbind(all.models, cbind(formula=as.character(model.fe.re), 
			ziformula=paste(as.character(model.zi), collapse=""), 
			dispformula=paste(paste(as.character(model.disp), collapse=""), drop.from.disp, sep="-")))
	}
	all.models=lapply(1:nrow(all.models), function(x){
		return(c(formula=as.vector(all.models[x, "formula"]), ziformula=as.vector(all.models[x, "ziformula"]), 
			dispformula=as.vector(all.models[x, "dispformula"])))
	})
	names(all.models)=rep(c("cond", "zi", "disp"), times=c(length(model.terms), length(drop.from.zi), length(drop.from.disp)))
	##define function doing the model comparison:
	model.fun<-function(i.models, idata, contr., xfamily, xweights){
		#i.model=as.formula(i.models)
		i.res=try(glmmTMB(formula=as.formula(i.models["formula"]), ziformula=as.formula(i.models["ziformula"]), 
			dispformula=as.formula(i.models["dispformula"]), data=idata, family=xfamily, weights=XXXweights, 
			control=contr), silent=T)
		return(i.res)
	}
	##run all models:
  if(para){
    require(parallel)
    n.cores=n.cores[1]
    cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores=="all"){
			n.cores=length(cl)
		}else if(n.cores=="all-1"){
			if(n.cores=="all-1"){n.cores=length(cl)-1}#else{n.cores=length(cl)}
		}
		n.cores=min(c(n.cores, length(all.models)))
		cl=cl[1:n.cores]
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      return(invisible(library(glmmTMB)))
    })
    x.all.res=parLapply(cl=cl, X=all.models, fun=model.fun, idata=ii.data, contr.=contr, xfamily=xfamily, xweights="xweights")
    parLapply(cl=cl, X=1:length(cl), fun=function(x){invisible(rm(list=ls()))})
    stopCluster(cl)
  }else{
		x.all.res=lapply(X=all.models, FUN=model.fun, idata=ii.data, contr.=contr, xfamily=xfamily, xweights="weights")
  }
  xclass=unlist(lapply(x.all.res, function(x){class(x)[[1]]}))
  failed.models=model.terms[xclass=="try-error"]
  x.all.res=x.all.res[xclass!="try-error"]
	if(length(x.all.res)>0){
		all.tests=lapply(1:length(x.all.res), function(x){
			#xx=as.data.frame(anova(x.all.res[[x]]$value, x.all.res[[1]]$value))
			ires=c(Chisq=-2*(as.vector(logLik(x.all.res[[x]]))-as.vector(logLik(model.res))), 
				"Chi Df"=length(unlist(fixef(model.res)))-length(unlist(fixef(x.all.res[[x]]))))
			ires=c(ires, "Pr(>Chisq)"=as.vector(pchisq(q=ires["Chisq"], df=ires["Chi Df"], lower.tail=F)))
			#browser()
			#ires$AIC=xx[1, "AIC"]
			return(list(logLik=as.vector(logLik(x.all.res[[x]])), AIC=extractAIC(x.all.res[[x]])[2], 
				Chisq=as.vector(ires["Chisq"]), "Chi Df"=as.vector(ires["Chi Df"]), 
				"Pr(>Chisq)"=ires["Pr(>Chisq)"], n=x.all.res[[x]]$modelInfo$nobs,
				conv=x.all.res[[x]]$sdr$pdHess))
		})
		aic=extractAIC(model.res)[2]
		####
		#browser()
		all.tests=data.frame(logLik=c(as.vector(logLik(model.res)), unlist(lapply(all.tests, function(x){x["logLik"]}))),
			AIC=c(aic, unlist(lapply(all.tests, function(x){x["AIC"]}))),
			Chisq=c(NA, unlist(lapply(all.tests, function(x){x["Chisq"]}))),
			"Chi Df"=c(NA, unlist(lapply(all.tests, function(x){x["Chi Df"]}))),
			"Pr(>Chisq)"=c(NA, unlist(lapply(all.tests, function(x){x["Pr(>Chisq)"]}))),
			n=c(model.res$modelInfo$nobs, unlist(lapply(all.tests, function(x){x["n"]}))),
			conv=c(model.res$sdr$pdHess, unlist(lapply(all.tests, function(x){x["conv"]}))))
		rownames(all.tests)=paste(c("none", c(model.terms, drop.from.zi, drop.from.disp)[xclass!="try-error"]), c("none", names(all.models)), sep="@")
		rownames(all.tests)=gsub(rownames(all.tests), pattern="none@none", replacement="none", fixed=T)
		if(length(failed.models)>0){
			xx=matrix(NA, ncol=ncol(all.tests), nrow=length(failed.models))
			colnames(xx)=names(all.tests)
			rownames(xx)=failed.models
			all.tests=rbind(all.tests, xx)
		}
		#all.tests[ncol(all.tests)]=as.character(all.tests[, ncol(all.tests)])
		if(return.model.results){
			to.return=list(drop1.res=all.tests, model.results=x.all.res)
		}else{
			to.return=list(drop1.res=all.tests, model.results=NULL)
		}
	}else{
		to.return=list(drop1.res=NULL, model.results=NULL)
	}
  return(to.return)
}

drop1p.glmmtmb.new<-function(m, para=F, data=NULL, n.cores=c("all-1", "all"), to.del=NULL,
  return.mults=F){
  #written Nov. 2018
  ##last changed May 2021 (added dropping from the zi- and dispformula; but not thoroughly tested yet)
  ##determine terms that can be dropped:
  xx=as.character(m$call)
  names(xx)=names(m$call)
  model.fe.re=xx["formula"]
  model.zi=as.formula(xx["ziformula"])
  model.disp=as.formula(xx["dispformula"])

  to.drop.fun<-function(mf){
    #xcall=gsub(x=xcall, pattern="\n", replacement=T)
    model.terms=attr(terms(as.formula(mf)), "term.labels")
    model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
    model.terms=strsplit(model.terms, split=":", fixed=T)
    to.del=unlist(lapply(1:length(model.terms), function(m1){
      xx=unlist(lapply((1:length(model.terms))[-m1], function(m2){
        length(intersect(model.terms[[m1]], model.terms[[m2]]))==length(model.terms[[m1]])
      }))
      sum(xx)>0
    }))
    model.terms=model.terms[!to.del]
    model.terms=unlist(lapply(model.terms, paste, collapse=":"))
    return(model.terms)
  }
  xcall=as.character(model.fe.re)
  if(is.null(to.del)){
    model.terms=to.drop.fun(mf=model.fe.re)
  }else{
    model.terms=to.del
  }
  drop.from.zi=NULL
  drop.from.disp=NULL
  if(!all(as.character(model.zi)%in%c("~", "0", "1"))){
    drop.from.zi=to.drop.fun(mf=as.character(model.zi))
  }
  if(!all(as.character(model.disp)%in%c("~", "0", "1"))){
    drop.from.disp=to.drop.fun(mf=as.character(model.disp))
  }

  model.terms2=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
  ##prepare data:
  all.models=cbind(formula=paste(as.character(model.fe.re), model.terms, sep="-"),
                   ziformula=paste(as.character(model.zi), collapse=""), dispformula=paste(as.character(model.disp), collapse=""))
  if(!is.null(drop.from.zi)){
    all.models=rbind(all.models, cbind(formula=as.character(model.fe.re),
                                       ziformula=paste(paste(as.character(model.zi), collapse=""), drop.from.zi, sep="-"),
                                       dispformula=paste(as.character(model.disp), collapse="")))
  }
  if(!is.null(drop.from.disp)){
    all.models=rbind(all.models, cbind(formula=as.character(model.fe.re),
                                       ziformula=paste(as.character(model.zi), collapse=""),
                                       dispformula=paste(paste(as.character(model.disp), collapse=""), drop.from.disp, sep="-")))
  }
  all.models=lapply(1:nrow(all.models), function(x){
    return(c(formula=as.vector(all.models[x, "formula"]), ziformula=as.vector(all.models[x, "ziformula"]),
             dispformula=as.vector(all.models[x, "dispformula"])))
  })
  names(all.models)=rep(c("cond", "zi", "disp"), times=c(length(model.terms), length(drop.from.zi), length(drop.from.disp)))
  ##define function doing the model comparison:
  model.fun<-function(i.models, m){
    #i.model=as.formula(i.models)
    i.res=try(update(m, formula=as.formula(i.models["formula"]), ziformula=as.formula(i.models["ziformula"]),
                      dispformula=as.formula(i.models["dispformula"])), silent=T)
    return(i.res)
  }
  ##run all models:
  if(para){
    require(parallel)
    n.cores=n.cores[1]
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    if(n.cores=="all"){
      n.cores=length(cl)
    }else if(n.cores=="all-1"){
      if(n.cores=="all-1"){n.cores=length(cl)-1}#else{n.cores=length(cl)}
    }
    n.cores=min(c(n.cores, length(all.models)))
    cl=cl[1:n.cores]
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      return(invisible(library(glmmTMB)))
    })
    x.all.res=parLapply(cl=cl, X=all.models, fun=model.fun, m=m)
    parLapply(cl=cl, X=1:length(cl), fun=function(x){invisible(rm(list=ls()))})
    stopCluster(cl)
  }else{
    x.all.res=lapply(X=all.models, FUN=model.fun, m=m)
  }
  xclass=unlist(lapply(x.all.res, function(x){class(x)[[1]]}))
  failed.models=model.terms[xclass=="try-error"]
  x.all.res=x.all.res[xclass!="try-error"]
  if(length(x.all.res)>0){
    all.tests=lapply(1:length(x.all.res), function(x){
      #xx=as.data.frame(anova(x.all.res[[x]]$value, x.all.res[[1]]$value))
      ires=c(Chisq=-2*(as.vector(logLik(x.all.res[[x]]))-as.vector(logLik(m))),
             "Chi Df"=length(unlist(fixef(m)))-length(unlist(fixef(x.all.res[[x]]))))
      ires=c(ires, "Pr(>Chisq)"=as.vector(pchisq(q=ires["Chisq"], df=ires["Chi Df"], lower.tail=F)))
      #browser()
      #ires$AIC=xx[1, "AIC"]
      return(list(logLik=as.vector(logLik(x.all.res[[x]])), AIC=extractAIC(x.all.res[[x]])[2],
                  Chisq=as.vector(ires["Chisq"]), "Chi Df"=as.vector(ires["Chi Df"]),
                  "Pr(>Chisq)"=ires["Pr(>Chisq)"], n=x.all.res[[x]]$modelInfo$nobs,
                  conv=x.all.res[[x]]$sdr$pdHess))
    })
    aic=extractAIC(m)[2]
    ####
    #browser()
    all.tests=data.frame(logLik=c(as.vector(logLik(m)), unlist(lapply(all.tests, function(x){x["logLik"]}))),
                         AIC=c(aic, unlist(lapply(all.tests, function(x){x["AIC"]}))),
                         Chisq=c(NA, unlist(lapply(all.tests, function(x){x["Chisq"]}))),
                         "Chi Df"=c(NA, unlist(lapply(all.tests, function(x){x["Chi Df"]}))),
                         "Pr(>Chisq)"=c(NA, unlist(lapply(all.tests, function(x){x["Pr(>Chisq)"]}))),
                         n=c(m$modelInfo$nobs, unlist(lapply(all.tests, function(x){x["n"]}))),
                         conv=c(m$sdr$pdHess, unlist(lapply(all.tests, function(x){x["conv"]}))))
    rownames(all.tests)=paste(c("none", c(model.terms, drop.from.zi, drop.from.disp)[xclass!="try-error"]), c("none", names(all.models)), sep="@")
    rownames(all.tests)=gsub(rownames(all.tests), pattern="none@none", replacement="none", fixed=T)
    if(length(failed.models)>0){
      xx=matrix(NA, ncol=ncol(all.tests), nrow=length(failed.models))
      colnames(xx)=names(all.tests)
      rownames(xx)=failed.models
      all.tests=rbind(all.tests, xx)
    }
    #all.tests[ncol(all.tests)]=as.character(all.tests[, ncol(all.tests)])
    if(return.mults){
      to.return=list(drop1.res=all.tests, mults=x.all.res)
    }else{
      to.return=list(drop1.res=all.tests, mults=NULL)
    }
  }else{
    to.return=list(drop1.res=NULL, mults=NULL)
  }
  return(to.return)
}
