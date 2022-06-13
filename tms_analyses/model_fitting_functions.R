

cache.get.fname <- function(varname, base){
  if(is.null(base))
    bname<-tools::file_path_sans_ext(basename(this.file.name()))
  else
    bname<-tools::file_path_sans_ext(basename(base))
  fname<-file.path('cache', 'vars', sprintf("%s_%s.RData",bname,varname))
  fname  
}
if.cached.load <- function(vname, expr,base=bname){
  if(!is.cached.var(vname, base=bname)){
    val=eval.parent(expr)
    assign(vname, val, envir=.GlobalEnv)
    cache.var(vname, bname)
  } else {
    val <- load.cache.var(vname,bname)
  }
  return(val)
}

printf <- function(s, ...){
  cat(sprintf(s, ...))
}
is.cached.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  file.exists(fname)
}
cache.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  print(fname)
  printf("CACHE> saving %s to %s\n", varname, fname)
  dir.create(dirname(fname), showWarnings = FALSE)
  save(list = varname, envir = .GlobalEnv, file = fname)
}

uncache.all <- function(base=NULL){ # remove all cached files. you can also go straight to the directory and delete them
  if(is.null(base))
    bname<-tools::file_path_sans_ext(basename(this.file.name()))
  else
    bname<-tools::file_path_sans_ext(basename(base))
  fnames=list.files(file.path('cache', 'vars'), pattern = sprintf("%s*", bname), full.names=T)
  for(fname in fnames){
    cat(sprintf("Deleting %s\n", fname))
    unlink(fname)
  }
}

load.cache.var <- function(varname, base=NULL){ # load cached objects
  fname<-cache.get.fname(varname, base)
  printf("CACHE> loading %s from %s\n", varname, fname)
  load(fname)
  return(eval(parse(text=varname)))
}

loo_wrapper <- function(...) { # 
  dots <- list(...)
  if (!"x" %in% names(dots)) {
    names(dots)[1] <- "x"
  }
  do.call(brms::loo, dots)
}
model_weights_wrapper <- function(...) {
  dots <- list(...)
  if (!"x" %in% names(dots)) {
    names(dots)[1] <- "x"
  }
  do.call(brms::model_weights, dots)
}


fit_and_plot <- function(mod.name,frm,load.only=T,plot.only.new=F,init="random", dataset, family = cumulative("probit")){
  #mod.name = formula.name.fname(frm)
  is.new=TRUE
  if(!is.cached.var(mod.name, base=bname)){
    mod <- brm(frm, data = dataset, family =family, init=init, iter = 5000, control = list(adapt_delta = 0.99, max_treedepth = 12)) %>% #control = list(adapt_delta = 0.99, max_treedepth = 12)
      add_criterion(c("loo"))
    assign(mod.name, mod, envir=.GlobalEnv)
    cache.var(mod.name, bname)
  } else {
    mod <- load.cache.var(mod.name,bname)
    is.new=FALSE
  }
  if(!load.only & ((is.new & plot.only.new) | (!plot.only.new))  ){
    pdf(plot.filename(sprintf("diag_%s.pdf", mod.name),bname), width=5, height=5)
    mcmc_rhat(brms::rhat(mod)) %>% print
    mcmc_neff(brms::neff_ratio(mod)) %>% print
    #dev.off()
    
    
    mcmc_intervals_data(as.matrix(mod), prob_outer = 0.95) %>%
      filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
      ggplot(aes(y=m, ymin=ll,ymax=hh,x=parameter))+
      geom_pointrange(position=position_dodge(width=0.2))+
      coord_flip()
    ggsave(plot.filename(sprintf("coef_%s.pdf",mod.name),bname), width=9,height=6)
    
    fit=mod
    nrep=100
    pred=predict(fit)
    
    dataset %>%
      cbind(
        replicate(n=nrep, apply(pred, 1, function(x){sample(1:4,1, prob=x)})) 
      )   %>%
      gather(sim.n,sim.response, 13:(13+nrep-1)) %>%
      group_by(randomization, condition, sim.n) %>%
      do({
        tibble(response=1:4,n=tabulate(.$sim.response, nbins=4))
      }) -> dataset.pred
    
    d.tab=dataset %>% group_by(randomization, condition) %>%
      do({
        v=as.numeric(data.frame(.)[,"probe.response"])
        tibble(response=1:4,n=tabulate(v, nbins=4))
      })
    
    dataset.pred %>% ungroup %>% 
      ggplot(aes(x=factor(response),y=n,color=condition))+
      geom_bar(data=d.tab, mapping=aes(fill=condition), stat="identity",position = position_dodge(width=1), alpha=0.2)+
      #geom_violin(aes(group=interaction(stim_setting,response),color=NULL),fill="grey",color=0, alpha=1, position=position_dodge(width=1))+
      stat_summary(fun.data = mean_qi,  position=position_dodge(width=1), geom="pointrange") +
      #facet_wrap(~question,ncol=1) +
      facet_grid(condition~.) +
      labs(x="Response",y="Number of subjects",
           title=sprintf("%s: Posterior predictive", fit$formula$resp), 
           subtitle=toString(capture.output(fit$formula)))
    
    #ggsave(plot.filename(sprintf("ppred_%s.pdf",mod.name),bname), width=9,height=6)
  }
  return(mod)
}

