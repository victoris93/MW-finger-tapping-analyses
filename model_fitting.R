library(brms)
library(bayesplot)
library(tidybayes)
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())

all_data_toni_nback <- all_data_toni_nback %>%
  mutate(condition = case_when(str_detect(condition, "baseline") ~ "baseline",
                               TRUE ~ as.character(condition)))

bname="model_fitting"
dir.create(bname)
#bname="model_fitting"
#options(mc.cores=parallel::detectCores())
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options) {
  uncache.all(base = bname)
}

this.file.name <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

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
uncache.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  cat(sprintf("Deleting %s\n", fname))
  unlink(fname)
}

uncache.all <- function(base=NULL){
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

load.cache.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  printf("CACHE> loading %s from %s\n", varname, fname)
  load(fname)
  return(eval(parse(text=varname)))
}

printf <- function(s, ...){
  cat(sprintf(s, ...))
}

loo_wrapper <- function(...) {
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

fit_and_plot <- function(mod.name,frm,load.only=T,plot.only.new=F,init="random", dataset){
  #mod.name = formula.name.fname(frm)
  is.new=TRUE
  if(!is.cached.var(mod.name, base=bname)){
    mod <- brm(frm, data = dataset, family =cumulative("probit"), init=init) %>%
      add_loo() %>% add_waic()
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
    dev.off()
    
    
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
      group_by(condition, sim.n) %>%
      do({
        tibble(response=1:4,n=tabulate(.$sim.response, nbins=4))
      }) -> dataset.pred
    
    d.tab=dataset %>% group_by(condition) %>%
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
    
    ggsave(plot.filename(sprintf("ppred_%s.pdf",mod.name),bname), width=9,height=6)
  }
  return(mod)
}


models <- list(
  formula(probe.response ~ zbv:zlog.apen + zbv + zlog.apen + condition + probeix + (1|subj/condition)),
  formula(probe.response ~ zbv:zlog.apen + zbv + zlog.apen + probeix:condition+ condition + probeix  + (1|subj/condition)),
  formula(probe.response ~ zbv:zlog.apen + zbv + zlog.apen + probeix:condition+ condition + probeix  + block_num + (1|subj/condition)),
  formula(probe.response ~ zbv:zlog.apen + zbv + zlog.apen + probeix:condition+ condition + probeix + block_num + block_num:condition + block_num:probeix + (1|subj/condition))
)

descriptions=c( "BV x AE + BV + AE + condition + trial",
                "BV x AE + BV + AE + condition x trial + condition x trial",
                "BV x AE + BV + AE + condition x trial + condition + trial + block",
                "BV x AE + BV + AE + condition x trial + condition + trial + block + block x condition + block x trial"
                )

names(models) <- sprintf("model_%i", 0:(length(models)-1))
models.wrap <- map2(names(models), models, ~ list(mod.name=.x, mod=.y))
models.fitted=lapply(models.wrap, function(lmod){ fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = all_data_toni_nback)})
names(models.fitted) <- names(models)
loos=if.cached.load("loos",
                    invoke(loo_wrapper, .x = models.fitted, model_names = names(models.fitted)),
                    base=bname)
mod.weights = if.cached.load("mod.weights",
                             map_df(c("loo", "waic", "loo2"), function(strat) {
                               r = invoke(
                                 model_weights_wrapper,
                                 .x = models.fitted,
                                 weights = strat,
                                 model_names = names(models.fitted)
                               )
                               bind_cols(strategy = strat, data.frame(t(r)))
                             }), bname)
r2s=lapply(models.fitted, bayes_R2, cl=22)
#print(loo::compare(x=loos$loos))
#as.data.frame(loos$ic_diffs__) %>% rownames_to_column() %>% 
#  mutate(z=LOOIC/SE) %>% print


print(mod.weights)

mod.desc=data.frame(mod=names(models.fitted), descriptions)
  map_df(c("loo","loo2","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions) %>%
  #mutate(loo2="",waic="") %>%
  gather(strategy,descriptions,loo,loo2,waic) -> mod.desc

### MODEL SELECTION  
mod.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","loo2")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="loo2" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy!="WAIC") %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions), y=0, hjust="left")+
  labs(x="",y="Posterior Probability")+
  facet_wrap(~strategy)+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "inside") 

## LOOIC
loo_compare(x = loos$loos) %>% data.frame %>% rownames_to_column(var = "mod") %>%
  arrange(mod) %>% 
  mutate(rel.looic=looic-first(looic)) %>%
  ggplot(aes(mod,rel.looic))+
  geom_bar(stat="identity", fill="grey")+
  geom_text(
    aes(fill = NULL, label = frm),
    y = 0,
    data = data.frame(mod = names(models), frm =
                        map_chr(
                          models, ~ paste(format(.x), sep = "", collapse = "")
                        )),
    hjust = "right"
  ) +
  #geom_pointrange(aes(ymin=rel.looic-se_looic,ymax=rel.looic+se_looic))+
  coord_flip()

mod=models.fitted[["model_1"]]
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  ggplot(aes(y=m, ymin=ll,ymax=hh,x=parameter))+
  geom_pointrange(position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red")+
  labs(y="Coefficient") +
  theme()

#========================
## inspection of the model
#========================
stop()
test.models = c("mod10", "mod09", "mod13", "mod11", "mod12", "mod14", "mod15", "mod16")
map2_df(test.models,
        models.fitted[test.models],
        ~ bind_cols(
          mod = .x,
          hypothesis(.y, "partstimulation:conditionreal<0", alpha = .05)$hypothesis
        ))


## all models that have the interaction
map2_df(names(models.fitted),
        models.fitted,
        function(modname,mod){
          if("partstimulation:conditionreal" %in% row.names(fixef(mod))){
            return(bind_cols(mod=modname,hypothesis(mod, "partstimulation:conditionreal<0", alpha = .05)$hypothesis))
          } else{
            return(NULL);
          }
        }) -> d.ia
d.ia

d.ia$Evid.Ratio %>% summary

summary(mod22)$fixed %>% data.frame %>% rownames_to_column() %>%
  mutate(summary=sprintf("$b=%.2f\\ [%.2f, %.2f]$", Estimate, l.95..CI, u.95..CI)) %>%
  select(rowname,summary)

hypothesis(mod22e, c("zbv>0", "zlog.apen<0", "partstimulation>0", "probeix>0", 
                     "conditionreal>0", "zbv:zlog.apen>0", "partstimulation:conditionreal<0"))
pred=with(models.fitted, predict(mod22))
library(ggforce)
d.nback %>% bind_cols(data.frame(pred) %>% setNames(c("pprob1","pprob2","pprob3","pprob4"))) %>% 
  #filter(subj==4) %>%
  arrange(part,probeix) %>%
  mutate(probeix=if_else(part=="stimulation", probeix+9, probeix)) %>%
  group_by(condition,part,probeix) %>% 
  mutate(cond.subj=1:n()) %>% ungroup %>%
  arrange(subj) -> d.tmp


d.tmp %>%
  ggplot(aes(probeix, probe.response, color = part)) +
  geom_point(aes(y = ppred, size f= probpred, color=NULL),
             color="grey", alpha=0.2,
             data = d.tmp %>%
               gather(ppred, probpred, starts_with("pprob")) %>%
               separate(ppred,into = c("blub","ppred"), sep = 5) %>%
               mutate(ppred=as.integer(ppred))
  )+
  geom_point(aes(y = bestpred),
             color="orange", alpha=0.2, size=10,
             data = bind_cols(d.tmp, bestpred=apply(d.tmp[13:16], 1, which.max))
  )+
  geom_point()+theme_bw() -> p

npages=n_pages(p+facet_grid_paginate(cond.subj~condition, nrow=5, ncol=2,page=1))


pdf(plot.filename("ppred_subj_mod22.pdf", bname))
map(1:npages,
    ~ print(
      p + facet_grid_paginate(
        cond.subj ~ condition,
        nrow = 5,
        ncol = 2,
        page = .x
      )
    ))
dev.off()