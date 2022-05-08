library(brms)
library(bayesplot)
library(tidybayes)
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(tidyverse)
theme_set(theme_classic())
setwd("/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/tms_analyses")

bname="tms_analyses"
#options(mc.cores=parallel::detectCores())
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options) {
  uncache.all(base = bname)
}


############################################### FITTING LIST OF MODELS ON TASK PROBE

models_task <- list(
  formula(probe.response ~ zbv * zlog.apen + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + randomization + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + randomization + (1|subj/condition)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + randomization + (1|subj/condition/visit)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + randomization + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition +  condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num * condition + (1|subj/condition))
)


descriptions_task=c(
  "BV x AE", 
  "BV x AE + trial", 
  "BV x AE + trial + block", 
  "BV x AE + trial + block + condition",
  "BV x AE + trial + block + condition + randomization",
  "BV x AE + trial + block + condition + visit",
  "BV x AE + trial + block + condition + visit + randomization",
  #"BV x AE + trial + block + condition + visit + randomization + visit(subj.nested)",
  #"BV x AE + trial + block + condition + visit + visit(subj.nested)",
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  "BV x AE + trial + block + condition + visit + condition : AE + condition : BV"#,
  #"BV x AE + trial + block + condition + randomization + condition : AE + condition : BV"#,
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  #"BV x AE + trial + block x condition"
  )

conditional_effects(models_task.fitted[7]$mod_task06)


pp_check(model_task_8_fit, "ecdf_overlay") #posterior predictive check
posterior_model_task_5 <- as.array(models_task.fitted[6]$mod_task05)

mcmc_pairs(posterior_model_task_5, pars = c("b_zbv", "b_zlog.apen", "b_conditionactive_rhTMS") , np = np_task_model_5,
           off_diag_args = list(size = 0.75))

mcmc_scatter(
  posterior_model_task_5, 
  pars = c("b_conditionactive_rhTMS", "b_zbv"), 
  np = np_task_model_5, 
  size = 1
)


posterior_model_task_5


log_posterior_model_task_6 <- log_posterior(models_task.fitted[6]$mod_task05)
np_task_model_5 <- nuts_params(models_task.fitted[6]$mod_task05)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_model_task_5, np = np_task_model_5)


pairs(models_task.fitted[6], variable = c("condition"))

names(models_task) <- sprintf("mod_task%02i", 0:(length(models_task)-1))

models_task.wrap <- map2(names(models_task), models_task, ~ list(mod.name=.x, mod=.y))
models_task.fitted=lapply(models_task.wrap, function(lmod){ fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = tms_data.nback)})
names(models_task.fitted) <- names(models_task)

loos_task=if.cached.load("loos",
                         invoke(loo_wrapper, .x = models_task.fitted, model_names = names(models_task.fitted)),
                         base=bname)

r2s=lapply(models_task.fitted, bayes_R2, cl=22)
mod_intention.weights = if.cached.load("mod_task.weights",
                                       map_df(c("loo", "waic", "stacking"), function(strat) {
                                         r = invoke(
                                           model_weights_wrapper,
                                           .x = models_task.fitted,
                                           weights = strat,
                                           model_names = names(models_task.fitted)
                                         )
                                         bind_cols(strategy = strat, data.frame(t(r)))
                                       }), bname)

print(loo_compare(x=loos_task$loos))
as.data.frame(loos_task$ic_diffs__) %>% rownames_to_column() %>% 
  mutate(z=LOOIC/SE) %>% print

print(mod_task.weights)

mod.desc=data.frame(mod=names(models_task.fitted), descriptions_task)
map_df(c("loo","stacking","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions_task) %>%
  #mutate(loo2="",waic="") %>%
  #pivot_longer(names_to = strategy, values_to = descriptions_task)
  gather(strategy,descriptions_task,loo,stacking,waic) -> mod.desc

mod_task.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","stacking")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="stacking" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy!="WAIC" & strategy!="pseudo-BMA") %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions_task), y=0, hjust="left")+
  labs(x="",y="Posterior Probability")+
  facet_wrap(~strategy)+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "inside") 

#ggsave("model_task_weights.png", plot = last_plot(), width=10,height=5)

############################################### MODELLING TASK PERFORMANCE

models_bv <- list(
  formula(zbv ~ probeix + (1|subj/condition)),
  formula(zbv ~ probeix + block_num + (1|subj/condition)),
  formula(zbv ~ probeix + block_num + condition + (1|subj/condition)),
  formula(zbv ~ probeix + block_num + condition + visit + (1|subj/condition))
)


descriptions_bv=c(
  "trial", 
  "trial + block", 
  "trial + block + condition",
  "trial + block + condition + visit"
)


names(models_bv) <- sprintf("mod_bv%02i", 0:(length(models_bv)-1))

models_bv.wrap <- map2(names(models_bv), models_bv, ~ list(mod.name=.x, mod=.y))
models_bv.fitted=lapply(models_bv.wrap, function(lmod){ fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = tms_data.nback, family = student(link = "identity", link_sigma = "log", link_nu = "logm1"))})
names(models_bv.fitted) <- names(models_bv)

loos_bv=if.cached.load("loos_bv",
                       invoke(loo_wrapper, .x = models_bv.fitted, model_names = names(models_bv.fitted)),
                       base=bname)

lapply(models_bv.fitted, bayes_R2, cl=22)
mod_bv.weights = if.cached.load("mod_bv.weights",
                                       map_df(c("loo", "waic", "stacking"), function(strat) {
                                         r = invoke(
                                           model_weights_wrapper,
                                           .x = models_bv.fitted,
                                           weights = strat,
                                           model_names = names(models_bv.fitted)
                                         )
                                         bind_cols(strategy = strat, data.frame(t(r)))
                                       }), bname)

print(loo_compare(x=loos_bv$loos))
as.data.frame(loos_bv$ic_diffs__) %>% rownames_to_column() %>% 
  mutate(z=LOOIC/SE) %>% print

print(mod_bv.weights)

mod.desc=data.frame(mod=names(models_bv.fitted), descriptions_bv)
map_df(c("loo","stacking","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions_bv) %>%
  #mutate(loo2="",waic="") %>%
  #pivot_longer(names_to = strategy, values_to = descriptions_bv)
  gather(strategy,descriptions_bv,loo,stacking,waic) -> mod.desc

mod_bv.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","stacking")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="stacking" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy!="WAIC" & strategy!="pseudo-BMA") %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions_bv), y=0, hjust="left")+
  labs(x="",y="Posterior Probability")+
  facet_wrap(~strategy)+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "inside") 

## BV: MCMC Intervals

color_scheme_set("viridisE")

mcmc_intervals(models_bv.fitted$mod_bv03, pars = c("b_block_num",
                                                       "b_probeix",
                                                       "b_visit2",
                                                       "b_conditionsham_rhTMS",
                                                       "b_conditionsham_arrhTMS",
                                                       "b_conditionactive_rhTMS",
                                                       "b_conditionactive_arrhTMS"
)) +
  ggplot2::scale_y_discrete(labels = c("b_block_num" = "Block",
                                       "b_probeix" = "Probe number",
                                       "b_visit2" = "Visit",
                                       "b_conditionsham_rhTMS" = "Sham rhTMS",
                                       "b_conditionsham_arrhTMS" = "Sham arrhTMS",
                                       "b_conditionactive_rhTMS" = "Active rhTMS",
                                       "b_conditionactive_arrhTMS" = "Active arrhTMS"))


################################### FINTTING MODELS ON AE

models_apen <- list(
  formula(zlog.apen ~ probeix + (1|subj/condition)),
  formula(zlog.apen ~ probeix + block_num + (1|subj/condition)),
  formula(zlog.apen ~ probeix + block_num + condition + (1|subj/condition)),
  formula(zlog.apen ~ probeix + block_num + condition + visit + (1|subj/condition))
)


descriptions_apen=c(
  "trial", 
  "trial + block", 
  "trial + block + condition",
  "trial + block + condition + visit"
)


names(models_apen) <- sprintf("mod_apen%02i", 0:(length(models_apen)-1))

models_apen.wrap <- map2(names(models_apen), models_apen, ~ list(mod.name=.x, mod=.y))
models_apen.fitted=lapply(models_apen.wrap, function(lmod){ fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = tms_data.nback, family = student(link = "identity", link_sigma = "log", link_nu = "logm1"))})
names(models_apen.fitted) <- names(models_apen)

loos_apen=if.cached.load("loos_apen",
                         invoke(loo_wrapper, .x = models_apen.fitted, model_names = names(models_apen.fitted)),
                         base=bname)

lapply(models_apen.fitted, bayes_R2, cl=22)
mod_apen.weights = if.cached.load("mod_apen.weights",
                                  map_df(c("loo", "waic", "stacking"), function(strat) {
                                    r = invoke(
                                      model_weights_wrapper,
                                      .x = models_apen.fitted,
                                      weights = strat,
                                      model_names = names(models_apen.fitted)
                                    )
                                    bind_cols(strategy = strat, data.frame(t(r)))
                                  }), bname)

print(loo_compare(x=loos_apen$loos))
as.data.frame(loos_apen$ic_diffs__) %>% rownames_to_column() %>% 
  mutate(z=LOOIC/SE) %>% print

print(mod_apen.weights)

mod.desc=data.frame(mod=names(models_apen.fitted), descriptions_apen)
map_df(c("loo","stacking","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions_apen) %>%
  #mutate(loo2="",waic="") %>%
  #pivot_longer(names_to = strategy, values_to = descriptions_apen)
  gather(strategy,descriptions_apen,loo,stacking,waic) -> mod.desc

mod_apen.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","stacking")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="stacking" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy!="WAIC" & strategy!="pseudo-BMA") %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions_apen), y=0, hjust="left")+
  labs(x="",y="Posterior Probability")+
  facet_wrap(~strategy)+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "inside") 

## AE: MCMC Intervals

color_scheme_set("viridisE")

mcmc_intervals(models_apen.fitted$mod_apen03, pars = c("b_block_num",
                                                       "b_probeix",
                                                       "b_visit2",
                                                       "b_conditionsham_rhTMS",
                                                       "b_conditionsham_arrhTMS",
                                                       "b_conditionactive_rhTMS",
                                                       "b_conditionactive_arrhTMS"
)) +
  ggplot2::scale_y_discrete(labels = c("b_block_num" = "Block",
                                       "b_probeix" = "Probe number",
                                       "b_visit2" = "Visit",
                                       "b_conditionsham_rhTMS" = "Sham rhTMS",
                                       "b_conditionsham_arrhTMS" = "Sham arrhTMS",
                                       "b_conditionactive_rhTMS" = "Active rhTMS",
                                       "b_conditionactive_arrhTMS" = "Active arrhTMS"))




## LOOIC
loo_compare(x = loos$loos) %>% data.frame %>% rownames_to_column(var = "mod") %>%
  arrange(mod) %>% 
  mutate(rel.looic=looic-first(looic)) %>%
  ggplot(aes(mod,rel.looic))+
  geom_bar(stat="identity", fill="grey")+
  geom_text(
    aes(fill = NULL, label = frm),
    y = 0,
    data = data.frame(mod = names(models_task), frm =
                        map_chr(
                          models_task, ~ paste(format(.x), sep = "", collapse = "")
                        )),
    hjust = "right"
  ) +
  #geom_pointrange(aes(ymin=rel.looic-se_looic,ymax=rel.looic+se_looic))+
  coord_flip()

ggsave("task_rel_looic.png", last_plot(), width=12,height=9)



#========================
## inspection of the model
#========================
stop()
test.models = c("mod_task05", "mod_task07")
map2_df(test.models,
        models_task.fitted[test.models],
        ~ bind_cols(
          mod = .x,
          hypothesis(.y, "conditionactive_rhTMS > 0", alpha = .05)$hypothesis
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
  geom_point(aes(y = ppred, size f = probpred, color=NULL),
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

####
mod=models_task.fitted$mod_task05
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  ggplot(aes(y=m, ymin=ll,ymax=hh,x=parameter))+
  geom_pointrange(position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red")+
  labs(y="Coefficient")

## MCMC Intervals
color_scheme_set("viridisE")
mcmc_intervals(models_task.fitted$mod_task05, pars = c("b_zbv",
                         "b_zlog.apen",
                         "b_block_num",
                         "b_probeix",
                         "b_zbv:zlog.apen",
                         "b_visit2",
                         "b_conditionsham_rhTMS",
                         "b_conditionsham_arrhTMS",
                         "b_conditionactive_rhTMS",
                         "b_conditionactive_arrhTMS"
                         )) +
  ggplot2::scale_y_discrete(labels = c("b_zbv" = "BV",
                                       "b_zlog.apen" = "AE",
                                       "b_block_num" = "Block",
                                       "b_probeix" = "Probe number",
                                       "b_zbv:zlog.apen" = "BV x AE",
                                       "b_visit2" = "Visit",
                                       "b_conditionsham_rhTMS" = "Sham rhTMS",
                                       "b_conditionsham_arrhTMS" = "Sham arrhTMS",
                                       "b_conditionactive_rhTMS" = "Active rhTMS",
                                       "b_conditionactive_arrhTMS" = "Active arrhTMS"))
 #labs(y="Coefficient",x="Predictor") +
#geom_hline(yintercept = 0.0, color="red", size=2, alpha=0.2)


######Model with interaction measures x stim

mcmc_intervals(models_task.fitted$mod_task07, pars = c("b_zbv",
                                                       "b_zlog.apen",
                                                       "b_block_num",
                                                       "b_probeix",
                                                       "b_zbv:zlog.apen",
                                                       #"b_visit2",
                                                       "b_conditionsham_rhTMS",
                                                       "b_conditionsham_arrhTMS",
                                                       "b_conditionactive_rhTMS",
                                                       "b_conditionactive_arrhTMS",
                                                       "b_zlog.apen:conditionactive_rhTMS",
                                                       "b_zlog.apen:conditionsham_rhTMS",
                                                       "b_zlog.apen:conditionactive_arrhTMS",
                                                       "b_zlog.apen:conditionsham_arrhTMS",
                                                       "b_zbv:conditionactive_rhTMS",
                                                       "b_zbv:conditionsham_rhTMS",
                                                       "b_zbv:conditionactive_arrhTMS",
                                                       "b_zbv:conditionsham_arrhTMS"
                                                       )) +
  ggplot2::scale_y_discrete(labels = c("b_zbv" = "BV",
                                       "b_zlog.apen" = "AE",
                                       "b_block_num" = "Block",
                                       "b_probeix" = "Probe number",
                                       "b_zbv:zlog.apen" = "BV x AE",
                                       #"b_visit2" = "Visit",
                                       "b_conditionsham_rhTMS" = "Sham rhTMS",
                                       "b_conditionsham_arrhTMS" = "Sham arrhTMS",
                                       "b_conditionactive_rhTMS" = "Active rhTMS",
                                       "b_conditionactive_arrhTMS" = "Active arrhTMS",
                                       "b_zlog.apen:conditionactive_rhTMS" = "AE x Active rhTMS",
                                       "b_zlog.apen:conditionsham_rhTMS" = "AE x Sham rhTMS",
                                       "b_zlog.apen:conditionactive_arrhTMS" = "AE x Active arrhTMS",
                                       "b_zlog.apen:conditionsham_arrhTMS" = "AE x Sham arrhTMS",
                                       "b_zbv:conditionactive_rhTMS" = "BV x Active rhTMS",
                                       "b_zbv:conditionsham_rhTMS"= "BV x Sham rhTMS",
                                       "b_zbv:conditionactive_arrhTMS" = "BV x Active arrhTMS",
                                       "b_zbv:conditionsham_arrhTMS" = "BV x Sham arrhTMS"))
               
## nya SfN
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9 ) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj"), parameter != "disc") %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  mutate(parameter=fct_relevel(parameter, "b_zbv:zlog.apen", after=5)) %>%
  mutate(parameter=fct_recode(parameter, 
                              `BV`="b_zbv",
                              `AE`="b_zlog.apen",
                              `Block`="b_block_num",
                              `Probe number`="b_probeix",
                              `BV x AE`="b_zbv:zlog.apen",
                              `Visit` = "b_visit2",
                              `Sham rhTMS` = "b_conditionsham_rhTMS",
                              `Sham arrhTMS` = "b_conditionsham_arrhTMS",
                              `Active rhTMS` = "b_conditionactive_rhTMS",
                              `Baseline` = "b_conditionbaseline")) %>%
  ggplot(aes(y=m, x=parameter))+
  geom_pointrange(aes(ymin=l,ymax=h), position=position_dodge(width=0.2), color="black", size=2,fatten=0.8)+
  geom_pointrange(aes(ymin=ll,ymax=hh), position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red", size=2, alpha=0.2)+
  labs(y="Coefficient",x="Predictor")+
  #annotate("text", y=-0.5, x=7.25, label="P(b<0)=0.97, Evidence Ratio=30.7",hjust = 0)+
  theme(axis.text.y = element_text(angle=0, hjust=1, size=14),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16)) -> p1
p1
ggsave("graphs/blockxstim.pdf", width=7, height=4)


####################################################### MODELLING INTENTION

models_intention <- list(
  formula(intention ~ zbv * zlog.apen + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + randomization + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + visit + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + visit + randomization + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + visit + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
  #formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + randomization + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
  #formula(intention ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(intention ~ zbv * zlog.apen + probeix + block_num + condition +  condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(intention ~ zbv * zlog.apen + probeix + block_num * condition + (1|subj/condition))
)


descriptions_intention=c(
  "BV x AE", 
  "BV x AE + trial", 
  "BV x AE + trial + block", 
  "BV x AE + trial + block + condition",
  "BV x AE + trial + block + condition + randomization",
  "BV x AE + trial + block + condition + visit",
  "BV x AE + trial + block + condition + visit + randomization",
  #"BV x AE + trial + block + condition + visit + randomization + visit(subj.nested)",
  #"BV x AE + trial + block + condition + visit + visit(subj.nested)",
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  "BV x AE + trial + block + condition + condition : AE + condition : BV",
  "BV x AE + trial + block + condition + visit + condition : AE + condition : BV"
  #,
  #"BV x AE + trial + block + condition + randomization + condition : AE + condition : BV"#,
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  #"BV x AE + trial + block x condition"
)

############################################### FITTING LIST OF MODELS ON TASK PROBE

names(models_intention) <- sprintf("mod_intention%02i", 0:(length(models_intention)-1))

models_intention.wrap <- map2(names(models_intention), models_intention, ~ list(mod.name=.x, mod=.y))
models_intention.fitted=lapply(models_intention.wrap, function(lmod){fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = tms_data.nback)})
names(models_intention.fitted) <- names(models_intention)

loos_intention=if.cached.load("loos_intention",
                              invoke(loo_wrapper, .x = models_intention.fitted, model_names = names(models_intention.fitted)),
                              base=bname)

r2s=lapply(models_intention.fitted, bayes_R2, cl=22)
mod_intention.weights = if.cached.load("mod.weights",
                                       map_df(c("loo", "waic", "stacking"), function(strat) {
                                         r = invoke(
                                           model_weights_wrapper,
                                           .x = models_intention.fitted,
                                           weights = strat,
                                           model_names = names(models_intention.fitted)
                                         )
                                         bind_cols(strategy = strat, data.frame(t(r)))
                                       }), bname)

print(loos_intention)
print(loo_compare(x=loos_intention$loos))
as.data.frame(loos_intention$ic_diffs__) %>% rownames_to_column() %>% 
  mutate(z=LOOIC/SE) %>% print

print(mod_intention.weights)

mod.desc=data.frame(mod=names(models_intention.fitted), descriptions_intention)
map_df(c("loo","stacking","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions_intention) %>%
  #mutate(loo2="",waic="") %>%
  gather(strategy,descriptions_intention,loo,stacking,waic) -> mod.desc

mod.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","stacking")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="stacking" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy==c("LOO")) %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions_intention), y=0, hjust="left")+
  labs(x="",y="Posterior Probability")+
  facet_wrap(~strategy)+
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "inside") 

ggsave("model_weights.png", plot = last_plot(), width=10,height=5)

mod.weights[,-1] %>% as.matrix %>% t %>% data.frame %>% 
  setNames(mod.weights$strategy) %>%
  rownames_to_column() -> modw.df
modw.df %>% arrange(desc(loo)) %>% head(2)
modw.df %>% arrange(desc(loo2)) %>% head(2)



####################################################### MODELLING ALERTNESS


models_somnolence <- list(
  formula(somnolence ~ zbv * zlog.apen + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + randomization + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + visit + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + visit + randomization + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + visit + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
  #formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + randomization + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
  #formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(somnolence ~ zbv * zlog.apen + probeix + block_num + condition +  condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(somnolence ~ zbv * zlog.apen + probeix + block_num * condition + (1|subj/condition))
)


descriptions_somnolence=c(
  "BV x AE", 
  "BV x AE + trial", 
  "BV x AE + trial + block", 
  "BV x AE + trial + block + condition",
  "BV x AE + trial + block + condition + randomization",
  "BV x AE + trial + block + condition + visit",
  "BV x AE + trial + block + condition + visit + randomization",
  #"BV x AE + trial + block + condition + visit + randomization + visit(subj.nested)",
  #"BV x AE + trial + block + condition + visit + visit(subj.nested)",
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  "BV x AE + trial + block + condition + condition : AE + condition : BV",
  "BV x AE + trial + block + condition + visit + condition : AE + condition : BV"
  #,
  #"BV x AE + trial + block + condition + randomization + condition : AE + condition : BV"#,
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  #"BV x AE + trial + block x condition"
)


names(models_somnolence) <- sprintf("mod_somnolence%02i", 0:(length(models_somnolence)-1))

models_somnolence.wrap <- map2(names(models_somnolence), models_somnolence, ~ list(mod.name=.x, mod=.y))
models_somnolence.fitted=lapply(models_somnolence.wrap, function(lmod){fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = tms_data.nback)})
names(models_somnolence.fitted) <- names(models_somnolence)

loos_somnolence=if.cached.load("loos_somnolence",
                               invoke(loo_wrapper, .x = models_somnolence.fitted, model_names = names(models_somnolence.fitted)),
                               base=bname)

r2s_somnolence=lapply(models_somnolence.fitted, bayes_R2, cl=22)
mod_somnolence.weights = if.cached.load("mod.weights",
                                        map_df(c("loo", "waic", "stacking"), function(strat) {
                                          r = invoke(
                                            model_weights_wrapper,
                                            .x = models_somnolence.fitted,
                                            weights = strat,
                                            model_names = names(models_somnolence.fitted)
                                          )
                                          bind_cols(strategy = strat, data.frame(t(r)))
                                        }), bname)

mod.desc=data.frame(mod=names(models_somnolence.fitted), descriptions_somnolence)
map_df(c("loo","stacking","waic"), ~ cbind(strategy=.x,mod.desc)) %>%
  spread(strategy,descriptions_somnolence) %>%
  #mutate(loo2="",waic="") %>%
  gather(strategy,descriptions_somnolence,loo,stacking,waic) -> mod.desc

mod.weights %>%
  gather(mod, prob, starts_with("mod")) %>% 
  full_join(mod.desc) %>%
  mutate(
    strategy=ordered(strategy, c("loo", "waic","stacking")),
    strategy=ordered(case_when(strategy=="loo" ~ "LOO",
                               strategy=="waic" ~ "WAIC",
                               strategy=="stacking" ~ "pseudo-BMA"),
                     c("LOO","WAIC","pseudo-BMA"))) %>%
  filter(strategy==c("LOO")) %>% droplevels %>%
  group_by(strategy) %>%
  mutate(win=if_else(prob==max(prob), T,F)) %>%
  ungroup %>%
  ggplot(aes(mod, prob, fill = win)) +
  geom_bar(stat = "identity", position =
             position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("lightblue", "orange"))+
  geom_text(mapping=aes(fill=NULL, label=descriptions_somnolence), y=0, hjust="left")+
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

ggsave("rel_looic.png", last_plot(), width=12,height=9)



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
  ggplot(aes(probeix, intention, color = part)) +
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

####
mod=models.fitted[["mod22"]]
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9 ) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  ggplot(aes(y=m, ymin=ll,ymax=hh,x=parameter))+
  geom_pointrange(position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red")+
  labs(y="Coefficient")

## nya SfN
mcmc_intervals_data(as.matrix(mod), prob_outer = 0.9 ) %>%
  filter(parameter!="lp__", !str_detect(parameter, "subj")) %>%
  filter(!str_detect(parameter, "Intercept")) %>%
  mutate(parameter=fct_relevel(parameter, "b_zbv:zlog.apen", after=5)) %>%
  mutate(parameter=fct_recode(parameter, 
                              `Variability`="b_zbv",
                              `Entropy`="b_zlog.apen",
                              `Block`="b_partstimulation",
                              `Trial`="b_probeix",
                              `Stimulation`="b_conditionreal",
                              `Variability x Entropy`="b_zbv:zlog.apen",
                              `Block x Stimulation`="b_partstimulation:conditionreal")) %>%
  ggplot(aes(y=m, x=parameter))+
  geom_pointrange(aes(ymin=l,ymax=h), position=position_dodge(width=0.2), color="black", size=2,fatten=0.8)+
  geom_pointrange(aes(ymin=ll,ymax=hh), position=position_dodge(width=0.2))+
  coord_flip()+geom_hline(yintercept = 0, color="red", size=2, alpha=0.2)+
  labs(y="Coefficient",x="Predictor")+
  #annotate("text", y=-0.5, x=7.25, label="P(b<0)=0.97, Evidence Ratio=30.7",hjust = 0)+
  theme(axis.text.y = element_text(angle=0, hjust=1, size=14),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16)) -> p1
p1
ggsave("graphs/blockxstim.pdf", width=7, height=4)
