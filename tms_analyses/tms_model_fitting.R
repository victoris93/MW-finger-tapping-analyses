library(brms)
library(bayesplot)
library(tidybayes)
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
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
  formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + randomization + (1|subj/condition))#,
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + (1|subj/condition)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + randomization + (1|subj/condition))#,
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + randomization + (1|subj/condition/visit)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + condition:zlog.apen + condition:zbv + (1|subj/condition)),
  #formula(probe.response ~ zbv * zlog.apen + probeix + block_num + condition + visit + condition:zlog.apen + condition:zbv + (1|subj/condition))#,
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
  "BV x AE + trial + block + condition + randomization"#,
  #"BV x AE + trial + block + condition + visit",
  #"BV x AE + trial + block + condition + visit + randomization"#,
  #"BV x AE + trial + block + condition + visit + randomization + visit(subj.nested)",
  #"BV x AE + trial + block + condition + visit + visit(subj.nested)",
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  #"BV x AE + trial + block + condition + visit + condition : AE + condition : BV"#,
  #"BV x AE + trial + block + condition + randomization + condition : AE + condition : BV"#,
  #"BV x AE + trial + block + condition + condition : AE + condition : BV",
  #"BV x AE + trial + block x condition"
  )

conditional_effects(models_task.fitted[7]$mod_task06)

names(models_task) <- sprintf("mod_task%02i", 0:(length(models_task)-1))

models_task.wrap <- map2(names(models_task), models_task, ~ list(mod.name=.x, mod=.y))
models_task.fitted=lapply(models_task.wrap, function(lmod){ fit_and_plot(lmod$mod.name, lmod$mod, load.only=T,plot.only.new=F, dataset = tms_data.nback)})
names(models_task.fitted) <- names(models_task)

loos_task=if.cached.load("loos_task",
                         invoke(loo_wrapper, .x = models_task.fitted, model_names = names(models_task.fitted)),
                         base=bname)

r2s_task <- lapply(models_task.fitted, bayes_R2, cl=22)
mod_task.weights = if.cached.load("mod_task.weights",
                                       map_df(c("loo", "waic", "stacking"), function(strat) {
                                         r = invoke(
                                           model_weights_wrapper,
                                           .x = models_task.fitted,
                                           weights = strat,
                                           model_names = names(models_task.fitted)
                                         )
                                         bind_cols(strategy = strat, data.frame(t(r)))
                                       }), bname)

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


## MCMC Intervals: 
color_scheme_set("viridisE")
mcmc_intervals(models_task.fitted$mod_task04, pars = c("b_zbv",
                                                       "b_zlog.apen",
                                                       "b_block_num",
                                                       "b_probeix",
                                                       "b_zbv:zlog.apen",
                                                       #"b_visit2",
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
                                       #"b_visit2" = "Visit",
                                       "b_conditionsham_rhTMS" = "Sham rhTMS",
                                       "b_conditionsham_arrhTMS" = "Sham arrhTMS",
                                       "b_conditionactive_rhTMS" = "Active rhTMS",
                                       "b_conditionactive_arrhTMS" = "Active arrhTMS"))
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
