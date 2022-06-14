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
setwd(getwd())

bname="tms_analyses"
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options) {
  uncache.all(base = bname)
}

dir.create("cache/vars", recursive = TRUE)
############################################### MODELING MW

models_task <- list( # models commented out didn't fit
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


names(models_task) <- sprintf("mod_task%02i", 0:(length(models_task)-1)) # assign names to models
models_task.wrap <- map2(names(models_task), models_task, ~ list(mod.name=.x, mod=.y)) # wrap them in a list where model names and models are distinct variables
models_task.fitted=lapply(models_task.wrap, function(lmod){ fit_and_plot(lmod$mod.name, 
                                                                         lmod$mod,
                                                                         load.only=T,
                                                                         plot.only.new=F,
                                                                         dataset = tms_data.nback)}) # fit models
names(models_task.fitted) <- names(models_task) # name the brmsfit objects. Check the result by typing models_task.fitted in console

loos_task=if.cached.load("loos_task",  # load LOO-diagnostics if already computed. Otherwise, compute LOO
                         invoke(loo_wrapper, 
                                .x = models_task.fitted,
                                model_names = names(models_task.fitted)),
                         base=bname) # NB! If |elpd_diff| is < 4, se_diff is unreliable. Any model can be selected.

lapply(models_task.fitted, bayes_R2, cl=22) # Return R^2 for each model

## MW: MCMC Intervals
color_scheme_set("viridisE")
mcmc_intervals(models_task.fitted$mod_task04, pars = c("b_zbv", #plotting 95% HDIs for these coefficients
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

############################################### MODELING TASK PERFORMANCE
# The same procedure as above

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
models_bv.fitted=lapply(models_bv.wrap, function(lmod){ fit_and_plot(lmod$mod.name, # Student-t distribution
                                                                     lmod$mod, 
                                                                     load.only=T,
                                                                     plot.only.new=F,
                                                                     dataset = tms_data.nback,
                                                                     family = student(link = "identity", link_sigma = "log", link_nu = "logm1"))})
names(models_bv.fitted) <- names(models_bv)
loos_bv=if.cached.load("loos_bv",
                       invoke(loo_wrapper,
                              .x = models_bv.fitted,
                              model_names = names(models_bv.fitted)),
                       base=bname)

lapply(models_bv.fitted, bayes_R2, cl=22)


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


################################### FITTING MODELS ON AE

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
models_apen.fitted=lapply(models_apen.wrap, function(lmod){ fit_and_plot(lmod$mod.name,
                                                                         lmod$mod,
                                                                         load.only=T,
                                                                         plot.only.new=F,
                                                                         dataset = tms_data.nback,
                                                                         family = student(link = "identity", link_sigma = "log", link_nu = "logm1"))})
names(models_apen.fitted) <- names(models_apen)

loos_apen=if.cached.load("loos_apen",
                         invoke(loo_wrapper, .x = models_apen.fitted, model_names = names(models_apen.fitted)),
                         base=bname)

lapply(models_apen.fitted, bayes_R2, cl=22)

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
