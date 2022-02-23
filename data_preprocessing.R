library(dplyr)
library(purrr)
library(Rcpp)
library(ggplot2)
library(readr)
library(forcats)
library(reshape2)
library(data.table)
library(scales)
library(stringr)
theme_set(theme_bw())

sourceCpp("/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/apen.cpp")

data.path_guillaume="/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/data/guillaume"
data.path_corentin="/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/data/corentin"
data.path_toni="/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/data/toni"


data_files_toni=list.files(path=data.path_toni, full.names = T )
data_files_guillaume=list.files(path=data.path_guillaume, full.names = T )
data_files_corentin=list.files(path=data.path_corentin, full.names = T )



data_toni <- do.call(rbind,
                          lapply(data_files_toni, function(fname){
                            read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                          }))
data_guillaume <- do.call(rbind,
                           lapply(data_files_guillaume, function(fname){
                             read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                           }))

data_corentin <- do.call(rbind,
                          lapply(data_files_corentin, function(fname){
                            read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                          }))


data_guillaume$response[data_guillaume$response == "s"] <- "lctrl"
data_guillaume$response[data_guillaume$response == "l"] <- "rctrl"

data_toni$response[data_toni$response == "s"] <- "lctrl"
data_toni$response[data_toni$response == "l"] <- "rctrl"

data_corentin$response[data_corentin$response == "s"] <- "lctrl"
data_corentin$response[data_corentin$response == "l"] <- "rctrl"

#all_data <- all_data %>% mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)))
#all_data <- all_data %>% mutate(distraction = case_when(stimulus == "probe_task" ~ lead(response, 2)))

#all_data$intention[all_data$intention == 0] <- "intentional"
#all_data$intention[all_data$intention == 1] <- "spontaneous"

get.nback_guillaume <- function(d, nback=20, which.apen.m=3, on.task.crit=1){
  if( !("condition" %in% names(d))){
    d$condition=1 ## this is for pilot2
  }
  d %>% filter(stimulus == "probe_task") %>% 
    mutate(attention=if_else(as.integer(response)>=on.task.crit, "on-task", "off-task")) %>%
    group_by(subj,condition) %>% 
    do({
      df=.
      df %>% mutate(probeix=1:n()) -> df
      dd=d %>% filter(subj==df$subj[1], condition==df$condition[1])
      probeix=which(dd$stimulus=="probe_task")
      nprobes=length(probeix)
      dd$probeix=rep(c(1:nprobes, -1), diff(c(0,probeix,dim(dd)[1])))
      
      trials=map(1:nback, function(x) df$trial-x) %>% unlist
      dd %>%
        filter(trial %in% trials) %>%
        filter(stimulus != "probe_task" & 
                 stimulus != "probe_intention" & 
                 stimulus != "probe_content" &
                 stimulus != "probe_somnolence") %>%
        mutate(focus=map_chr(probeix, function(t){
          df$attention[t]}),
          probe.response=map_int(probeix, function(t){
            as.integer(df$response[t])
          })+1)
    }) %>% ungroup -> d.nback
  
  
  d.nback %<>% group_by(condition, subj) %>% 
    filter(stimulus=="tap") %>% mutate(tap=case_when(response=="rctrl" ~ 1,
                                                     response=="lctrl" ~ 0,
                                                     TRUE ~ -1)) %>% 
    filter(tap>=0) %>% mutate(tap=as.integer(tap)) %>% ungroup 
  
  
  d.nback %>% 
    group_by(condition, subj, focus, probeix,probe.response) %>%
    summarize(apen=apen_int(tap,3)[which.apen.m+1]) %>%
    mutate(log.apen=log(log(2)-apen)) %>% ungroup -> d.nback.apen
  
  
  d.nback %>% group_by(subj,condition,focus,probeix) %>%
    mutate(timediff=c(0,diff(time))) %>% filter(timediff>0) %>%
    summarise(bv=sd(timediff)) %>% ungroup -> d.nback.bv
  
  d.nback.apen %>% 
    full_join(d.nback.bv) %>% ungroup -> d.nback.full
  
  d.nback.full %>% ungroup %>% mutate(zlog.apen=-(log.apen-mean(log.apen, na.rm=T))/sd(log.apen,na.rm=T),
                                      zbv=(bv-mean(bv,na.rm=T))/sd(bv,na.rm=T)) %>%
    mutate(focus=factor(focus)) %>%
    mutate(off.focus=as.integer(focus=="off-task")) -> d.nback.full.z
  
  
  
  return(d.nback.full.z)
}

get.nback_toni <- function(d, nback=20, which.apen.m=3, on.task.crit=1){
  if( !("condition" %in% names(d))){
    d$condition=1 ## this is for pilot2
  }
  d %>% filter(stimulus=="probe_task") %>% 
    mutate(attention=if_else(as.integer(response)<=on.task.crit, "on-task", "off-task")) %>%
    group_by(subj,condition) %>% 
    do({
      df=.
      df %>% mutate(probeix=1:n()) -> df
      dd=d %>% filter(subj==df$subj[1], condition==df$condition[1])
      probeix=which(dd$stimulus=="probe_task")
      nprobes=length(probeix)
      dd$probeix=rep(c(1:nprobes, -1), diff(c(0,probeix,dim(dd)[1])))
      
      trials=map(1:nback, function(x) df$trial-x) %>% unlist
      dd %>%
        filter(trial %in% trials) %>%
        filter(stimulus != "probe_task" & stimulus != "probe_intention") %>%
        mutate(focus=map_chr(probeix, function(t){
          df$attention[t]}),
          probe.response=map_int(probeix, function(t){
            as.integer(df$response[t])
          })+1)
    }) %>% ungroup -> d.nback
  
  
  d.nback %<>% group_by(condition, subj) %>% 
    filter(stimulus=="tap") %>% mutate(tap=case_when(response=="rctrl" ~ 1,
                                                     response=="lctrl" ~ 0,
                                                     TRUE ~ -1)) %>% 
    filter(tap>=0) %>% mutate(tap=as.integer(tap)) %>% ungroup 
  
  
  d.nback %>% 
    group_by(condition, subj, focus, probeix,probe.response) %>%
    summarize(apen=apen_int(tap,3)[which.apen.m+1]) %>%
    mutate(log.apen=log(log(2)-apen)) %>% ungroup -> d.nback.apen
  
  
  d.nback %>% group_by(subj,condition,focus,probeix) %>%
    mutate(timediff=c(0,diff(time))) %>% filter(timediff>0) %>%
    summarise(bv=sd(timediff)) %>% ungroup -> d.nback.bv
  
  d.nback.apen %>% 
    full_join(d.nback.bv) %>% ungroup -> d.nback.full
  
  d.nback.full %>% ungroup %>% mutate(zlog.apen=-(log.apen-mean(log.apen, na.rm=T))/sd(log.apen,na.rm=T),
                                      zbv=(bv-mean(bv,na.rm=T))/sd(bv,na.rm=T)) %>%
    mutate(focus=factor(focus)) %>%
    mutate(off.focus=as.integer(focus=="off-task")) -> d.nback.full.z
  
  
  
  return(d.nback.full.z)
}

get.nback_corentin <- function(d, nback=20, which.apen.m=3, on.task.crit=2){
  if( !("condition" %in% names(d))){
    d$condition=1 ## this is for pilot2
  }
  d %>% filter(stimulus=="probe_task") %>% 
    mutate(attention=if_else(as.integer(response)>=on.task.crit, "on-task", "off-task")) %>%
    group_by(subj,condition) %>% 
    do({
      df=.
      df %>% mutate(probeix=1:n()) -> df
      dd=d %>% filter(subj==df$subj[1], condition==df$condition[1])
      probeix=which(dd$stimulus=="probe_task")
      nprobes=length(probeix)
      dd$probeix=rep(c(1:nprobes, -1), diff(c(0,probeix,dim(dd)[1])))
      
      trials=map(1:nback, function(x) df$trial-x) %>% unlist
      dd %>%
        filter(trial %in% trials) %>%
        filter(stimulus != "probe_task" & stimulus != "probe_intention" & stimulus != "probe_distraction") %>%
        mutate(focus=map_chr(probeix, function(t){
          df$attention[t]}),
          probe.response=map_int(probeix, function(t){
            as.integer(df$response[t])
          })+1)
    }) %>% ungroup -> d.nback
  
  
  d.nback %<>% group_by(condition, subj) %>% 
    filter(stimulus=="tap") %>% mutate(tap=case_when(response=="rctrl" ~ 1,
                                                     response=="lctrl" ~ 0,
                                                     TRUE ~ -1)) %>% 
    filter(tap>=0) %>% mutate(tap=as.integer(tap)) %>% ungroup 
  
  
  d.nback %>% 
    group_by(condition, subj, focus, probeix,probe.response) %>%
    summarize(apen=apen_int(tap,3)[which.apen.m+1]) %>%
    mutate(log.apen=log(log(2)-apen)) %>% ungroup -> d.nback.apen
  
  
  d.nback %>% group_by(subj,condition,focus,probeix) %>%
    mutate(timediff=c(0,diff(time))) %>% filter(timediff>0) %>%
    summarise(bv=sd(timediff)) %>% ungroup -> d.nback.bv
  
  d.nback.apen %>% 
    full_join(d.nback.bv) %>% ungroup -> d.nback.full
  
  d.nback.full %>% ungroup %>% mutate(zlog.apen=-(log.apen-mean(log.apen, na.rm=T))/sd(log.apen,na.rm=T),
                                      zbv=(bv-mean(bv,na.rm=T))/sd(bv,na.rm=T)) %>%
    mutate(focus=factor(focus)) %>%
    mutate(off.focus=as.integer(focus=="off-task")) -> d.nback.full.z
  
  
  
  return(d.nback.full.z)
}


data_guillaume_nback <- get.nback_guillaume(data_guillaume, nback = 25, which.apen.m=2)
data_toni_nback <- get.nback_toni(data_toni, nback = 25, which.apen.m=2)
data_corentin_nback <- get.nback_corentin(data_corentin, nback = 25, which.apen.m=2)


probe_data_guillaume <- data_guillaume %>% 
  filter(stimulus %in% c("probe_task", "probe_intention", "probe_content", "probe_somnolence")) %>% 
  mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)), 
         content = case_when(stimulus == "probe_task" ~ lead(response, 2)),
         somnolence = case_when(stimulus == "probe_task" ~ lead(response, 3)))

probe_data_guillaume$intention[probe_data_guillaume$intention <= 1] <- "spontaneous"
probe_data_guillaume$intention[probe_data_guillaume$intention >= 2] <- "intentional"

probe_data_guillaume <- probe_data_guillaume %>% 
  filter(stimulus == "probe_task") %>% 
  arrange(condition, response)

data_guillaume_nback <- data_guillaume_nback %>% 
  arrange(condition, probe.response) %>%
  mutate(intention = probe_data_guillaume$intention,
         content = probe_data_guillaume$content,
         somnolence = probe_data_guillaume$somnolence,
         block_num = probe_data_guillaume$block_num)

data_guillaume_nback$intention <- as.factor(data_guillaume_nback$intention)

#TONI: PROBE PROCESSING

data_toni_nback$probe.response[data_toni_nback$probe.response == 4] <- "Score: 1"
data_toni_nback$probe.response[data_toni_nback$probe.response == 1] <- "Score: 4"
data_toni_nback$probe.response[data_toni_nback$probe.response == 3] <- "Score: 2"
data_toni_nback$probe.response[data_toni_nback$probe.response == 2] <- "Score: 3"

data_toni_nback$probe.response[data_toni_nback$probe.response == "Score: 1"] <- 1
data_toni_nback$probe.response[data_toni_nback$probe.response == "Score: 4"] <- 4
data_toni_nback$probe.response[data_toni_nback$probe.response == "Score: 2"] <- 2
data_toni_nback$probe.response[data_toni_nback$probe.response == "Score: 3"] <- 3 

data_toni_nback$probe.response <- as.numeric(data_toni_nback$probe.response)

probe_data_toni <- data_toni %>% 
  filter(stimulus %in% c("probe_task", "probe_intention")) %>% 
  mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)))

probe_data_toni$intention[probe_data_toni$intention == 1] <- "no TUT"
probe_data_toni$intention[probe_data_toni$intention == 2] <- "spontaneous TUT"
probe_data_toni$intention[probe_data_toni$intention == 3] <- "intentional TUT"

probe_data_toni <- probe_data_toni %>% 
  filter(stimulus == "probe_task") %>% 
  arrange(condition,response)

data_toni_nback <- data_toni_nback %>% 
  arrange(condition, probe.response) %>%
  mutate(intention = probe_data_toni$intention)

data_toni_nback$intention <- as.factor(data_toni_nback$intention)
data_toni_nback$intention <- relevel(data_toni_nback$intention, "no TUT")


#CORENTIN: PROBE PROCESSING
probe_data_corentin <- data_corentin %>% 
  filter(stimulus %in% c("probe_task", "probe_intention", "probe_distraction")) %>% 
  mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)),
         distraction = case_when(stimulus == "probe_task" ~ lead(response, 2)))

probe_data_corentin$intention[probe_data_corentin$intention == 0] <- "intentional"
probe_data_corentin$intention[probe_data_corentin$intention == 1] <- "spontaneous"


probe_data_corentin <- probe_data_corentin %>% 
  filter(stimulus == "probe_task") %>% 
  arrange(condition, response)

data_corentin_nback <- data_corentin_nback %>% 
  arrange(condition, probe.response) %>%
  mutate(intention = probe_data_corentin$intention)

data_corentin_nback$intention <- as.factor(data_corentin_nback$intention)

# MERGE DATA

all_data <- bind_rows(data_corentin_nback, data_toni_nback, data_guillaume_nback)
all_data <- all_data %>%
  mutate(condition = case_when(str_detect(condition, "baseline") ~ "baseline",
                               str_detect(condition, "active_rhTMS") ~ "active_rhTMS",
                               str_detect(condition, "sham_rhTMS") ~ "sham_rhTMS",
                               str_detect(condition, "active_arrhTMS") ~ "active_arrhTMS",
                               str_detect(condition, "sham_arrhTMS") ~ "sham_arrhTMS",
                               TRUE ~ as.character(condition)))

# PLOTS

## MW distribution

all_data %>%
  ggplot(aes(x = probe.response)) +
  geom_histogram(binwidth = 0.5)

## MW distribution by intention

all_data %>%
  ggplot(aes(x = probe.response)) +
  geom_histogram(binwidth = 0.5) + 
  facet_wrap(~ intention)

## Cor(AE x MW)

all_data %>%
  ggplot(aes(x = apen, y = probe.response)) +
  geom_point() + 
  geom_smooth(method = "lm")

cor.test(x = all_data$apen,
    y = all_data$probe.response, "two.sided", method = "spearman")

## Cor(BV x MW)

all_data %>%
  ggplot(aes(x = bv, y = probe.response)) +
  geom_point() + 
  geom_smooth(method = "lm")

cor.test(x = all_data$bv,
         y = all_data$probe.response, "two.sided", method = "spearman")

## AE x BV interaction

all_data_z <- all_data %>% 
  reshape2::melt(id.vars = "focus", measure.vars = c("zlog.apen", "zbv"), varnames = c("Variable", "Z-score"))

setnames(all_data_z, old = c('focus','variable', "value"), new = c('Focus','Measure', "Z-score"))

all_data_z %>%
  ggplot(aes(x= `Focus`, y =`Z-score`,  group = `Measure`, color=`Measure`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red"))
#all_data_toni_nback <- arrange(probe.response

##############################################################################################


dataframe_TMS <- all_data_toni_nback %>% 
  filter(str_detect(condition, "TMS"))

df_ttest <- all_data_toni_nback %>% 
  filter(str_detect(condition, "baseline")) %>% 
  mutate(condition = "Baseline") %>% 
  bind_rows(dataframe_TMS)

df_ttest %>% 
  ggplot(aes(x = probe.response)) +
  geom_bar() +
  facet_wrap(~ condition)

df_ttest_long_response_condition <-  df_ttest %>% 
  reshape2::melt(id.vars = "condition", 
                 measure.vars = c("probe.response"))

setnames(df_ttest_long_response_condition, 
         old = c('condition','variable', "value"), 
         new = c('Condition','Measure', "MW Score"))

df_ttest_long_response_condition %>%
  filter(Condition == "Baseline"|!str_detect(Condition, "arrhTMS")) %>%
  ggplot(aes(x= `Condition`, y =`MW Score`, group = `Condition`, color=`Condition`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2))# +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red"))

df_ttest_long_response_condition %>%
  ggplot(aes(x= `Condition`, y =`MW Score`,  color=`Condition`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  theme_bw()

library(gplots)
df_ttest_long_response_condition %>%
  filter(Condition == "Baseline"|str_detect(Condition, "arrhTMS")) %>%
  plotmeans(`MW Score` ~ `Condition`, frame = FALSE,
          mean.labels = TRUE, connect = FALSE)
  #geom_boxplot(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  #geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2))# +
#scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red"))
