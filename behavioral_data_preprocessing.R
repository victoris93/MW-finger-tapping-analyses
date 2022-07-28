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

setwd(getwd())
sourceCpp("apen.cpp")

data.path="data/behavioral"

data_files=list.files(path=data.path, full.names = T )
data <- do.call(rbind,
                     lapply(data_files, function(fname){
                       read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                     }))

data <- data %>%
  filter(stimulus != "probe_content")

data$response[data$response == "s"] <- "lctrl"
data$response[data$response == "l"] <- "rctrl"

get.nback <- function(d, nback=20, which.apen.m=3, on.task.crit=1){
  if( !("condition" %in% names(d))){
    d$condition=1 
  }
  d %>% filter(stimulus == "probe_task") %>% 
    mutate(attention=if_else(as.integer(response)>on.task.crit, "on-task", "off-task")) %>%
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
                 stimulus != "probe_somnolence" &
                 stimulus != "pulse") %>%
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


#data_probes_long <- data %>% 
#  reshape2::reshape2::melt(id.vars = c("focus", "intention", "subj"), measure.vars = c("zlog.apen", "zbv"), varnames = c("Variable", "Z-score"))

data.nback <- get.nback(data, nback = 25, which.apen.m=2)
probe_data <- data %>% 
  filter(stimulus %in% c("probe_task", "probe_intention", "probe_somnolence")) %>% 
  mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)), 
         somnolence = case_when(stimulus == "probe_task" ~ lead(response, 2)))

probe_data <- probe_data %>% 
  filter(stimulus == "probe_task") %>% 
  arrange(condition, response)

data.nback <- data.nback %>% 
  arrange(condition, probe.response) %>%
  mutate(intention = probe_data$intention,
         somnolence = probe_data$somnolence,
         block_num = probe_data$block_num,
         age = probe_data$age,
         sex =probe_data$sex)

data.nback$intention <- as.integer(data.nback$intention) + 1
data.nback$intention_factor[data.nback$intention >= 3] <- "intentional"
data.nback$intention_factor[data.nback$intention <= 2] <- "spontaneous"
data.nback$sex[data.nback$sex == FALSE] <- 1
data.nback$sex[data.nback$sex == "M"] <- 0
data.nback$intention_factor <- as.factor(data.nback$intention_factor)
data.nback$sex <- as.factor(data.nback$sex)

data.nback <- data.nback %>%
  mutate(condition = case_when(str_detect(condition, "baseline") ~ "baseline",
                               str_detect(condition, "active_rhTMS") ~ "active_rhTMS",
                               str_detect(condition, "sham_rhTMS") ~ "sham_rhTMS",
                               str_detect(condition, "active_arrhTMS") ~ "active_arrhTMS",
                               str_detect(condition, "sham_arrhTMS") ~ "sham_arrhTMS",
                               TRUE ~ as.character(condition)))

####################################################### SANITY CHECK


data.nback_z <- data.nback %>%
  reshape2::melt(id.vars = c("focus", "subj"), measure.vars = c("zlog.apen", "zbv"))

data.table::setnames(data.nback_z, old = c("focus", "subj",'variable', "value"), new = c('Focus', "Subject",'Measure', "Z-score"))

data.nback_z %>%
  ggplot(aes(x= `Focus`, y =`Z-score`,  group = `Measure`, color=`Measure`)) + 
  geom_pointrange( stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red"))

data.nback_z %>%
  ggplot(aes(x= `Focus`, y =`Z-score`,  group = `Measure`, color=`Measure`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red")) +
  facet_wrap(~ `Subject`)

data.nback_z %>% filter(`Subject` != "polya") %>%
  ggplot(aes(x= `Focus`, y =`Z-score`,  group = `Measure`, color=`Measure`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red")) 


