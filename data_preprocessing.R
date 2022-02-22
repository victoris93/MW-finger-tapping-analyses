library(knitr)
library(dplyr)
library(purrr)
library(Rcpp)
library(ggplot2)
library(readr)
library(forcats)
sourceCpp("apen.cpp")

data.path="/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/data"

data_files_baseline=list.files(path=data.path, pattern=".*_baseline_.*\\.csv", full.names = T )
data_files_rhTMS=list.files(path=data.path, pattern=".*_rhTMS_.*\\.csv", full.names = T )
data_files_randTMS=list.files(path=data.path, pattern=".*_randTMS_.*\\.csv", full.names = T )


baseline_v <- do.call(rbind,
                      lapply(data_files_baseline, function(fname){
                        read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                      }))
rhTMS <- do.call(rbind,
                 lapply(data_files_rhTMS, function(fname){
                   read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                 }))
randTMS <- do.call(rbind,
                   lapply(data_files_randTMS, function(fname){
                     read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                   }))


cbind(baseline_v, condition="baseline") %>% 
  rbind(cbind(rhTMS,condition="rhTMS")) %>%
  rbind(cbind(randTMS,condition="randTMS")) -> all_data

all_data$response[all_data$response == "s"] <- "lctrl"
all_data$response[all_data$response == "l"] <- "rctrl"

all_data <- all_data %>% mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)))
all_data <- all_data %>% mutate(distraction = case_when(stimulus == "probe_task" ~ lead(response, 2)))

all_data$intention[all_data$intention == 0] <- "intentional"
all_data$intention[all_data$intention == 1] <- "spontaneous"
all_data[20877, "intention"] <- "intentional"
all_data[20877, "distraction"] <- 0

intention_distraction_values <- all_data[,c("intention", "distraction")] %>% 
  filter(!is.na(intention) & !is.na(distraction))

all_data %>% filter(stimulus == "probe_task")

#all_data$response[all_data$response == 3] <- 4
#all_data$response[all_data$response == 2] <- 3
#all_data$response[all_data$response == 1] <- 2
#all_data$response[all_data$response == 0] <- 1

get.nback <- function(d, nback=20, which.apen.m=3, on.task.crit=2){
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

all_data_nback <- get.nback(all_data, nback = 20, which.apen.m=2)


#
# trial-wise re-arrangement of data.
# use as: 
#
# pilot1 %>% mutate(ISI=as.factor(ISI)) %>% group_by(subj, ISI) %>% 
#   do( rearrange.df(.) ) %>% ungroup %>% 
#   mutate(reltime=resp_time-stim_time) %>% data.frame -> pilot1.bytrial
#
#
rearrange.ftrngt.bytrial.df.onesubj <- function(df){
  stim_times=with(df, time[stimulus==1])
  resp_times=with(df, time[stimulus==0])
  resp_trial=unlist(lapply(resp_times, function(x){ which.min(abs(stim_times+.1-x))}))
  d<-data.frame(
    subj=df$subj[1],
    ISI=df$ISI[1],
    trial=1:length(stim_times),
    stim_time=stim_times,
    resp_time=NA,
    response=NA,
    nresponses=0
  )
  d<-within(d,{
    resp_time[resp_trial]=resp_times;
    response[resp_trial]=with(df, as.character(response[stimulus==0]));
    nresponses=unlist(lapply(1:dim(d)[1], function(i){sum(i==resp_trial)}));
  })
}
