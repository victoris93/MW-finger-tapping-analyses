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

tms_data_v1.path="/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/tms_analyses/data/G2/Visit_1"
tms_data_v2.path="/Users/VictoriaShevchenko/Documents/STAGE_M2/Analyses/tms_analyses/data/G2/Visit_2"
randomization <- read.csv("randomization.csv")

tms_data_files_v1 =list.files(path=tms_data_v1.path, full.names = T )
tms_data_files_v2 =list.files(path=tms_data_v2.path, full.names = T )

tms_data_v1 <- do.call(rbind,
                          lapply(tms_data_files_v1, function(fname){
                            read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                          }))

tms_data_v2 <- do.call(rbind,
                       lapply(tms_data_files_v2, function(fname){
                         read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                       }))


tms_data_v1$response[tms_data_v1$response == "s"] <- "lctrl"
tms_data_v1$response[tms_data_v1$response == "l"] <- "rctrl"

tms_data_v2$response[tms_data_v2$response == "s"] <- "lctrl"
tms_data_v2$response[tms_data_v2$response == "l"] <- "rctrl"


get.nback <- function(d, nback=20, which.apen.m=3, on.task.crit=1){
  if( !("condition" %in% names(d))){
    d$condition=1 ## this is for pilot2
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
    mutate(log.apen=log(log(2)-apen)) %>% left_join(randomization,by="subj") %>% ungroup -> d.nback.apen
  
  
  d.nback %>% group_by(subj,condition,focus,probeix) %>%
    mutate(timediff=c(0,diff(time))) %>% filter(timediff>0) %>%
    summarise(bv=sd(timediff)) %>%
    left_join(randomization,by=c("subj")) %>% ungroup -> d.nback.bv
  
  d.nback.apen %>% 
    full_join(d.nback.bv) %>% ungroup -> d.nback.full
  
  d.nback.full %>% ungroup %>% mutate(zlog.apen=-(log.apen-mean(log.apen, na.rm=T))/sd(log.apen,na.rm=T),
                                      zbv=(bv-mean(bv,na.rm=T))/sd(bv,na.rm=T)) %>%
    mutate(focus=factor(focus)) %>%
    mutate(off.focus=as.integer(focus=="off-task"),
           randomization=fct_relevel(randomization,"rh")) -> d.nback.full.z
  
  
  
  return(d.nback.full.z)
}

tms_data_v1$visit <- 1
tms_data_v2$visit <- 2

tms_data_v1.nback <- get.nback(tms_data_v1, nback = 25, which.apen.m=2)
tms_data_v2.nback <- get.nback(tms_data_v2, nback = 25, which.apen.m=2)

tms_data_v1.nback$visit <- 1
tms_data_v2.nback$visit <- 2

tms_data <- rbind(tms_data_v1, tms_data_v2)
tms_data.nback <- rbind(tms_data_v1.nback, tms_data_v2.nback)

#PROBE PROCESSING: 3 PROBES
probe_data_tms <- tms_data %>% 
  filter(stimulus %in% c("probe_task", "probe_intention", "probe_somnolence")) %>% 
  mutate(intention = case_when(stimulus == "probe_task" ~ lead(response)), 
         somnolence = case_when(stimulus == "probe_task" ~ lead(response, 2)))


probe_data_tms <- probe_data_tms %>% 
  filter(stimulus == "probe_task") %>% 
  arrange(visit, condition, response)

tms_data.nback <- tms_data.nback %>% 
  arrange(visit,condition, probe.response) %>%
  mutate(intention = probe_data_tms$intention,
         somnolence = probe_data_tms$somnolence,
         block_num = probe_data_tms$block_num,
         age = probe_data_tms$age,
         sex =probe_data_tms$sex)

tms_data.nback$intention <- as.integer(tms_data.nback$intention) + 1
tms_data.nback$intention_factor[tms_data.nback$intention >= 3] <- "intentional"
tms_data.nback$intention_factor[tms_data.nback$intention <= 2] <- "spontaneous"
tms_data.nback$sex[tms_data.nback$sex == FALSE] <- 1
tms_data.nback$sex[tms_data.nback$sex == "M"] <- 0
tms_data.nback$intention_factor <- as.factor(tms_data.nback$intention_factor)
tms_data.nback$sex <- as.factor(tms_data.nback$sex)
tms_data.nback$visit <- as.factor(tms_data.nback$visit)
tms_data.nback$somnolence <- as.integer(tms_data.nback$somnolence) + 1
tms_data.nback <- tms_data.nback %>%
  mutate(condition = as.factor(case_when(str_detect(condition, "baseline") ~ "baseline",
                               str_detect(condition, "active_rhTMS") ~ "active_rhTMS",
                               str_detect(condition, "sham_rhTMS") ~ "sham_rhTMS",
                               str_detect(condition, "active_arrhTMS") ~ "active_arrhTMS",
                               str_detect(condition, "sham_arrhTMS") ~ "sham_arrhTMS",
                               TRUE ~ as.character(condition))))
tms_data.nback$condition <- relevel(tms_data.nback$condition, ref = "baseline")

#Save data
save(tms_data.nback, file = "tms_data_preprocessed.Rdata")

# PLOTS

## MW distribution

tms_data.nback %>%
  ggplot(aes(x = probe.response)) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~ subj)


tms_data.nback %>%
  ggplot(aes(x = apen)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~ subj)

tms_data.nback %>%
  ggplot(aes(x = bv)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~ subj)


## MW distribution by intention

tms_data.nback %>%
  ggplot(aes(x = probe.response)) +
  geom_histogram(binwidth = 0.5) + 
  facet_wrap(~ intention_factor)

## Cor(AE x MW)

tms_data.nback %>%
  ggplot(aes(x = apen, y = probe.response)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~condition)

cor.test(x = tms_data.nback$apen,
    y = tms_data.nback$probe.response, "two.sided", method = "spearman")

## Cor(BV x MW)

tms_data.nback %>%
  ggplot(aes(x = bv, y = probe.response)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~condition)

cor.test(x = tms_data.nback$bv,
         y = tms_data.nback$probe.response, "two.sided", method = "spearman")

### Non-parametric ANOVA (Kruskal-Wallis rank sum test)

shapiro.test(tms_data.nback$probe.response)
shapiro.test(tms_data.nback$zlog.apen)
shapiro.test(tms_data.nback$bv)

kruskal.test(probe.response ~ condition, data = tms_data.nback) # groups are different
tms_data.nback$apen_ranks <- rank(tms_data.nback$probe.response)
by(tms_data.nback$apen_ranks, tms_data.nback$condition, mean) # where the difference lies
pgirmess::kruskalmc(probe.response ~ condition, data = tms_data.nback, cont = "two-tailed")
pgirmess::kruskalmc(probe.response ~ condition, data = tms_data.nback)

wilcox.test(filter(tms_data.nback_z, `Condition` == "sham_arrhTMS")$`MW Score`,
            filter(tms_data.nback_z, `Condition` == "active_rhTMS")$`MW Score`,
            paired = FALSE)

kruskal.test(zlog.apen ~ condition, data = tms_data.nback) # groups are different
tms_data.nback$apen_ranks <- rank(tms_data.nback$zlog.apen)
by(tms_data.nback$apen_ranks, tms_data.nback$condition, mean) # where the difference lies
pgirmess::kruskalmc(zlog.apen ~ condition, data = tms_data.nback)

kruskal.test(zbv ~ condition, data = tms_data.nback) # groups are different
tms_data.nback$bv_ranks <- rank(tms_data.nback$zbv)
by(tms_data.nback$bv_ranks, tms_data.nback$condition, mean) # where the difference lies
pgirmess::kruskalmc(zbv ~ condition, data = tms_data.nback)


### LONG DATA FRAME
tms_data.nback_z <- tms_data.nback %>%
  reshape2::melt(id.vars = c( "focus", "probe.response", "intention", "somnolence", "condition", "visit", "subj"), measure.vars = c("zlog.apen", "zbv"), varnames = c("Variable", "Score"))

data.table::setnames(tms_data.nback_z, old = c("focus", "probe.response", "intention", "somnolence", "condition", "visit", "subj",'variable', "value"), new = c('Focus', "MW Score", "Intention", "Alertness", "Condition", "Visit", "Subject",'Measure', "Score"))

tms_data.nback_z %>%
  ggplot(aes(x= `Focus`, y =`Score`,  group = `Measure`, color=`Measure`)) + 
  geom_pointrange(aes(shape = `Condition`),stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red")) +
  facet_wrap(~ `Condition`)

 tms_data.nback_z %>%
  ggplot(aes(x= `Condition`, y =`Score`, group = `Measure`, shape = `Condition`, color = `Measure`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("#F8766D", "#00BFC4")) +
  theme_set(theme_bw())

tms_data.nback_z %>%
  ggplot(aes(x= `MW Score`, y =`Score`, group = `Measure`, color = `Measure`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  facet_wrap(~`Condition`) +
  scale_color_manual(labels = c("AE", "BV"), values = c("#F8766D", "#00BFC4")) +
  facet_wrap(~`Subject`)

tms_data.nback_z %>% filter(`Subject` != "polya") %>%
  ggplot(aes(x= `Focus`, y =`Score`,  group = `Measure`, color=`Measure`)) + 
  geom_pointrange(stat="summary", fun.data=mean_se, fun.args = list(mult=1), position=position_dodge(0.05)) +
  geom_line(stat="summary", fun.data=mean_se, fun.args = list(mult=2)) +
  scale_color_manual(labels = c("AE", "BV"), values = c("blue", "red")) 
##############################################################################################


