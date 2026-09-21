################################################################################################################################
# Intrinsic credibility of statistically significant outcomes in Anaesthesia and Pain Medicine trials
# 
# Markus Huber*, Patrick Wüthrich and Alexander Fuchs
# Department of Anaesthesiology and Pain Medicine, Inselspital, Bern University Hospital, University of Bern, Bern, Switzerland.
# *Corresponding author: Markus Huber (markus.huber@insel.ch)

# Code Version: 21 Sep 2026
################################################################################################################################

#######
# setup
#######

rm(list=ls())
library(meta)
library(metafor)
library(tidyverse)
library(ggExtra)
library(compareGroups)
library(BSDA)

`%notin%` <- Negate(`%in%`)
options(scipen=9999)

#####################################################################
# Helper function for illustrative Figure 1
#
# the function is based on the
# Analysis of Credibility: https://doi.org/10.1098/rsos.171047
#
# input:  upper (U) and lower (L) bounds of a 95%-confidence interval
# output: a dataframe with the boundaries of:
#         (i)   a sceptical prior distribution
#         (ii)  the observed data
#         (iii) the posterior distribution
#####################################################################

intervals <- function(L,U){
  # the so-called sceptical limits (SL)
  SL        = (U-L)^2/(4*sqrt(U*L))
  # prior
  sd.prior  = SL/1.96
  var.prior = sd.prior^2
  mu.prior  = mean(c(-SL,SL))
  # data
  mu.data   = mean(c(U,L))
  sd.data   = abs(U-mu.data)/1.96
  var.data  = sd.data^2
  # posterior
  var.posterior  = 1/((1/var.data)+(1/var.prior))
  sd.posterior   = sqrt(var.posterior)
  mu.posterior   = var.posterior*((mu.data/var.data)+(mu.prior/var.prior))
  return(
    rbind(
      data.frame(type = "Sceptical prior",y = mu.prior,ymin = min(c(-SL,SL)), ymax = max(c(-SL,SL))),
      data.frame(type = "Data",y = mu.data,ymin = min(c(mu.data-1.96*sd.data,mu.data+1.96*sd.data)), ymax = max(c(mu.data-1.96*sd.data,mu.data+1.96*sd.data))),
      data.frame(type = "Posterior",y = mu.posterior,ymin = min(c(mu.posterior-1.96*sd.posterior,mu.posterior+1.96*sd.posterior)), ymax = max(c(mu.posterior-1.96*sd.posterior,mu.posterior+1.96*sd.posterior)))))
}

############################################################################
# load the data from the Cochrane Database of Systematic Reviews (CDSR)
# note the data from each Systematic Review can be downloaded (if available)
# directly at: https://www.cochranelibrary.com/cdsr/reviews

# number of reviews
n.reviews = length(list.files("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/data", recursive = FALSE, pattern=".rm5"))
  
# set directory
setwd("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/data")

# define data frames
df.results <- data.frame()
df.dens    <- data.frame()
df.sample  <- data.frame()
df.credi <- data.frame()

# read the data of each systematic review and meta analyses where available
for (i in 1:n.reviews){
  
  Sys.sleep(1)
  # these SR/MA did not provide outcome data / provided corrupted files
  if (i %notin% c(32,46,49,61,74,76,85,89,93,97,121,129,134,135,154,169,226,230,261,266,278,285,289,294,327,332,333)){
  print(i)
  
  # read the data with the dedicated "read.rm5" function in R
  rm(review)
  review = (read.rm5(file = list.files("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/data", recursive = FALSE, pattern=".rm5")[i],numbers.in.labels = FALSE))
  rm(cochrane.id);cochrane.id = list.files("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/data", recursive = FALSE, pattern=".rm5")[i] %>% str_replace("StatsDataOnly.rm5","")
  
  review <- review %>% as.data.frame()
  
  # get the the number of outcomes / analyses
  rm(structure)
  structure <- review %>% select(comp.no,outcome.no) %>% unique()
  
  # extract all data for each comparison ("comp.no") and outcome ("outcome.no")
  for (myrow in 1:nrow(structure)){
    
    #################
    # binary outcomes
    rm(data.studies)
    data.studies = review %>%
      filter(comp.no==structure$comp.no[myrow] & outcome.no==structure$outcome.no[myrow]) %>% 
      select(c("studlab", "event.e", "n.e", "event.c", "n.c","outclab")) %>% 
      unique() %>% 
      na.omit() %>% 
      mutate(
        ai = event.e,
        bi = n.e-event.e,
        ci = event.c,
        di = n.c-event.c
      )
    
    if (nrow(data.studies)>0){
      
      for (mystudy in 1:nrow(data.studies)){
        
        # for binary outcomes, derive the P values with Fisher's exact test
        rm(myfisher)
        myfisher <- fisher.test(matrix(c(data.studies[mystudy,"ai"], data.studies[mystudy,"bi"], data.studies[mystudy,"ci"], data.studies[mystudy,"di"]), nrow = 2, byrow = FALSE))
        
        df.credi <- rbind(
          df.credi,
          data.frame(
            meta = i,
            file = list.files("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/data", recursive = FALSE, pattern=".rm5")[i],
            comp.no     = structure$comp.no[myrow],
            outcome.no  = structure$outcome.no[myrow],
            study       = mystudy,
            studyid     = data.studies$studlab[mystudy],
            outcome     = data.studies$outclab[mystudy],
            type        = "binary",
            or.mean     = myfisher$estimate,
            or.lower    = myfisher$conf.int[1],
            logor.lower = myfisher$conf.int[1] %>% log(),
            or.upper    = myfisher$conf.int[2],
            logor.upper = myfisher$conf.int[2] %>% log(),
            # P values
            p           = myfisher$p.value,
            # t statistic based on the P value (see equations in https://doi.org/10.1098/rsos.181534))
            t     = qnorm(1-myfisher$p.value/2)) %>% 
            mutate(
              p.ic = 2*(1-pnorm(t/(sqrt(2)))) # P value for intrinsic credibility
            ) %>% 
            mutate(
              p.rep = 1-(p.ic/2) # Probability of replication (https://doi.org/10.1098/rsos.181534))
            ))
      }
    }
    
    ######################
    # continuous outcomes
    rm(data.studies)
    data.studies = review %>%
      filter(comp.no==structure$comp.no[myrow] & outcome.no==structure$outcome.no[myrow]) %>% 
      select(c("studlab", "mean.e", "sd.e","n.e","mean.c", "sd.c","n.c","outclab")) %>% 
      unique() %>% 
      na.omit()
    
    if (nrow(data.studies)>0){
      
      for (mystudy in 1:nrow(data.studies)){
        
        # Summarized t-test based on summary measures (from the package BSDA)
        rm(mytest)
        mytest <- with(data.studies[mystudy,],
                       tsum.test(
                         mean.x = mean.e,
                         s.x = sd.e,
                         n.x = n.e,
                         mean.y = mean.c,
                         s.y = sd.c,
                         n.y = n.c, 
                         var.equal = FALSE))
        
        df.credi <- rbind(
          df.credi,
          data.frame(
            meta       = i,
            file       = list.files("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/data", recursive = FALSE, pattern=".rm5")[i],
            comp.no    = structure$comp.no[myrow],
            outcome.no = structure$outcome.no[myrow],
            study      = mystudy,
            studyid    = data.studies$studlab[mystudy],
            outcome    = data.studies$outclab[mystudy],
            type       = "continuous",
            # set odds ratio to NA for continuous outcomes
            # (note these were only required for illustrative Figure 1)
            or.mean     = NA,
            or.lower    = NA,
            logor.lower = NA,
            or.upper    = NA,
            logor.upper = NA,
            # P values etc as above
            p     = mytest$p.value,
            t     = qnorm(1-mytest$p.value/2)) %>% 
            mutate(
              p.ic = 2*(1-pnorm(t/(sqrt(2))))
            ) %>% 
            mutate(
              p.rep = 1-(p.ic/2)
            ))
      }
    }    
    
  }
  }
}

########
########
# Figure
########
########

# (A) significant, but not instrinsically

rm(dummy)
dummy <- df.credi %>% 
  filter(type=="binary") %>% 
  filter(p>0.03 & p<0.04) %>% 
  filter(studyid=="POISE 2008" & comp.no==1 & outcome.no==2)

dummy$p

fig1A <- intervals(L=dummy$logor.lower,U=dummy$logor.upper) %>% 
  mutate(
    y=exp(y),
    ymin=exp(ymin),
    ymax=exp(ymax)
  ) %>% 
  mutate(type = factor(type,levels = c("Sceptical prior","Data","Posterior"))) %>% 
  ggplot(aes(x=type,y=y,ymin=ymin,ymax=ymax,color=type))+
  geom_hline(yintercept=1,linetype="dashed")+
  geom_point(size=4)+
  geom_errorbar(width=0.2)+
  theme_classic()+
  theme(
    text = element_text(size=16),
    legend.position = "none"
  )+
  scale_y_continuous(breaks = seq(0.2,2.2,by=0.2),limits = c(0.1,2.2))+
  ylab("Odds ratio (OR)")+
  ggsci::scale_color_jama()+
  xlab("")+
  annotate("text", x = "Data", y = 0.9, label = "P=0.031")+
  geom_segment(aes(x = "Posterior", y = 0.7, xend = "Data", yend = 0.7),arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "black",
               size = 1)+
  annotate("text", x = 2.5, y = 0.55, label = "Reverse Bayes")

# (B) significant and  instrinsically

rm(dummy)
dummy <- df.credi %>% 
  filter(type=="binary") %>% 
  filter(p>0.001 & p<0.002) %>% 
  filter(studyid=="Rumbak 2004" & comp.no==1 & outcome.no==3)

dummy$p

fig1B <- intervals(L=dummy$logor.lower,U=dummy$logor.upper) %>% 
  mutate(
    y=exp(y),
    ymin=exp(ymin),
    ymax=exp(ymax)
  ) %>% 
  mutate(type = factor(type,levels = c("Sceptical prior","Data","Posterior"))) %>% 
  ggplot(aes(x=type,y=y,ymin=ymin,ymax=ymax,color=type))+
  geom_hline(yintercept=1,linetype="dashed")+
  geom_point(size=4)+
  geom_errorbar(width=0.2)+
  theme_classic()+
  theme(
    text = element_text(size=16),
    legend.position = "none"
  )+
  scale_y_continuous(breaks = seq(0.2,2.2,by=0.2),limits = c(0.1,2.2))+
  ylab("Odds ratio (OR)")+
  ggsci::scale_color_jama()+
  xlab("")+
  annotate("text", x = "Data", y = 0.75, label = "P=0.002")+
  geom_segment(aes(x = "Posterior", y = 1.1, xend = "Data", yend = 1.1),arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "black",
               size = 1)+
  annotate("text", x = 2.5, y = 1.3, label = "Reverse Bayes")


figAB <- ggpubr::ggarrange(fig1A,fig1B,labels = c("A","B"),nrow=1)


#########
# overall
#########

fig1C <- df.credi %>% 
  pivot_longer(c(p,p.ic)) %>% 
  mutate(
    name.new = case_when(
      name=="p"~"Traditional",
      name=="p.ic"~"Related to intrinsic credibility"
    )
  ) %>% 
  mutate(
    name.new = factor(name.new,levels = c("Traditional","Related to intrinsic credibility"))) %>% 
  ggplot(aes(x = value, fill = name.new))+
  geom_histogram(alpha=0.3,position = "identity",breaks = seq(0,1,by=0.05))+
  aes(y = 2*after_stat(count)/sum(after_stat(count)))+
  scale_x_continuous(breaks = seq(0,1,by=.05),guide = guide_axis(n.dodge = 2),
                     labels = c("0","0.05","0.1","0.15","0.2","0.25",
                                "0.3","0.35","0.4","0.45","0.5",
                                "0.55","0.6","0.65","0.7","0.75",
                                "0.8","0.85","0.9","0.95","1"),
                     expand = c(0,0))+
  scale_y_continuous(labels = scales::percent,breaks = seq(0,0.3,by=0.04),expand = c(0,0),limits = c(0,0.3))+
  theme_classic()+
  theme(
    text = element_text(size=18),
    legend.position = c(0.5,0.7),
    legend.title = element_blank(),
    legend.text=element_text(size=18)
  )+
  geom_vline(xintercept=0.05, linetype="dashed",size=1)+
  ggsci::scale_color_aaas()+
  ggsci::scale_fill_aaas()+
  guides(fill = guide_legend(nrow = 1))+
  xlab(expression(italic(P)~value))+
  ylab("Distribution")

ggpubr::ggarrange(
  figAB,
  fig1C,
  labels = c("","C"),
  heights = c(1,1),
  nrow=2
)

ggsave("T:/Research KAS/M. Huber/INTERNAL CREDIBILITY/BJA/revisions/figure1.jpeg",dpi=600,width=8,height=7)

##################################
# Table for additional information

df.table <- 
  df.credi %>% 
  transmute(
    p.category = case_when(
      p<0.0056~"<=0.0056",
      (p>=0.0056 & p<=0.05)~"0.0056 - 0.05",
      TRUE~">0.05"))
                           
createTable(
  compareGroups(
    ~.,data=df.table,method=NA
  )
)
                           
createTable(
  compareGroups(
    ~.,data=df.credi,method=NA
  )
)

##################################################################################################
# These tables provide additional information on the percentage of intrinsically credible findings

df.table <- 
  df.credi %>% 
  transmute(
    p,
    p.ic,
    p.significance = case_when(p<=0.05~"Significant",TRUE~"Non significant"),
    p.ic.significance = case_when(p.ic<=0.05~"Significant",TRUE~"Non significant"))
    
createTable(
  compareGroups(
    ~.,data=df.table,method=NA
  )
)

createTable(
  compareGroups(
    ~.,data=df.table %>% 
      filter(p<=0.05),method=NA
  )
)


   