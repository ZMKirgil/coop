library(haven)
library(readr)
library(tidyverse)
library(ggplot2)
library(nFactors)
library(GPArotation)
library(psych)

# Read files -----------------------------------------------------

coop <- read_sav("data/01_coop-data.sav")
save(coop, file="data/02_coop-data.R")

# Recoding & Scaling -----------------------------------------------------

z_std <- function(observed) {
  
  result <- (observed - mean(observed)) / sd(observed)
}

coop_clean <- coop %>%
  mutate(min = sec/60) %>%
  select(-(1:3), -c(35:37)) %>%
  rename(gid = gID,
         pid = pID) %>%
  zap_labels() %>%
  mutate_if(is.character,as.numeric) %>%
  mutate(mean_cont = (C1 + C2 + C3 + C4 + C5 + C6) / 6) %>% 
  mutate(ID1 = recode(ID1, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         ID6 = recode(ID6, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         ID7 = recode(ID7, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         ID8 = recode(ID8, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         Strat1 = recode(Strat1, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         Strat4 = recode(Strat4, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         App1 = recode(App1, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         App3 = recode(App3, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         App5 = recode(App5, `1` = "7", `2` = "6", `3` = "5", `4` = "4", `5` = "3", `6` = "2", `7` = "1"),
         TL = ifelse(condition %% 2, "1", "0"), # 0 = Loose (even), 1 = Tight (odd); 1 "TightSym", 2 "LooseSym",3 "TightMon",4 "LooseMon"
         MS = ifelse(condition > 2, "1", "0")) %>%
  mutate_if(is.character,as.numeric) %>%
  select_if(is.numeric) %>% 
  mutate_all(funs(z = z_std(.))) %>%
  select(-(32:33), -(41:49), -(64:65), -(70:80)) #checkout error message

save(coop_clean, file="data/04_coop-data-cleaned.R")

# Factor Analysis -----------------------------------------------------
# SID
ID <- coop_clean %>%
  select(c(39:46))

ev <- eigen(cor(ID, use = "complete.obs")) # get eigenvalues
ap <- parallel(subject = nrow(ID),var = ncol(ID),
               rep = 100,cent = .05)
nS <- nScree(x = ev$values, aparallel = ap$eigen$qevpea); plotnScree(nS) # -> suggests 1 factors
fa <- fa(ID, nfactors = 1, rotate = "promax", fm = "ML", scores = TRUE) ; fa

I<-factor.scores(ID, fa, Phi = NULL, method = c("Bartlett"),rho = NULL,impute = "none")
coop_clean$ID <- I[["scores"]] #extract Bartlett factor score
coop_clean$SID <- coop_clean$ID[,1] 

psych::alpha(ID[c("ID1_z","ID2_z","ID3_z","ID4_z","ID5_z", "ID6_z", "ID7_z", "ID8_z")], check.keys=TRUE) # -> alpha=0.91

# Collective Intentionality
CI <- coop_clean %>%
  select(c(53:55))

ev <- eigen(cor(CI, use="complete.obs")) # get eigenvalues
ap <- parallel(subject=nrow(CI),var=ncol(CI),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea); plotnScree(nS) 

fa<-fa(CI, nfactors=1, rotate="promax", fm="ML", scores = TRUE); fa

CI<-factor.scores(CI, fa, Phi = NULL, method = c("Bartlett"),rho=NULL,impute="none")
coop_clean$CI<-CI[["scores"]]
coop_clean$CI <- coop_clean$CI[,1]

alpha(coop_clean[c("Strat1_z", "Strat2_z", "Strat3_z")], check.keys=TRUE) #0.81

# Social Appropriateness
AP <- coop_clean %>%
  select(c(47:51))

psych::alpha(AP, check.keys=TRUE) # 0.47

# Missing Data -------------------------------------------------------------------------

coop_clean[!complete.cases(coop_clean),] 
mice::md.pattern(coop_clean) # -> two missing data points in Pre1a and Pre1ab


# Long Data Format ---------------------------------------------------------------------
coop_clean_long <- coop_clean %>%
  select(gid, pid, 4:9, 1:3, 24:25, 47:51, 30:32, 36:38, 47:52, 56, 58:59) %>%
  rename(sanction = Strat4_z) %>%
  gather(time, coop, C1:C6, factor_key=TRUE) %>%
  mutate(time = recode(time, `C1` = "1", `C2` = "2", `C3` = "3", `C4` = "4", `C5` = "5", `C6` = "6"),
         time = as.numeric(time))
  
save(coop_clean, file="data/05_coop-data-clean-long.R")

# Preparation for dynamic latent class analysis ------------------------------------------
coop_long_LG <- coop_clean_long %>%  # Time dummie
  mutate(timeD = recode(time, `1` = "0", `2` = "1", `3` = "1", `4` = "1", `5` = "1", `6` = "1"))

coop_long_LG$lag <- c(NA, coop_long_LG$coop[1:(length(coop_long_LG$coop)-1)])  # Lag variable 

library(stats); require(graphics) # Lag by group
lagg <- function(x)c(NA,x[1:(length(x)-1)])
coop_long_LG$laggrp <- ave(coop_long_LG$coop, coop_long_LG$gid, FUN = lagg)

write_sav(coop_long_LG, "coop_long_LG.sav")
