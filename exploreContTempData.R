library(tidyverse)
library(patchwork)

# Data Sets ----------------------------------------

## reference 10/40 ----
fix10 <- read_csv("temp10Cons.csv") %>% 
  mutate(fw = factor(fw)) %>%
  mutate(temp = temp - 273.15)
fix40 <- read_csv("temp40Cons.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

## Lin is linear change ----
Lin <- read_csv("tempLinRun.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

## Seasonal change at 20, 30, 40 ----

Season10 <- read_csv("tempSeason10.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

Season25 <- read_csv("tempSeason25.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

Season40 <- read_csv("tempSeason40.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

## Season + Linear ----

LinSeason <- read_csv("tempLinSeason.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

## Lin + Variation ----

LinVar_n <- read_csv("tempLinVar_n.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  # fws = 25; 1000 =  50reps of 20 temps
  mutate(replicate = rep(1:25, each = 1000))

LinVar_lo <- read_csv("tempLinVar_lo.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  # fws = 25; 1000 =  50reps of 20 temps
  mutate(replicate = rep(1:25, each = 1000))

LinVar_hi <- read_csv("tempLinVar_hi.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  # fws = 25; 1000 =  50reps of 20 temps
  mutate(replicate = rep(1:25, each = 1000))

## Lin + Season + Variation ----

LinVarSeason_n <- read_csv("tempLinVarSeason_n.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  # fws = 25; 1000 =  50reps of 20 temps
  mutate(replicate = rep(1:25, each = 1000))

LinVarSeason_lo <- read_csv("tempLinVarSeason_lo.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  # fws = 25; 1000 =  50reps of 20 temps
  mutate(replicate = rep(1:25, each = 1000))

LinVarSeason_hi <- read_csv("tempLinVarSeason_hi.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  # fws = 25; 1000 =  50reps of 20 temps
  mutate(replicate = rep(1:25, each = 1000))


## Combine all without variation ----
# each = 25fws*20times
df_novar <- bind_rows(Lin, Season10, Season25, Season40, LinSeason) %>% 
  mutate(tempSeq = rep(c("Lin",
                         "Season10", "Season25", "Season40",
                         "LinSeason"), 
                       each = 500))

## Combine all with variation

# summarise ------------------------------------------

# note use of step rather than temperature as x-axis
# for linear it is 20 -> 40
# for season it is 21.5 -> 19.5 repeated for 10 cycles

## Reference Temperatures
ref_df <- bind_rows(fix10, fix40) %>% 
  mutate(fixedTemp = rep(c(10,40), each = 25)) # nwebs

## reference data to annotate sequence plots ----
sumRef <- ref_df %>% 
  group_by(temp) %>% 
  summarise(
    meanRich = mean(richness),
    meanBiomass = mean(biomass)
  )

## summarise all with non var to mean from among n networks ----
sum_novar <- df_novar %>% 
  group_by(tempSeq, step) %>% 
  summarise(
    meanRichness = mean(richness),
    seRichness = 1.96*(sd(richness)/sqrt(n())),
    meanBiomass = mean(biomass),
    seBiomass = 1.96*(sd(biomass)/sqrt(n()))
  )

## summarise random var data to have single mean among 10 webs per random run ----
sum_var <- LinVar %>% 
  group_by(replicate, step) %>% 
  summarise(meanRichness = mean(richness),
            meanBiomass = mean(biomass))

## summarise random var+Season data to have single mean among 10 webs per random run ----
sumT_varsSeason <- LinVarSeason %>% 
  group_by(replicate, step) %>% 
  summarise(meanRichness = mean(richness),
            meanBiomass = mean(biomass))

# plots, faceted by tempSeq type. -------------------------------

## reference plot ----
p1 <- ggplot(ref_df, aes(x = factor(temp), y = biomass, fill = factor(temp)))+
  geom_boxplot()+
  scale_fill_manual(values = c("10" = "blue", "40" = "red"))+
  theme_bw()+
  theme(legend.position = "none")

p2 <- ggplot(ref_df, aes(x = factor(temp), y = richness, fill = factor(temp)))+
  geom_boxplot()+
  scale_fill_manual(values = c("10" = "blue", "40" = "red"))+
  theme_bw()+
  theme(legend.position = "none")

p1+p2



# temperature sequeucen treatment plots ----

p3 <- ggplot(df, aes(x = step, y = biomass, col = fw))+
  geom_point()+
  geom_line()+
  geom_line(data = sumT, aes(x = step, y = meanBiomass), col = "black")+
  geom_ribbon(data = sumT, aes(x = step, y = meanBiomass,
                               ymax = meanBiomass+seBiomass,
                               ymin = meanBiomass-seBiomass), col = "grey80", alpha = 0.5)+
  # add constant (single temp) means
  # need to adjust just these colours
  geom_point(data = sumRef, aes(x = 20, y = meanBiomass, col = factor(meanBiomass)), size = 5)+
  geom_hline(data = sumRef, aes(yintercept = meanBiomass), linetype = 'dashed')+
  facet_wrap(~ tempSeq)+
  theme_bw()+
  theme(legend.position = "none")

p3+p1

p4 <- ggplot(df, aes(x = step, y = richness, col = fw))+
  geom_point()+
  geom_line()+
  # add mean + CI ribbon
  geom_line(data = sumT, aes(x = step, y = meanRichness), col = "black")+
  geom_ribbon(data = sumT, aes(x = step, y = meanRichness,
                               ymax = meanRichness+seRichness,
                               ymin = meanRichness-seRichness), col = "grey80", alpha = 0.5)+
  # add constant (single temp) means
  # need to adjust just these colours
  geom_point(data = sumRef, aes(x = 20, y = meanRich, col = factor(meanRich)), size = 5)+
  geom_hline(data = sumRef, aes(yintercept = meanRich), linetype = 'dashed')+
  # facet by sequence
  facet_wrap(~ tempSeq)+
  # theming
  theme_bw()+
  theme(legend.position = "none")

p4+p1

## Temp Seq Plots with stoch, so each line is a unique set of randomness ----

p5 <- ggplot(sumT_vars, aes(x = step, y = meanBiomass, colour = factor(replicate)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")

p6 <- ggplot(sumT_vars, aes(x = step, y = meanRichness, colour = factor(replicate)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")

p5+p6


## Temp Seq Plots with Season and stoch, so each line is a unique set of randomness ----

p7 <- ggplot(sumT_varsSeason, aes(x = step, y = meanBiomass, colour = factor(replicate)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")

p8 <- ggplot(sumT_varsSeason, aes(x = step, y = meanRichness, colour = factor(replicate)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")

p7+p8

# ALL Plots In ----
(p1|p2|p3|p4)/(p5|p6|p7|p8)


# Checking that Temp 20 in sequence == reference 20 mean ----
## TRUE ----
seq20_step1 <- Lin %>% 
  filter(temp == 20) %>%
  summarise(meanRichness = mean(richness),
            meanBiomass = mean(biomass))

seq20_step1
sumT %>% filter(step == 1 & tempSeq == "Lin") %>% 
  select(meanRichness, meanBiomass)
