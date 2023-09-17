library(tidyverse)
library(patchwork)

# reference 20/40
fix20 <- read_csv("temp20Cons") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)
fix40 <- read_csv("temp40Cons") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

# Lin is linear change
Lin <- read_csv("tempLinRun.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

# Seasonal change at 20, 30, 40
Season20 <- read_csv("tempSeason20.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

Season30 <- read_csv("tempSeason30.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

Season40 <- read_csv("tempSeason30.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

# Season + Linear
LinSeason <- read_csv("tempLinSeason.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

# Lin + Variation
LinVar <- read_csv("tempLinVar.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15) %>% 
  mutate(replicate = rep(1:10, each = 200))

# Combine
df <- bind_rows(Lin, Season20, Season30, Season40, LinSeason) %>% 
  mutate(tempSeq = rep(c("Lin",
                         "Season20", "Season30", "Season40",
                         "LinSeason"), 
                       each = 200))

# summarise
# note use of step rather than temperature as x-axis
# for linear it is 20 -> 40
# for season it is 21.5 -> 19.5 repeated for 10 cycles

sumT <- df %>% 
  group_by(tempSeq, step) %>% 
  summarise(
    meanRichness = mean(richness),
    seRichness = 1.96*(sd(richness)/sqrt(n())),
    meanBiomass = mean(biomass),
    seBiomass = 1.96*(sd(biomass)/sqrt(n()))
  )

# plots, faceted by tempSeq type.

# reference
ref_df <- bind_rows(fix20, fix40) %>% 
  mutate(fixedTemp = rep(c(20,40), each = 10)) # nwebs
p1 <- ggplot(ref_df, aes(x = factor(temp), y = richness, fill = factor(temp)))+
  geom_boxplot()+
  scale_fill_manual(values = c("20" = "blue", "40" = "red"))+
  theme_bw()+
  theme(legend.position = "none")

p2 <- ggplot(ref_df, aes(x = factor(temp), y = biomass, fill = factor(temp)))+
  geom_boxplot()+
  scale_fill_manual(values = c("20" = "blue", "40" = "red"))+
  theme_bw()+
  theme(legend.position = "none")

p1+p2

sumRef <- ref_df %>% 
  group_by(temp) %>% 
  summarise(
    meanRich = mean(richness),
    meanBiomass = mean(biomass)
  )

# treatments

p3 <- ggplot(df, aes(x = step, y = richness, col = fw))+
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

p3+p1

p4 <- ggplot(df, aes(x = step, y = biomass, col = fw))+
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

p4+p2



