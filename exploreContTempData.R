library(tidyverse)
library(patchwork)

# Lin is linear change
Lin <- read_csv("tempLinRun.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

# Seasonal change at 20
Season <- read_csv("tempSeasonRun20.csv") %>% 
  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

# Combine
df <- bind_rows(Lin, Season) %>% 
  mutate(tempSeq = rep(c("Lin","Season"), each = 200))

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

p1 <- ggplot(df, aes(x = step, y = richness, col = fw))+
  geom_point()+
  geom_line()+
  geom_line(data = sumT, aes(x = step, y = meanRichness), col = "black")+
  geom_ribbon(data = sumT, aes(x = step, y = meanRichness,
                               ymax = meanRichness+seRichness,
                               ymin = meanRichness-seRichness), col = "grey", alpha = 0.5)+
  facet_wrap(~ tempSeq)+
  theme_bw()+
  theme(legend.position = "none")

p1

p2 <- ggplot(df, aes(x = step, y = biomass, col = fw))+
  geom_point()+
  geom_line()+
  geom_line(data = sumT, aes(x = step, y = meanBiomass), col = "black")+
  geom_ribbon(data = sumT, aes(x = step, y = meanBiomass,
                               ymax = meanBiomass+seBiomass,
                               ymin = meanBiomass-seBiomass), col = "grey", alpha = 0.5)+
  facet_wrap(~ tempSeq)+
  theme_bw()+
  theme(legend.position = "none")

p1/p2
