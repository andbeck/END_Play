library(tidyverse)
library(patchwork)

df <- read_csv("tempRun.csv") %>% 

  mutate(fw = factor(fw)) %>% 
  mutate(temp = temp - 273.15)

sumT <- df %>% 
  group_by(temp) %>% 
  summarise(
    meanRichness = mean(richness),
    seRichness = 1.96*(sd(richness)/sqrt(n())),
    meanBiomass = mean(biomass),
    seBiomass = 1.96*(sd(biomass)/sqrt(n()))
  )

p1 <- ggplot(df, aes(x = temp, y = richness, col = fw))+
  geom_point()+
  geom_line()+
  geom_line(data = sumT, aes(x = temp, y = meanRichness), col = "black")+
  geom_ribbon(data = sumT, aes(x = temp, y = meanRichness,
                               ymax = meanRichness+seRichness,
                               ymin = meanRichness-seRichness), col = "grey", alpha = 0.5)+
  theme_bw()+
  theme(legend.position = "none")

p2 <- ggplot(df, aes(x = temp, y = biomass, col = fw))+
  geom_point()+
  geom_line()+
  geom_line(data = sumT, aes(x = temp, y = meanBiomass), col = "black")+
  geom_ribbon(data = sumT, aes(x = temp, y = meanBiomass,
                               ymax = meanBiomass+seBiomass,
                               ymin = meanBiomass-seBiomass), col = "grey", alpha = 0.5)+
  theme_bw()+
  theme(legend.position = "none")

p1/p2
