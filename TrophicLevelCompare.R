# Trophic Level Compare

library(tidyverse)
library(stringr)

TLfunc <- read_csv("TL.csv")
TL_temp <- t(TLfunc) %>% as.data.frame() %>% rownames_to_column()
node <- as.numeric(str_remove(TL_temp$rowname, "s"))
TL_use <- data.frame(TL_temp, node) %>% select(-rowname)
names(TL_use) <- c("TL","node")

d2pfunc <- read_csv("d2p.csv")
d2p_temp <- t(d2pfunc) %>% as.data.frame() %>% rownames_to_column()
node <- as.numeric(str_remove(d2p_temp$rowname, "s"))
d2p_use <- data.frame(d2p_temp, node) %>% select(-rowname)
names(d2p_use) <- c("TL","node")

TLeva <- read_csv("EvaTL.csv")
TLeva_use <- data.frame(TLeva, node = 1:50) 
names(TLeva_use) <- c("TL","node")

compareTL <- left_join(TL_use, d2p_use, by = "node") %>% 
  left_join(., TLeva_use, by = "node") %>% 
  select(node, TL_base = TL.x, TL_d2p = TL.y, TL_eva = TL) %>% 
  pivot_longer(-node, names_to = "TLMethod", values_to = "TL")


ggplot(compareTL, aes(x = node, y = TL, fill = TLMethod))+
  geom_col(position = "dodge")+
  theme_bw()+
  coord_flip()

           