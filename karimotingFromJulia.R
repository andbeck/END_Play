library(tidyverse)
library(patchwork)

df <- read_csv("karimoto_basics.csv")

df

C_effect <- df %>% 
  filter(B == 0)
B_effect <- df %>% 
  filter(C == 0)

p1 <- ggplot(C_effect, aes(x = C, y = LoreauSync))+
  geom_line()

p2 <- ggplot(B_effect, aes(x = B, y = LoreauSync))+
  geom_line()

p1+p2
