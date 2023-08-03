library(tidyverse)
library(patchwork)

df10 <- read_csv("Binzer_2016_z10.csv")
df100 <- read_csv("Binzer_2016_z100.csv")
df10_2 <-read_csv("Binzer_2016_z10_finer.csv") 
df100_2 <-read_csv("Binzer_2016_z100_finer.csv") 


col = rev(RColorBrewer::brewer.pal(9, "Greens"))

p1 <- ggplot(df10, aes(x = eutrophication, y = temp, fill = persistence))+
  geom_tile()+
  scale_fill_gradientn(colours = col)

p2 <- ggplot(df100, aes(x = eutrophication, y = temp, fill = persistence))+
  geom_tile()+
  scale_fill_gradientn(colours = col)

p3 <- ggplot(df10_2, aes(x = eutrophication, y = temp-273.5, fill = persistence))+
  geom_tile()+
  scale_fill_gradientn(colours = col)

p3
