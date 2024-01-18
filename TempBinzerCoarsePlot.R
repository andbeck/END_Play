library(tidyverse)
library(patchwork)

df10 <- read_csv("Binzer_2016_z10.csv")
df100 <- read_csv("Binzer_2016_z100.csv")
df10_2 <-read_csv("Binzer_2016_z10_finer.csv") 
df100_2 <-read_csv("Binzer_2016_z100_finer.csv") 

df10_2K <- read_csv("Binzer_2016_z10_coarseK.csv")
df100_apb <- read_csv("Binzer_2016_z100_APBAgain.csv")

df100_expbasimp <- read_csv("Binzer_2016_Z100_ExpBASimple.csv")

col = rev(RColorBrewer::brewer.pal(9, "Greens"))

basic <- df10_2K %>% 
  group_by(temp, eutrophication) %>% 
  summarise(
    meanPersistence = mean(persistence),
    sePersistence = sd(persistence)/sqrt(n())
  ) %>% 
  mutate(persistence2  = 1 - meanPersistence)

basic2 <- df100_expbasimp %>% 
  group_by(temp, eutrophication) %>% 
  summarise(
    meanPersistence = mean(persistence),
    sePersistence = sd(persistence)/sqrt(n())
  ) %>% 
  mutate(meanPersistence2 = 1-meanPersistence)


p0 <- ggplot(basic, aes(x = temp-273.15, y = meanPersistence, 
                        group = eutrophication, col = factor(eutrophication)))+
  geom_point()+
  geom_line()+
  labs(x = "Temperature (˚C)", y = "Mean Persistence (n = 30 webs)", 
       col = "K")+
  geom_errorbar(aes(ymin = (meanPersistence) - sePersistence,
                    ymax = (meanPersistence) + sePersistence), alpha = 0.5)+
  theme_bw(base_size = 15)

p0

ggplot(basic2, aes(x = eutrophication, y = temp-273.15, z = meanPersistence, fill = meanPersistence))+
  geom_tile()+
  ggtitle("100")+
  scale_fill_gradientn(colours = col)

## ----------

p1 <- ggplot(df10, aes(x = eutrophication, y = temp, fill = persistence))+
  geom_tile()+ggtitle("10")+
  scale_fill_gradientn(colours = col)

p2 <- ggplot(df100, aes(x = eutrophication, y = temp, fill = persistence))+
  geom_tile()+ggtitle("100")+
  scale_fill_gradientn(colours = col)

p1+p2

p3 <- ggplot(df10_2, aes(x = eutrophication, y = temp-273.5, fill = persistence))+
  geom_raster(interpolate = TRUE)+
  labs(y = "Temperature C˚", x = "Productivity", title = "10")+
  scale_fill_gradientn(colours = col)

p4 <- ggplot(df100_2, aes(x = eutrophication, y = temp-273.5, fill = persistence))+
  geom_raster(interpolate = TRUE)+
  labs(y = "Temperature C˚", x = "Productivity", title = "100")+
  scale_fill_gradientn(colours = col)

p3 + p4

