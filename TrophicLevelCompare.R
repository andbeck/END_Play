# Trophic Level Compare

library(tidyverse)
library(stringr)
library(ATNr)

## functions Lurgi ----
NormalizeMatrix <- function(M){
  colsum_M <- colSums(M);
  colsum_M[colsum_M==0] <- 1;
  return(t(t(M)/colsum_M));
  
}

TrophicPositions <- function(FW){
  S <- dim(FW)[1];
  if(S < 3) return(FW);
  M <- NormalizeMatrix(FW);
  if(det(diag(S) - t(M)) != 0){
    TP <- solve(diag(S)-t(M), rep(1,S));
  }else{
    tmp <- diag(S);
    for(i in 1:9){
      tmp <- tmp %*% t(M) + diag(S);
    }
    TP <- tmp %*% rep(1,S);
  }
  W <- M;
  TLs <- TrophicLevels(FW);
  return(list(TP,TLs,W));
}


TrophicLevels <- function(FW){
  S <- dim(FW)[1];
  
  if(S < 2) return(FW)
  
  TrLevels <- rep(-1, S)
  
  #we find the neighbours of each node in the network
  M_temp <- FW;
  diag(M_temp) <- 0;
  
  #all the species with trophic position 1 are basal species, and therefore
  #belong to the trophic level 0
  TrLevels[FW == 1.0] <- 0;
  producers <- which(TrLevels == 0);
  for(i in 1:S){
    if(TrLevels[i] == -1){
      herb <- TRUE;
      top <- TRUE;
      
      #for each species we verify two things:
      #first: if any of its prey does not belong to trophic level 0,
      #then it is not a herbivore
      if(sum(TrLevels[which(M_temp[,i] != 0)]) != 0) herb <- FALSE;
      
      #second: if it has any predators, then it is not a top predator
      if(sum(M_temp[i,]) > 0) top <- FALSE;
      
      #after we've found out whether it is a herbivore or a top
      #predator, or none of those, assign the value accordingly
      if(herb){
        TrLevels[i] = 1;
      }else if(top){
        TrLevels[i] = 3;
      }else{
        TrLevels[i] = 2;
      }
    }
  }
  return(TrLevels);
}


WebJulia <- read_csv("FiftyWeb.csv") %>% as.data.frame() %>% 
  mutate(across(everything(), ~ as.numeric(.x))) %>% 
  as.matrix()

# remember to transpoise matrix as Julia web is row = predator
# R codes are row = resource

TLlurgi <- TrophicPositions(t(WebJulia))[[1]]

# ATNr ----

TLatn <- TroLev(t(WebJulia))

# cbind(TLlurgi, TLatn)
# cbind(TLeva$TL, TLatn)

# Julia Processed Data ====

TLfunc <- read_csv("TL.csv")
TL_temp <- t(TLfunc) %>% as.data.frame() %>% rownames_to_column()
node <- as.numeric(str_remove(TL_temp$rowname, "s"))
TL_use <- data.frame(TL_temp, node) %>% select(-rowname)
names(TL_use) <- c("TL_base","node")

d2pfunc <- read_csv("d2p.csv")
d2p_temp <- t(d2pfunc) %>% as.data.frame() %>% rownames_to_column()
node <- as.numeric(str_remove(d2p_temp$rowname, "s"))
d2p_use <- data.frame(d2p_temp, node) %>% select(-rowname)
names(d2p_use) <- c("TL_d2p","node")

TLeva <- read_csv("EvaTL.csv")
TLeva_use <- data.frame(TLeva, node = 1:50) 
names(TLeva_use) <- c("TL_eva","node")

TLlurgi_use <- data.frame(TLlurgi, node = 1:50) 
names(TLlurgi_use) <- c("TL_lurgi","node")

TLatn_use <- data.frame(TLatn, node = 1:50) %>% remove_rownames()
names(TLatn_use) <- c("TL_atn","node")

# comparison ----

compareTL <- left_join(TL_use, d2p_use, by = "node") %>% 
  left_join(., TLeva_use, by = "node") %>% 
  left_join(., TLlurgi_use, by = "node") %>% 
  left_join(., TLatn_use, by = "node") %>% 
  select(node, TL_base, TL_d2p, TL_eva, TL_lurgi, TL_atn) %>% 
  pivot_longer(-node, names_to = "TLMethod", values_to = "TL")


ggplot(compareTL, aes(x = node, y = TL, fill = TLMethod))+
  geom_col(position = "dodge")+
  theme_bw()+
  coord_flip()

           