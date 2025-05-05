needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan", "readxl")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")

'%notin%' <- Negate('%in%')

dupdataR <- read_excel("dupdataR.xlsx") 

which_repeated <- dupdataR %>% group_by(Barcode) %>% summarize(num_entry = n()) %>% arrange(desc(num_entry))

dupdataR.t  <- dupdataR %>% filter(Barcode %in% which_repeated$Barcode[1:10])
dupdataR.tc <- dupdataR.t %>% group_by(Barcode) %>% slice(5:8) %>% mutate(Barcode = Barcode * 20)
dupdataR.tk <- dupdataR.t %>% group_by(Barcode) %>% slice(1:4)
dupdataR %<>% filter(Barcode %notin% which_repeated$Barcode[1:10])
dupdataR %<>% rbind(., dupdataR.tc, dupdataR.tk)

dupdataR %<>% 
  mutate(Barcode = as.factor(Barcode)) %>%
  mutate(Barcode = as.numeric(Barcode))

dupdataR %>% 
  filter(Swab == 2) %>%
  group_by(Barcode) %>% 
  summarize(
  lwr = min(estCopyPerSwab)
, upr = max(estCopyPerSwab)
) %>% mutate(dif = upr - lwr) %>% 
  mutate(equiv = ifelse(lwr == 0 & upr > 0, 1, 0)) %>%
  mutate(equiv = ifelse(lwr > 0 & upr > 0, 2, 0)) %>% 
  mutate(equiv = as.factor(equiv)) %>% {
  ggplot(., aes(x = Barcode)) + 
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = equiv)) +
    geom_point(data = dupdataR %>% 
        filter(BD == 1, Swab == 2)
      , aes(Barcode, estCopyPerSwab)) +
    guides(colour = F) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 5E6, 9E6)) +
      scale_colour_brewer(palette = "Dark2")
  }

dupdataR %>% group_by(Barcode, Swab) %>% 
  summarize(ms = mean(estCopyPerSwab)) %>%
  ungroup(Swab) %>% 
  pivot_wider(names_from = Swab, values_from = ms) %>%
  mutate(diff = abs(`2` - `1`)) %>% 
  ungroup() %>% 
  summarize(prop_equiv = length(which(diff > 0)) / n())

dupdataR %>% 
  filter(Swab == 2) %>%
  dplyr::select(Barcode, Swab, SwabRep, estCopyPerSwab) %>%
  group_by(Barcode) %>% 
  summarize(
    num_zero = length(which(round(estCopyPerSwab, 1) == 0))
  , num_pos  = length(which(round(estCopyPerSwab) > 0))
    ) %>% mutate(equiv = abs(num_zero - num_pos)) %>%
  mutate(equiv = ifelse(equiv == 3, 1, 0)) %>%
  ungroup() %>% 
  summarize(prop_equiv = 1 - length(which(equiv == 1)) / n())

dupdataR %>% 
  group_by(Barcode) %>% 
  summarize(
    si = (length(which(estCopyPerSwab > 0))/n())) %>% 
  left_join(.,
  dupdataR %>% 
  group_by(Barcode) %>% 
  filter(estCopyPerSwab > 0) %>%
  summarize(
    cs = mean(estCopyPerSwab)) 
    ) %>% mutate(spec_bar = ifelse(Barcode == 22, 1, 0) %>% as.factor()) %>% {
    ggplot(., aes(cs, si)) + geom_point(size = 3, aes(shape = spec_bar)) +
      guides(shape = F) +
      scale_x_continuous(
        trans = "pseudo_log"
     , breaks = c(1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 5E6, 9E6)) +
      xlab("Average copies per swab across all PCR runs") +
      ylab("Proportion of runs that found non-zero copies")
  }

dupdataR %>% 
  group_by(Barcode) %>% 
  filter(Swab == 2) %>%
  summarize(cs = mean(estCopyPerSwab)) %>%  
  left_join(.,
dupdataR %>% group_by(Barcode, Swab) %>% 
  summarize(ms = mean(estCopyPerSwab)) %>%
  ungroup(Swab) %>% 
  pivot_wider(names_from = Swab, values_from = ms) %>%
  mutate(diff = abs(`2` - `1`)) %>% 
  dplyr::select(Barcode, diff)
  ) %>% mutate(spec_bar = ifelse(Barcode == 22, 1, 0) %>% as.factor()) %>% {
      ggplot(., aes(cs, diff)) + geom_point(aes(shape = spec_bar), size = 3) +
      scale_x_continuous(trans = "pseudo_log", breaks = c(1E1, 1E2, 1E3, 1E4, 1E5, 1E6)) +
      scale_y_continuous(trans = "pseudo_log", breaks = c(1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 5E6, 9E6)) +
      xlab("Average copies per swab across all PCR replicates for swab 2") +
      ylab("Difference in copies per swab between 
swab 1 and the average of swab 2")
      }
