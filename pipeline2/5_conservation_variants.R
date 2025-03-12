library(tidyr)
library(dplyr)
library(readr)
library(plotrix)
library(ggplot2)

#fichier issu du script combiner_les_var_pour_graphe.py qui contient tous les variants sur les ORF
setwd("path/to/input")
file <- "all_combined.csv"

#############lecture du fichier
df <- read.csv(file, header = TRUE, sep = ",")
#Passages en vraie liste
df$Passage <- lapply(df$Passage, function(x) {
  x_clean <- gsub("\\[|\\]|'", "", x)
  passages_list <- strsplit(x_clean, ",\\s*")[[1]]
  return(passages_list)
})

#separer cold/hot
variants_cold <- df %>% 
  filter(condition == "cold")
variants_hot <- df %>% 
  filter(condition == "hot")

#fonction pour ajouter une "cle" pour chaque variant different selon {pos, ORF, svtype, alt}
combine_variants <- function(df) {
  df <- df %>%
    group_by(pos, svtype, alt, ORF) %>%
    mutate(variant_cle = cur_group_id()) %>%
    ungroup()
  df_combined <- df %>%
    group_by(variant_cle) %>%
    summarise(
      pos = first(pos),
      svtype = first(svtype),
      alt = first(alt),
      ORF = first(ORF),
      Passage = list(unique(unlist(Passage)))
    ) %>%
    ungroup()
  return(df_combined)
}

variants_cold_combined <- combine_variants(variants_cold)
variants_hot_combined <- combine_variants(variants_hot)

##########graph cold: 

variants_cold_long <- variants_cold_combined %>%
  unnest(Passage)

labels_cold <- variants_cold_long %>%
  group_by(variant_cle) %>%
  filter(Passage == max(Passage)) %>%  # Get last Passage for labeling
  ungroup() %>%
  mutate(label = paste0(ORF))

###
pdf(file='conservation_cold.pdf',bg = 'transparent',width=12,height=10)
ggplot(variants_cold_long, aes(x = Passage, y = variant_cle, group = variant_cle)) +
  geom_line(data = subset(variants_cold_long, svtype == "INS"), 
            aes(x = Passage, y = variant_cle, group = variant_cle), 
            color = "deeppink4", size = 1, show.legend = FALSE) +
  geom_line(data = subset(variants_cold_long, svtype == "DEL"), 
            aes(x = Passage, y = variant_cle, group = variant_cle), 
            color = "darkblue", size = 1, show.legend = FALSE) +
  geom_point(aes(color = svtype), size = 2) + 
  geom_text(data = labels_cold, aes(label = label), hjust = -0.1, vjust = 0.1, size = 1.9, color = "black") +
  scale_color_manual(
    values = c("INS" = "deeppink2", "DEL" = "blue"),
    name = "Type"
  ) +
  labs(title = "Conservation des variants sur les ORF au cours des passages en condition froid",
       x = "Passages",
       y = " ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



##########graph hot: 

variants_hot_long <- variants_hot_combined %>%
  unnest(Passage)

labels_hot <- variants_hot_long %>%
  group_by(variant_cle) %>%
  filter(Passage == max(Passage)) %>%  # Get last Passage for labeling
  ungroup() %>%
  mutate(label = paste0(ORF))

###
pdf(file='conservation_hot.pdf',bg = 'transparent',width=12,height=10)
ggplot(variants_hot_long, aes(x = Passage, y = variant_cle, group = variant_cle)) +
  geom_line(data = subset(variants_hot_long, svtype == "INS"), 
            aes(x = Passage, y = variant_cle, group = variant_cle), 
            color = "deeppink4", size = 1, show.legend = FALSE) +
  geom_line(data = subset(variants_hot_long, svtype == "DEL"), 
            aes(x = Passage, y = variant_cle, group = variant_cle), 
            color = "darkblue", size = 1, show.legend = FALSE) +
  geom_point(aes(color = svtype), size = 2) + 
  geom_text(data = labels_hot, aes(label = label), hjust = -0.1, vjust = 0.1, size = 1.9, color = "black") +
  scale_color_manual(
    values = c("INS" = "deeppink2", "DEL" = "blue"),
    name = "Type"
  ) +
  labs(title = "Conservation des variants sur les ORF au cours des passages en condition chaud",
       x = "Passages",
       y = " ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


