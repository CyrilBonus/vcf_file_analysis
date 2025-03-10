library(tidyr)
library(dplyr)
library(readr)
library(plotrix)
library(ggplot2)


#######################################lecture des fichiers

#répertoire contenant tous les fichiers de stats générés par le script 3_infos_gen_variants.py
directory <- "/path/to/stats_dir"

#lecture des fichiers (csv)
csv_files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
variant_counts_list <- list()

#boucler sur chaque fichier et extraire les counts de INS et DEL
for (file in csv_files) { 
  data <- read.csv(file)
  file_name <- gsub(".csv", "", basename(file)) #extraire le nom du fichier
  generation <- as.integer(gsub(".*?(\\d+).*", "\\1", file_name)) #extraire la generation a partir du filename
  condition <- ifelse(grepl("hot", file_name, ignore.case = TRUE), "hot",
                      ifelse(grepl("cold", file_name, ignore.case = TRUE), "cold", 
                             NA)
                      )
  gen_cond <- paste0(generation, "_", condition)
  ins_count <- as.integer(data[data$Metric == "INS", "Value"])
  del_count <- as.integer(data[data$Metric == "DEL", "Value"])
  #.append les informations
  variant_counts_list[[length(variant_counts_list) + 1]] <- data.frame(File = file_name, 
                                                                       Generation = generation, 
                                                                       Condition = condition, 
                                                                       Gen_Cond = gen_cond,
                                                                       Variant = "INS", 
                                                                       Count = ins_count)
  variant_counts_list[[length(variant_counts_list) + 1]] <- data.frame(File = file_name, 
                                                                       Generation = generation, 
                                                                       Condition = condition, 
                                                                       Gen_Cond = gen_cond,
                                                                       Variant = "DEL", 
                                                                       Count = del_count)
}

variant_counts_df <- bind_rows(variant_counts_list) #combiner les resultats
print(variant_counts_df)


#######################################plot 

plot1<- ggplot(variant_counts_df, aes(x = Gen_Cond, y = Count, fill = Variant)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "nombre de variants de type INS et DEL par Generation et Condition",
       x = "Generation et condition",
       y = "Count") +
  scale_fill_manual(values = c("INS" = "lightskyblue", "DEL" = "lightpink")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot1)
#######################################tests de significativité

#ins del sep : test if diff significative between hot/cold for ins and del

# separer INS et DEL
INS_df <- variant_counts_df %>%
  filter(Variant == "INS") %>%
  select(Generation, Condition, Count) %>%
  pivot_wider(names_from = Condition, values_from = Count, names_prefix = "")
 
DEL_df <- variant_counts_df %>%
  filter(Variant == "DEL") %>%
  select(Generation, Condition, Count) %>%
  pivot_wider(names_from = Condition, values_from = Count, names_prefix = "")

#df2 : sans P15, pas de comparaison de temperature
INS_df2<-INS_df[ !(INS_df$Generation %in% c(15)), ]
DEL_df2<-DEL_df[ !(DEL_df$Generation %in% c(15)), ]

t.test(INS_df2$cold, INS_df2$hot, paired = FALSE)
t.test(DEL_df2$cold, DEL_df2$hot, paired = FALSE)
t.test(INS_df2$cold, DEL_df2$cold, paired = TRUE)
t.test(INS_df2$hot, DEL_df2$hot, paired = TRUE)
t.test(x=c(INS_df$cold, INS_df$hot), y=c(DEL_df$cold, DEL_df$hot), paired = TRUE)


a1 <- INS_df2$cold
a2 <- INS_df2$hot
b1 <- DEL_df2$cold
b2 <- DEL_df2$hot
 
comparaisons <- list(
  list(a1, a2, TRUE, "INS_df2$cold", "INS_df2$hot"),
  list(b1, b2, TRUE, "DEL_df2$cold", "DEL_df2$hot"),
  list(a1, b1, FALSE, "INS_df2$cold", "DEL_df2$cold"),
  list(a2, b2, FALSE, "INS_df2$hot", "DEL_df2$hot"),
)

for (comp in comparaisons) {
  data1 <- comp[[1]]
  data2 <- comp[[2]]
  ispaired <- comp[[3]] 
  nom1 <- comp[[4]]
  nom2 <- comp[[5]]
  test <- t.test(data1, data2, paired = ispaired)
  if (test$p.value < 0.05) {
    print(paste("Il y a une différence significative entre", nom1, "et", nom2, ", p-value = ", test$p.value))
  } else {
    print(paste("Pas de différence significative entre", nom1, "et", nom2, ", p-value = ", test$p.value))
  }
}

##########################


