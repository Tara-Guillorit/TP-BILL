# Charger les bibliothèques
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Charger le fichier dans un dataframe
df <- read.table("extract.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
colnames(df) <- c("ORF", "Coordonnees", "Start", "End", "Len", "Type")

# Convertir les colonnes numériques
df$Start <- as.numeric(df$Start)
df$End <- as.numeric(df$End)
df$Len <- as.numeric(df$Len)

# Extraire les coordonnées
df_coords <- df %>%
  mutate(Coordonnees = str_extract_all(Coordonnees, "\\d+")) %>%  # Extraire uniquement les nombres
  unnest(Coordonnees) %>%                                         # Convertir liste en colonnes
  mutate(Coordonnees = as.numeric(Coordonnees))                   # Convertir en numérique

# Séparer les types pour le graph
df_del <- df %>% filter(Type == "DEL")  # Suppressions
df_ins <- df %>% filter(Type == "INS")  # Insertions

# Tracer le graphique
ggplot() +
  # Tracer les segments des coordonnées (en BLEU)
  geom_segment(data = df_coords, aes(x = Coordonnees, xend = lead(Coordonnees),
                                     y = ORF, yend = ORF),
               color = "blue", size = 1, na.rm = TRUE) +

  # Tracer les DEL (lignes rouges entre Start et End)
  geom_segment(data = df_del, aes(x = Start, xend = End,
                                  y = ORF, yend = ORF),
               color = "red", size = 3) +

  # Tracer les INS (points verts sur Start)
  geom_point(data = df_ins, aes(x = Start, y = ORF),
             color = "green", size = 3) +

  # Ajouter des labels aux DEL (lignes rouges)
  geom_text(data = df_del, aes(x = (Start + End) / 2, y = ORF,
                               label = paste(ORF, Type, Len, sep = " , ")),
            color = "red", size = 4, vjust = -1) +  # Texte au-dessus des lignes

  # Ajouter des labels aux INS (points verts)
  geom_text(data = df_ins, aes(x = Start, y = ORF,
                               label = paste(ORF, Type, Len, sep = " , ")),
            color = "green", size = 4, vjust = -1) +  # Texte au-dessus des points

  # Améliorer l'affichage
  theme_classic() +
  labs(title = "Représentation des ORF et des mutations",
       x = "Position", y = "ORF") +
  theme(axis.text.y = element_text(size = 8))  # Taille du texte des ORFs

