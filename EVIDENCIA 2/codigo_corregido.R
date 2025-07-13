# Cargar librerias
library(ade4)
library(ape)
library(adegenet)
library(Biostrings)
library(ggplot2)

# Define las variantes y sus nombres
variantes <- c("Civet.fasta", "HCoV-229E.fna", "HCoV-HKU1.fna", "HCoV-NL63.fna", "HCoV-OC43.fna", "MERS-CoV.fna", "RaTG13.fasta", "RpYN06.fasta", "SARS-CoV.fna", "SARS-CoV-2.fna")
nombre_variantes <- c("Civet SARS-CoV", "Human Coronavirus 229E", "Human Coronavirus HKU1", "Human Coronavirus NL63", "Human Coronavirus OC43", "Middle East Respiratory Syndrome", "Bat Coronavirus RaTG13", "Betacoronavirus RpYN06", "Severe Acute Respirstory Syndrome", "Severe Acute Respirstory Syndrome 2")

# Función para contar bases de ADN
contar_bases <- function(ADN) {
  a_count <- str_count(ADN, pattern = "A")
  t_count <- str_count(ADN, pattern = "T")
  g_count <- str_count(ADN, pattern = "G")
  c_count <- str_count(ADN, pattern = "C")
  n_count <- str_count(ADN, pattern = "N")
  data.frame(Base = c("A", "T", "G", "C", "N"), Count = c(a_count, t_count, g_count, c_count, n_count))
}

# Crear una lista para almacenar los datos y nombres de las variantes
datos <- list()
nombres <- list()

# Main
for (i in 1:length(variantes)) {
  # Imprime el nombre de la variante
  cat("Nombre de la variante:", nombre_variantes[i], "\n")
  
  # Lee el archivo de ADN y cuenta las bases
  ADN_set <- readDNAStringSet(variantes[i])
  adn <- toString(ADN_set)
  base_counts <- contar_bases(adn)
  
  # Agregar datos y nombres a las listas
  datos[[i]] <- base_counts
  nombres[[i]] <- nombre_variantes[i]
}

# Combinar datos en un único dataframe
datos_combinados <- do.call(rbind, datos)

# Añadir información de la variante
datos_combinados$Variante <- rep(nombre_variantes, each = 5)

# Crear el gráfico combinado
ggplot(datos_combinados, aes(x = Base, y = Count, fill = Base)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Variante, scales = "free_y") +
  labs(title = "Bases de ADN en diferentes variantes",
       x = "Base", y = "Count") +
  theme_minimal()

sessionInfo()

dna <- fasta2DNAbin(file="usflu.fasta")
# veamos que tiene nuestra variable dna
dna

annot <- read.csv("usflu.annot.csv", header=TRUE, row.names=1)
annot

D <- dist.dna(dna, model = "TN93")
length(D)
mean_dist <- mean(D, na.rm = TRUE)
D[is.na(D)] <- mean_dist

MatrizDG <- as.data.frame(as.matrix(D))
table.paint(MatrizDG, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)

tre <- nj(D)  # Use njs instead of nj
plot(tre, cex = 0.6)
title("Árbol de tipo NJ")

h_cluster <- hclust(D, method = "average", members = NULL)
plot(h_cluster, cex = 0.6)

# Define function to impute missing values with mean distance
impute_mean <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_x
  return(x)
}

myBoots <- boot.phylo(tre, dna, function(e) root(nj(impute_mean(dist.dna(e, model = "TN93"))), 1))

myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
plot(tre, edge.width = 1, cex = 0.7)
title("NJ tree + bootstrap values")
tiplabels(frame = "none", pch = 20, col = transp(num2col(annot$year, col.pal = myPal), 0.7), cex = 3, fg = "transparent")
temp <- pretty(1993:2008, 5)
legend("bottom", fill = transp(num2col(temp, col.pal = myPal), 0.7), leg = temp, ncol = 2, cex = 0.6)
nodelabels(myBoots, cex = 0.6)