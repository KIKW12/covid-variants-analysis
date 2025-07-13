# Cargar librerias
library(stringr)
library(Biostrings)
library(ggplot2)

# Define las variantes y sus nombres
variantes <- c("wuhan.fna", "ALPHA.fasta", "beta.fasta", "GAMMA.fasta", "Omnicron.fasta")
nombre_variantes <- c("Coronavirus", "Alpha", "Beta", "Gamma", "Omicron")

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
