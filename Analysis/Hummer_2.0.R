setwd("C:/Users/jcriosor/OneDrive/Proyectos/Hummingbirds_RN/Analisis_2024")


########################### MORFOLOGIA ################################################

###### PARENTESIS ####
## Primero creo un subset de 100 arboles para poder correr el codigo rapido mientras
## queda el definitivo. Al final asegurarse de que estoy corriendo el arbol completo (1000 arboles) en
## alltrees <- read.nexus.

library(ape)

# Cargar el archivo de árboles filogenéticos en formato .nex
alltrees <- read.nexus("Phylo/Strisores_Humm_all_1000/Strisores_humm_all_1000_Hacket_all.nex")

# Crear un subset de 100 árboles
subset_trees <- alltrees[1:100]

# Guardar el subset de 100 árboles en un nuevo archivo .nex
write.nexus(subset_trees, file = "Phylo/Subset_100_trees.nex")


##### FIN DEL PARENTESIS ###

# Cargar las librerías necesarias
library(ape)
library(phytools)
library(readr)
library(phangorn)
library(geiger)

# Cargar los datos morfológicos
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", 
                        header = TRUE, sep = ";", fileEncoding = "UTF-8")

# Renombrar la columna de especies (Species3) para trabajar con ella
colnames(morpho_data)[which(colnames(morpho_data) == "Species3")] <- "Species"

# Reemplazar espacios en los nombres de las especies por "_"
morpho_data$Species <- gsub(" ", "_", morpho_data$Species)
rownames(morpho_data) <- morpho_data$Species

# Convertir comas a puntos en las medidas morfométricas
morpho_measures <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", 
                     "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance", 
                     "Secondary1", "Hand.Wing.Index", "Tail.Length", "Mass")

# Convertir las columnas de medidas a numérico
morpho_data[morpho_measures] <- lapply(morpho_data[morpho_measures], 
                                       function(x) as.numeric(gsub(",", ".", x)))

# Verificar si hay valores NA en las medidas
if(any(is.na(morpho_data[morpho_measures]))) {
  cat("Existen valores NA en las medidas morfológicas. Eliminando filas con NA...\n")
  morpho_data <- na.omit(morpho_data)
}

# Cargar el árbol filogenético (subset de 100 árboles)
subset_trees <- read.nexus("Phylo/Subset_100_trees.nex")

# Crear un árbol de maxima credibilidad a partir del subset de los 100 árboles
sum_tree_all <- maxCladeCred(subset_trees)

# Verificar las especies que están en el árbol y en los datos morfológicos
cat("Verificando coincidencias entre las especies en los datos y el árbol...\n")
name_check_result <- name.check(sum_tree_all, morpho_data)

# Verificar si name.check devolvió una lista (con especies faltantes o adicionales)
if (is.list(name_check_result)) {
  # Especies en los datos pero no en el árbol
  if (length(name_check_result$data_not_tree) > 0) {
    cat("Especies en los datos pero no en el árbol:\n")
    print(name_check_result$data_not_tree)
  }
  
  # Especies en el árbol pero no en los datos
  if (length(name_check_result$tree_not_data) > 0) {
    cat("Especies en el árbol pero no en los datos:\n")
    print(name_check_result$tree_not_data)
  }
} else {
  cat("No hay discrepancias en los nombres de las especies entre los datos y el árbol.\n")
}

# Filtrar el dataframe morfológico para que contenga solo las especies presentes en el árbol
morpho_data <- morpho_data[rownames(morpho_data) %in% sum_tree_all$tip.label, ]

# Reordenar los datos morfológicos para que coincidan con el orden de las puntas del árbol
morpho_data <- morpho_data[match(sum_tree_all$tip.label, rownames(morpho_data)), ]

# Verificar si hay valores NA nuevamente después del filtrado
cat("Verificando si hay valores NA después del filtrado...\n")
if(any(is.na(morpho_data[morpho_measures]))) {
  cat("Existen valores NA después del filtrado. Eliminando filas con NA...\n")
  morpho_data <- na.omit(morpho_data)
}

# Verificar si las dimensiones coinciden entre los datos morfológicos y las puntas del árbol
cat("Número de especies en los datos morfológicos:", nrow(morpho_data), "\n")
cat("Número de especies en el árbol:", length(sum_tree_all$tip.label), "\n")

if (nrow(morpho_data) != length(sum_tree_all$tip.label)) {
  stop("El número de especies en los datos morfológicos no coincide con el número de especies en el árbol.")
}

# Realizar el PCA filogenético utilizando el modelo de Movimiento Browniano (BM)
pPCA_all <- phyl.pca(sum_tree_all, morpho_data[, morpho_measures], method = "BM", mode = "corr")

# Extraer los scores del PCA
PCAs_scores <- as.data.frame(pPCA_all$S)

# Agregar la columna de familia (Family3) a los scores
PCAs_scores$family <- morpho_data$Family3

# Guardar los scores en un archivo CSV
write.csv(PCAs_scores, "PCA_Scores.csv")

# Mostrar los resultados del PCA
print(summary(pPCA_all))

# Graficar los resultados del PCA filogenético
biplot(pPCA_all)

# Grafico no pude hacr funcionar la leyenda pero es solo indicativo
# Si no tienen GGally
# install.packages("GGally")
GGally::ggpairs(PCAs_scores, columns = 1:4,
                mapping = ggplot2::aes(color = family),
                upper  = list(continuous = "blank"),
                diag  = list(continuous = "blankDiag"))+
  theme_minimal()


### esta es otra forma de hacer una figura similar, pero esta si deja ver la leyenda
### anque es menos bonita.

# Cargar las librerías necesarias
library(ggplot2)

# Asegurarse de que la columna 'family' sea un factor
PCAs_scores$family <- as.factor(PCAs_scores$family)

# Definir una paleta de colores
color_palette <- rainbow(length(unique(PCAs_scores$family)))
names(color_palette) <- levels(PCAs_scores$family)

# Ajustar la disposición del gráfico para tener espacio para la leyenda
layout(matrix(1:2, ncol = 2), widths = c(3, 1))  # Reserva espacio para la leyenda en la segunda columna

# Generar el gráfico de pares (PC1 a PC4) utilizando pairs de R base
par(mar = c(4, 4, 2, 2))  # Ajuste de márgenes para el gráfico
pairs(PCAs_scores[, 1:4], 
      col = color_palette[PCAs_scores$family],  # Colorear por familia
      pch = 19,  # Tipo de punto
      main = "PCA Filogenético de Aves (PC1 a PC4)")

# Crear la leyenda en la segunda columna reservada por layout()
par(mar = c(0, 0, 0, 0))  # Sin márgenes para la leyenda
plot.new()  # Crear un nuevo gráfico vacío para la leyenda
legend("center", 
       legend = levels(PCAs_scores$family), 
       col = color_palette, 
       pch = 19, 
       title = "Familia", 
       bty = "n")  # Sin borde en la leyenda




### ahora una hermosa figura en 3D


# Definir una paleta de colores personalizada usando rainbow
color_palette <- rainbow(length(unique(PCAs_scores$family)))

# Crear un gráfico 3D con los primeros tres componentes principales
fig <- plot_ly(PCAs_scores, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~family,  # Colorear por familia
               colors = color_palette,  # Paleta personalizada
               type = "scatter3d",  # Tipo de gráfico 3D
               mode = "markers",  # Usar marcadores
               marker = list(size = 5))  # Tamaño de los puntos

# Añadir título y etiquetas a los ejes
fig <- fig %>%
  layout(title = "PCA Filogenético 3D (PC1, PC2, PC3)",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")),
         legend = list(title = list(text = 'Familia')))

# Mostrar el gráfico 3D
fig

### la figura muestra que Trochilidae se separa mucho en el espacio morfologico, no se
### mezcla con nada mas. Igual que podargidae, extrañamente.


##### Ahora quiero ver como se distribuye el morfoespacio pero solo en colibries
##### discriminando por los nueve clados descritos en McGuire 2014

# Cargar las librerías necesarias
library(ape)
library(phytools)
library(readr)
library(GGally)
library(ggplot2)
library(plotly)

# Cargar los datos morfológicos, incluidos los clades
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", 
                        header = TRUE, sep = ";", fileEncoding = "UTF-8")

# Renombrar la columna de especies (Species3) para trabajar con ella
colnames(morpho_data)[which(colnames(morpho_data) == "Species3")] <- "Species"

# Reemplazar espacios en los nombres de las especies por "_"
morpho_data$Species <- gsub(" ", "_", morpho_data$Species)
rownames(morpho_data) <- morpho_data$Species

# Filtrar los resultados del pPCA para incluir solo Trochilidae
trochilidae_species <- morpho_data[morpho_data$Family3 == "Trochilidae", "Species"]
PCAs_scores_trochilidae <- PCAs_scores[rownames(PCAs_scores) %in% trochilidae_species, ]

# Agregar la columna 'clade' a los resultados filtrados
PCAs_scores_trochilidae$clade <- morpho_data[rownames(PCAs_scores_trochilidae), "Clade"]

# Crear una paleta de colores basada en los clades únicos
clades_unicos <- unique(PCAs_scores_trochilidae$clade)
color_palette <- setNames(rainbow(length(clades_unicos)), clades_unicos)  # Asignar colores a cada clade

### Gráfico 2D con ggpairs para los primeros 4 componentes principales
ggpairs(PCAs_scores_trochilidae, columns = 1:4,  # Usar solo los primeros 4 componentes principales
        mapping = aes(color = clade),  # Colorear por clade
        upper = list(continuous = "blank"),  # Dejar la parte superior en blanco
        diag = list(continuous = "blankDiag")) +  # Dejar la diagonal en blanco
  scale_color_manual(values = color_palette) +  # Usar la misma paleta de colores
  theme_minimal() +
  labs(title = "PCA Filogenético de Trochilidae (PC1 a PC4)") +
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título

### Gráfico 3D utilizando plotly (PC1 vs PC2 vs PC3)
fig <- plot_ly(PCAs_scores_trochilidae, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~clade,  # Colorear por clado
               colors = color_palette,  # Usar la misma paleta de colores
               type = "scatter3d",  # Tipo de gráfico 3D
               mode = "markers",  # Usar marcadores
               marker = list(size = 5))  # Tamaño de los puntos

fig <- fig %>%
  layout(title = "PCA Filogenético 3D de Trochilidae (PC1, PC2, PC3)",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")),
         legend = list(title = list(text = 'Clado')))

# Mostrar el gráfico 3D
fig

### entre los linajes de colibries hay mucha siperposición, no es tan facil separar los
### grupos en el espacio morfologico.




### ####### ahora vamos a calcular las densidades de cada variable en el morfoespacio por familia

# Cargar las librerías necesarias
library(ggplot2)
library(reshape2)

# Cargar los datos morfológicos
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", 
                        header = TRUE, sep = ";", fileEncoding = "UTF-8")

# Renombrar la columna de especies (Species3) para trabajar con ella
colnames(morpho_data)[which(colnames(morpho_data) == "Species3")] <- "Species"

# Medidas morfométricas de interés
morpho_measures <- c("Beak.Depth", "Beak.Length_Culmen", "Beak.Length_Nares", 
                     "Beak.Width", "Hand.Wing.Index", "Kipps.Distance", 
                     "Secondary1", "Tail.Length", "Tarsus.Length", "Wing.Length")

# Convertir las columnas morfométricas a formato numérico
morpho_data[morpho_measures] <- lapply(morpho_data[morpho_measures], 
                                       function(x) as.numeric(gsub(",", ".", x)))

# Ver cuántos valores NA hay por cada variable morfométrica
cat("Número de valores faltantes (NA) por columna:\n")
print(colSums(is.na(morpho_data[morpho_measures])))

# Si una columna tiene muchos valores faltantes, podemos eliminarla
# También podemos eliminar filas que contengan muchos NA, según sea necesario
# En este caso, eliminamos filas que contengan NA en más de la mitad de las medidas morfométricas

# Establecer un umbral para decidir cuándo eliminar filas (en este caso, más de la mitad de las columnas)
umbral <- length(morpho_measures) / 2

# Filtrar filas con más de umbral de NA
morpho_data_filtrado <- morpho_data[rowSums(is.na(morpho_data[morpho_measures])) <= umbral, ]

# Verificar nuevamente después del filtrado
cat("Número de filas después de aplicar el umbral de NA:\n")
print(nrow(morpho_data_filtrado))

# Convertir los datos a formato largo para ggplot2
data_long <- melt(morpho_data_filtrado, id.vars = "Family3", measure.vars = morpho_measures)

# Verificar si el proceso de melt produjo valores en la columna 'variable'
cat("Valores únicos en la columna 'variable' después del melt:\n")
print(unique(data_long$variable))

# Crear el gráfico de densidades para cada variable morfométrica, separado por familia
if (length(unique(data_long$variable)) > 0) {
  ggplot(data_long, aes(x = value, fill = Family3, color = Family3)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~variable, scales = "free", ncol = 4) +  # Crear un panel por variable
    scale_x_log10() +  # Escala logarítmica en el eje X
    theme_minimal() +
    labs(title = "Densidad estimada de medidas morfométricas por familia (Escala logarítmica)", 
         x = "Valor", 
         y = "Densidad", 
         fill = "Familia", 
         color = "Familia") +
    theme(
      plot.title = element_text(hjust = 0.5),  # Centrar el título
      legend.position = "bottom",              # Posicionar la leyenda en la parte inferior
      legend.title = element_text(size = 12),  # Tamaño del título de la leyenda
      legend.text = element_text(size = 10),   # Tamaño del texto de la leyenda
      strip.text = element_text(size = 12))    # Tamaño del texto en los paneles
} else {
  cat("No hay suficientes valores en la columna 'variable' para facetar.\n")
}

















#######################  Codigo no exitoso para hacer el phylomorphospace  #################

PCAs_scores <- as.data.frame(pPCA_all$S)
PCAs_scores$family <- morpho_data$Family3
PCAs_scores$species <- rownames(PCAs_scores)


# Extraer las puntuaciones de PC1 y PC2 para las especies
pca_matrix <- as.matrix(PCAs_scores[, 1:2])

# Usar phylomorphospace para obtener las coordenadas de las relaciones filogenéticas
phylo_coords <- phylomorphospace(sum_tree_all, pca_matrix, plot = FALSE)

# Crear un dataframe con las coordenadas de las conexiones filogenéticas
coords_edges <- data.frame(
  xstart = phylo_coords$X[sum_tree_all$edge[, 1]],  # Coordenada X del nodo inicial
  ystart = phylo_coords$Y[sum_tree_all$edge[, 1]],  # Coordenada Y del nodo inicial
  xend = phylo_coords$X[sum_tree_all$edge[, 2]],    # Coordenada X del nodo final
  yend = phylo_coords$Y[sum_tree_all$edge[, 2]]     # Coordenada Y del nodo final
)

# Función para agregar los convex hulls
add_convex_hull <- function(df) {
  df[chull(df$PC1, df$PC2), ]
}

# Crear los polígonos convexos para las familias
convex_hulls <- PCAs_scores %>%
  group_by(family) %>%
  do(add_convex_hull(.))

# Graficar las relaciones filogenéticas, los puntos de especies, y los convex hulls
pca_plot <- ggplot() +
  geom_segment(data = coords_edges, aes(x = xstart, y = ystart, xend = xend, yend = yend),
               color = "gray", linetype = "solid") +  # Conexiones filogenéticas
  geom_point(data = PCAs_scores, aes(x = PC1, y = PC2, color = family), size = 3) +  # Puntos de especies
  geom_polygon(data = convex_hulls, aes(x = PC1, y = PC2, fill = family, group = family), 
               alpha = 0.2, color = "black", linetype = "dashed") +  # Polígonos convexos
  theme_minimal() +
  labs(title = "PCA Filogenético con Morfometría Externa", x = "PC1", y = "PC2")

# Mostrar el gráfico final con las conexiones filogenéticas, puntos y convex hulls
print(pca_plot)
