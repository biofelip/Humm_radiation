setwd("C:/Users/felip/Documents/Proyectos_Git/Humm_radiation/")

############ 0. Ajuste de modelos evolutivos #####################################

## antes de cualquier cosa necesitamos saber cual es el modelo de cambio que mejor se ajusta a nuestras variables
## con esta información podemos sabe como corregir la inercia filogenetica de los datos
## voy a probar los modelos mas clasicos: BM, OU, EB, ARD, ER y ACDC.
## una ves sepa con el Aic cual es el mas "ajustado" procederemos a usarlo en el calculo de PGLS

# Cargar las librerías necesarias
# cambio para juan
library(ape)
library(geiger)
library(dplyr)
library(phangorn)

# Cargar el archivo de datos morfológicos y el árbol filogenético
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", header = TRUE, sep = ";", fileEncoding = "UTF-8")
subset_trees <- read.nexus("Pyhlo/Subset_100_trees.nex")
sum_tree_all <- maxCladeCred(subset_trees)

# Renombrar la columna de especies
morpho_data$Species3 <- gsub(" ", "_", morpho_data$Species3)  # Reemplazar espacios por guiones bajos
rownames(morpho_data) <- morpho_data$Species3  # Usar Species3 como rownames

# Definir las variables morfométricas de interés
morpho_measures <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", 
                     "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance", 
                     "Secondary1", "Tail.Length", "Mass")

# Filtrar solo las columnas necesarias
morpho_data_filtered <- morpho_data[, c("Species3", morpho_measures)]

# Verificar si hay valores no numéricos y convertir correctamente
# Los errores anteriores pueden deberse a caracteres no numéricos como comas o caracteres especiales.
morpho_data_filtered[morpho_measures] <- lapply(morpho_data_filtered[morpho_measures], function(x) {
  as.numeric(gsub(",", ".", gsub("[^0-9.-]", "", x)))  # Reemplazar caracteres no numéricos y convertir a numérico
})

# Mostrar cuántos NAs existen en cada columna para revisar problemas
cat("Número de valores NA en cada columna después de la conversión:\n")
print(colSums(is.na(morpho_data_filtered[morpho_measures])))

# Eliminar filas con NA en las variables morfométricas
morpho_data_filtered <- na.omit(morpho_data_filtered)

# Filtrar las especies presentes en el árbol filogenético
species_in_tree <- sum_tree_all$tip.label
morpho_data_filtered <- morpho_data_filtered[rownames(morpho_data_filtered) %in% species_in_tree, ]

# Verificar si hay diferencias en los nombres entre el árbol y los datos
missing_in_data <- setdiff(species_in_tree, rownames(morpho_data_filtered))
missing_in_tree <- setdiff(rownames(morpho_data_filtered), species_in_tree)

if (length(missing_in_data) > 0) {
  cat("Especies en el árbol que faltan en los datos:\n")
  print(missing_in_data)
}

if (length(missing_in_tree) > 0) {
  cat("Especies en los datos que faltan en el árbol:\n")
  print(missing_in_tree)
}

# Asegurarse de que el orden de las especies en la base de datos y en el árbol coincida
morpho_data_filtered <- morpho_data_filtered[match(species_in_tree, rownames(morpho_data_filtered)), ]

# Verificar que las especies en los datos coincidan con las del árbol
if (!all(rownames(morpho_data_filtered) == species_in_tree)) {
  stop("El orden de las especies no coincide entre los datos y el árbol.")
}

# Función para ejecutar un solo modelo
run_fitContinuous <- function(model, phy, data) {
  results_list <- list()  # Lista para almacenar los resultados de cada variable
  
  for (measure in morpho_measures) {
    cat("\nEjecutando el modelo", model, "para la variable", measure, "\n")
    
    # Extraer los datos para la variable y asegurarse de que tienen nombres
    trait_data <- data[[measure]]
    names(trait_data) <- rownames(data)
    
    # Verificar que los datos sean numéricos
    if (!is.numeric(trait_data)) {
      cat("Error: los datos para la variable", measure, "no son numéricos\n")
      next
    }
    
    result <- tryCatch({
      fitContinuous(phy = phy, dat = trait_data, model = model)
    }, error = function(e) {
      cat("Error en el ajuste del modelo", model, "para la variable", measure, ": ", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      cat("Modelo:", model, "completado para", measure, "\n")
      results_list[[measure]] <- result
    } else {
      results_list[[measure]] <- NULL
    }
  }
  
  return(results_list)
}

# -------------------------------------------------------------------
# --- Ejecutar cada modelo para todas las variables ---
# -------------------------------------------------------------------

# Ejecutar el modelo Browniano (BM) para todas las variables
bm_results <- run_fitContinuous("BM", sum_tree_all, morpho_data_filtered)

# Ejecutar el modelo Ornstein-Uhlenbeck (OU) para todas las variables
ou_results <- run_fitContinuous("OU", sum_tree_all, morpho_data_filtered)

# Ejecutar el modelo Early Burst (EB) para todas las variables
eb_results <- run_fitContinuous("EB", sum_tree_all, morpho_data_filtered)

# Ejecutar el modelo ARD (tasas dependientes del estado) para todas las variables
ard_results <- run_fitContinuous("rate_trend", sum_tree_all, morpho_data_filtered)

# Ejecutar el modelo ER (tasa constante) para todas las variables
er_results <- run_fitContinuous("mean_trend", sum_tree_all, morpho_data_filtered)

# Ejecutar el modelo ACDC (Acelerado o Retardado) para todas las variables
acdc_results <- run_fitContinuous("lambda", sum_tree_all, morpho_data_filtered)


# -------------------------------------------------------------------
# --- Comparar los modelos evaluando los valores de AIC ---
# -------------------------------------------------------------------

# Crear una tabla vacía para almacenar los AIC con las variables como filas y los modelos como columnas
aic_table <- data.frame(
  Variable = morpho_measures,
  BM = NA, OU = NA, EB = NA, ARD = NA, ER = NA, ACDC = NA
)

# Rellenar la tabla de AIC para cada modelo y cada variable
for (measure in morpho_measures) {
  aic_table[aic_table$Variable == measure, "BM"] <- if (!is.null(bm_results[[measure]])) bm_results[[measure]]$opt$aic else NA
  aic_table[aic_table$Variable == measure, "OU"] <- if (!is.null(ou_results[[measure]])) ou_results[[measure]]$opt$aic else NA
  aic_table[aic_table$Variable == measure, "EB"] <- if (!is.null(eb_results[[measure]])) eb_results[[measure]]$opt$aic else NA
  aic_table[aic_table$Variable == measure, "ARD"] <- if (!is.null(ard_results[[measure]])) ard_results[[measure]]$opt$aic else NA
  aic_table[aic_table$Variable == measure, "ER"] <- if (!is.null(er_results[[measure]])) er_results[[measure]]$opt$aic else NA
  aic_table[aic_table$Variable == measure, "ACDC"] <- if (!is.null(acdc_results[[measure]])) acdc_results[[measure]]$opt$aic else NA
}

# Mostrar la tabla de AIC
cat("\nTabla de AIC para todos los modelos y variables:\n")
print(aic_table)

# Guardar la tabla de AIC en un archivo CSV
write.csv(aic_table, "aic_model_results.csv", row.names = FALSE)


# -------------------------------------------------------------------
# --- Comparar los modelos evaluando los valores de AICc ---
# -------------------------------------------------------------------

# Crear una tabla vacía para almacenar los AICc con las variables como filas y los modelos como columnas
aic_table <- data.frame(
  Variable = morpho_measures,
  BM = NA, OU = NA, EB = NA, ARD = NA, ER = NA, ACDC = NA
)

# Rellenar la tabla de AICc para cada modelo y cada variable
for (measure in morpho_measures) {
  aic_table[aic_table$Variable == measure, "BM"] <- if (!is.null(bm_results[[measure]])) bm_results[[measure]]$opt$aicc else NA
  aic_table[aic_table$Variable == measure, "OU"] <- if (!is.null(ou_results[[measure]])) ou_results[[measure]]$opt$aicc else NA
  aic_table[aic_table$Variable == measure, "EB"] <- if (!is.null(eb_results[[measure]])) eb_results[[measure]]$opt$aicc else NA
  aic_table[aic_table$Variable == measure, "ARD"] <- if (!is.null(ard_results[[measure]])) ard_results[[measure]]$opt$aicc else NA
  aic_table[aic_table$Variable == measure, "ER"] <- if (!is.null(er_results[[measure]])) er_results[[measure]]$opt$aicc else NA
  aic_table[aic_table$Variable == measure, "lambda"] <- if (!is.null(acdc_results[[measure]])) acdc_results[[measure]]$opt$aicc else NA
}

# Mostrar la tabla de AICc
cat("\nTabla de AICc para todos los modelos y variables:\n")
print(aic_table)

# Guardar la tabla de AICc en un archivo CSV
write.csv(aic_table, "aicc_model_results.csv", row.names = FALSE)


######
##### Los resultados medicen que Lambda es el modelo mas adecuado para mis datos, por lo tanto es el que usare en el PGLS
#### tanto en e AIC como en el AICc los valores mas bajos son los de lamba y por una diferencia mayor a 4
#### Lambda es un modelo flexible que ajusta el grado de dependencia filogenética en los datos, proporcionando 
### un valor óptimo entre 0 (sin dependencia filogenética) y 1 (dependencia total, que sería masomenos el Browniano)
### nota: en el analiis aparece como ACDC porque me equivoque pero luego no queria cambiar todo el codigo para ajustar el nombre




############ AHORA SI LOS ANALISIS ###########

########## 1. Supuestos de Normalidad, Homocedasticidad y PGLS para controlar por el tamaño #####################

# Cargar las librerías necesarias
library(ape)
library(nlme)
library(phytools)
library(dplyr)
library(lmtest)

# Cargar los datos morfológicos
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", 
                        header = TRUE, sep = ";", fileEncoding = "UTF-8")

# Cargar el árbol filogenético (subset de 100 árboles)
subset_trees <- read.nexus("Pyhlo/Subset_100_trees.nex")

# Crear un árbol de máxima credibilidad a partir del subset de los 100 árboles
sum_tree_all <- maxCladeCred(subset_trees)

# Definir las medidas morfométricas que queremos transformar y analizar
morpho_measures <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", 
                     "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance", 
                     "Secondary1", "Tail.Length", "Mass")  # Agregamos 'Mass' como covariable

# Verificar qué columnas de morpho_measures existen en morpho_data
existing_columns <- morpho_measures %in% colnames(morpho_data)

# Filtrar solo las columnas existentes
valid_measures <- morpho_measures[existing_columns]

# Convertir las columnas que existen a formato numérico
morpho_data[valid_measures] <- lapply(morpho_data[valid_measures], function(x) as.numeric(gsub(",", ".", x)))

# Reemplazar espacios o guiones en los nombres de las especies para que coincidan entre los datos y el árbol
morpho_data$Species3 <- gsub(" ", "_", morpho_data$Species3)
morpho_data$Species3 <- gsub("-", "_", morpho_data$Species3)

# Log-transformación de las medidas morfométricas para mejorar la normalidad
log_transformed_data <- morpho_data %>%
  mutate(across(all_of(valid_measures), log))

# Filtrar solo las especies que están presentes en el árbol
log_transformed_data <- log_transformed_data[log_transformed_data$Species3 %in% sum_tree_all$tip.label, ]

# Reordenar los datos para que coincidan con el árbol
log_transformed_data <- log_transformed_data[match(sum_tree_all$tip.label, log_transformed_data$Species3), ]

# Recortar el árbol para que coincida con las especies presentes en los datos
tree_subset <- drop.tip(sum_tree_all, setdiff(sum_tree_all$tip.label, log_transformed_data$Species3))

# Verificar si el árbol fue correctamente recortado
if (is.null(tree_subset)) {
  stop("El árbol ha sido recortado a un tamaño nulo. Revisa que las especies coincidan correctamente entre el árbol y los datos.")
}

# ---------------------------------
# Evaluar normalidad y homocedasticidad
# ---------------------------------

# Crear un dataframe para almacenar los resultados de normalidad y homocedasticidad
supuestos_de_analisis <- data.frame(
  Variable = character(),
  p_valor_normalidad = numeric(),
  p_valor_homocedasticidad = numeric(),
  stringsAsFactors = FALSE
)

# Evaluar la normalidad usando la prueba Kolmogorov-Smirnov para cada medida
for (trait in valid_measures) {
  # Normalidad (Kolmogorov-Smirnov)
  ks_test <- ks.test(log_transformed_data[[trait]], "pnorm", 
                     mean(log_transformed_data[[trait]], na.rm = TRUE), 
                     sd(log_transformed_data[[trait]], na.rm = TRUE))
  
  # Prueba de homocedasticidad (Breusch-Pagan)
  if (trait != "Mass") {  # No tiene sentido probar homocedasticidad con Mass contra sí mismo
    lm_model <- lm(log_transformed_data[[trait]] ~ log_transformed_data$Mass)
    bp_test <- bptest(lm_model)
    
    # Almacenar resultados en la tabla
    supuestos_de_analisis <- rbind(supuestos_de_analisis, 
                                   data.frame(Variable = trait, 
                                              p_valor_normalidad = ks_test$p.value,
                                              p_valor_homocedasticidad = bp_test$p.value))
  }
}

# Mostrar la tabla con los resultados de normalidad y homocedasticidad
cat("Resultados de las pruebas de supuestos del análisis:\n")
print(supuestos_de_analisis)

# ---------------------------------
# PGLS después de validar normalidad y homocedasticidad
# ---------------------------------

# Definir una lista para almacenar los residuos obtenidos
residuals_list <- list()

# Ajustar un modelo PGLS para cada medida morfométrica, utilizando "Mass" como covariable
for (trait in valid_measures) {
  if (trait != "Mass") {  # No hacemos regresión de Mass contra sí mismo
    # Ajustar el modelo PGLS correctamente
    formula_pgls <- as.formula(paste(trait, "~ Mass"))
    
    # Asegurarnos de que las especies están alineadas con el árbol
    model_pgls <- gls(formula_pgls, 
                      correlation = corPagel(value = 0.5, phy = tree_subset, form = ~Species3),  # corPagel es el equivalente a "lambda" en PGLS
                      data = log_transformed_data, method = "ML")
    
    # Obtener los residuos y almacenarlos
    residuals_list[[trait]] <- residuals(model_pgls)
    
    # Mostrar un resumen del modelo
    cat("\nResumen del PGLS para", trait, ":\n")
    print(summary(model_pgls))
  }
}

# Guardar los residuos en un dataframe para análisis posteriores
residuals_df <- as.data.frame(residuals_list)

# Mostrar un resumen de los residuos
cat("\nResumen de los residuos obtenidos:\n")
print(head(residuals_df))

# Puedes guardar los residuos en un archivo CSV si lo deseas
write.csv(residuals_df, "residuos_pgls.csv")



#### como los resultados de p de "mass" en todos los modelos de PGLS estan dando sospechosamente cercanos a 0 absoluto
### voy a graficar cada variable contra mass para reviar si efectivamente la relación es extrema

# Graficar la relación entre Mass y cada variable morfométrica
for (trait in valid_measures) {
  if (trait != "Mass") {  # No hacemos regresión de Mass contra sí mismo
    plot(log_transformed_data$Mass, log_transformed_data[[trait]], 
         main = paste("Relación entre Mass y", trait),
         xlab = "Log(Mass)", ylab = paste("Log(", trait, ")", sep = ""),
         pch = 16)
    abline(lm(log_transformed_data[[trait]] ~ log_transformed_data$Mass), col = "red")
  }
}



##### tambien voy a calcular los intervalos de confianza para usarlos como herramienta para corroborar ese resultado ####

# Calcular intervalos de confianza para los coeficientes de PGLS
for (trait in valid_measures) {
  if (trait != "Mass") {  # No hacemos regresión de Mass contra sí mismo
    formula_pgls <- as.formula(paste(trait, "~ Mass"))
    
    # Asegurarnos de que las especies están alineadas con el árbol
    model_pgls <- gls(formula_pgls, 
                      correlation = corBrownian(phy = tree_subset, form = ~Species3), 
                      data = log_transformed_data, method = "ML")
    
    # Calcular intervalos de confianza
    confint_pgls <- confint(model_pgls)
    
    # Mostrar los intervalos de confianza
    cat("\nIntervalos de confianza para los coeficientes del PGLS para", trait, ":\n")
    print(confint_pgls)
  }
}


#### los intervalos de confianza son muy estrechos y no incluyen el cero, lo cual soporta la idea de que el 
### valor de p de "mass" en los modelos es muy cercano a cero


### Adicionalmente, voy a correr un modelo lineal simple para ver si se comporta igual
### es decir, un modelo sin tener en cuenta la filogenia 

# Ajustar un modelo OLS para cada variable morfométrica
for (trait in valid_measures) {
  if (trait != "Mass") {
    formula_ols <- as.formula(paste(trait, "~ Mass"))
    model_ols <- lm(formula_ols, data = log_transformed_data)
    
    # Mostrar el resumen del modelo OLS
    cat("\nResumen del modelo OLS para", trait, ":\n")
    print(summary(model_ols))
  }
}

###### efectivamente los resultados de este mdoelos lienal simple corroboran que el peso de "mass"
##### sobre todas als variables es muuy grande, muy cercano a cer en el valor de p, esto sumado a los
#### resultados del calculo de inetrvalos de confianza soprtan muy bien los resultados de PGLS, con lo cual
#### podemos estar tranquilos y asumir que los resultados del PGLS estan bien :)






#### ahora una figura bonita de los intervalos de confianza para soportar los resultados de PGLS en el control por tamaño


# Cargar las librerías necesarias
library(ggplot2)
library(nlme)

# Definir un dataframe para almacenar los coeficientes y sus intervalos de confianza
intervalos_confianza <- data.frame(
  Variable = character(),
  Coeficiente = numeric(),
  IC_inferior = numeric(),
  IC_superior = numeric(),
  stringsAsFactors = FALSE
)

# Calcular los intervalos de confianza para los coeficientes del modelo PGLS
for (trait in valid_measures) {
  if (trait != "Mass") {  # No hacemos regresión de Mass contra sí mismo
    formula_pgls <- as.formula(paste(trait, "~ Mass"))
    
    # Ajustar el modelo PGLS
    model_pgls <- gls(formula_pgls, 
                      correlation = corBrownian(phy = tree_subset, form = ~Species3), 
                      data = log_transformed_data, method = "ML")
    
    # Extraer los coeficientes e intervalos de confianza
    coef_pgls <- coef(model_pgls)
    confint_pgls <- confint(model_pgls)
    
    # Almacenar los resultados en el dataframe
    intervalos_confianza <- rbind(intervalos_confianza, 
                                  data.frame(Variable = trait, 
                                             Coeficiente = coef_pgls["Mass"], 
                                             IC_inferior = confint_pgls["Mass", 1], 
                                             IC_superior = confint_pgls["Mass", 2]))
  }
}

# Crear la figura con ggplot2
ggplot(intervalos_confianza, aes(x = Variable, y = Coeficiente)) +
  geom_point(size = 4, color = "blue") +  # Puntos para los coeficientes
  geom_errorbar(aes(ymin = IC_inferior, ymax = IC_superior), width = 0.2, color = "black") +  # Líneas para los intervalos de confianza
  theme_minimal() +  # Estilo limpio
  labs(title = "Intervalos de confianza de los coeficientes del PGLS",
       x = "Variable morfométrica", y = "Coeficiente de Mass") +
  coord_flip() +  # Voltear el eje para mejorar la legibilidad
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título









########################### 2.   MANOVA Normal ################################

## este analisis sirve para ver si las especies de distintas familias tienen fenotipos morfologicos multivariados diferentes
## vamos a usar como input los residuos de PGLS calculados en el analisis anterior. Como usamos los residuos como un input
## entonces no es necesario que este analisis sea un "MANOVA FILOGENETICO". Ya la inercia filogenetica fue corregida en el PGLS


# Cargar las librerías necesarias
library(dplyr)
library(ggplot2)
library(MASS)  # Para DFA (Linear Discriminant Analysis)
library(ggpubr)  # Para combinar gráficos
library(ape)  # Para trabajar con árboles filogenéticos

# Cargar los residuos obtenidos del PGLS y la base de datos taxonómica
residuos_pgls <- read.csv("residuos_pgls.csv", row.names = 1)
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", header = TRUE, sep = ";", fileEncoding = "UTF-8")

# Asegurarnos de que los nombres de las especies en el archivo morfológico coincidan con el formato del árbol y residuos
morpho_data$Species3 <- gsub(" ", "_", morpho_data$Species3)  # Reemplazar espacios por guiones bajos si es necesario
rownames(morpho_data) <- morpho_data$Species3  # Usar Species3 como rownames

# Seleccionar las columnas de interés
taxonomic_info <- morpho_data[, c("Species3", "Family3", "Order3", "Clade")]

# Combinar los residuos con la información taxonómica
residuos_pgls_combined <- cbind(taxonomic_info, residuos_pgls)

# -------------------------------------------------------------------
# --- MANOVA estándar con todas las especies del conjunto de datos ---
# -------------------------------------------------------------------

# Convertir "Family3" en un factor
residuos_pgls_combined$Family3 <- as.factor(residuos_pgls_combined$Family3)

# Verificar los niveles de Family3
cat("Niveles de Family3 en el conjunto de datos:\n")
print(levels(residuos_pgls_combined$Family3))

# Extraer las variables morfométricas (excluir la taxonomía)
Y <- as.matrix(residuos_pgls_combined[, -c(1:5)])  # Excluir las primeras 5 columnas de taxonomía

# Verificar si hay valores NA en los datos
cat("Número de valores faltantes en Y:\n")
print(colSums(is.na(Y)))

# Eliminar cualquier fila con valores NA
residuos_pgls_combined <- residuos_pgls_combined[complete.cases(Y), ]

# Volver a extraer las variables morfométricas sin NA
Y <- as.matrix(residuos_pgls_combined[, -c(1:5)])

# Verificar que haya al menos dos niveles en Family3
if (length(unique(residuos_pgls_combined$Family3)) < 2) {
  stop("El análisis MANOVA requiere al menos dos niveles en 'Family3'.")
}

# Ejecutar el análisis MANOVA con los residuos
manova_result <- manova(Y ~ Family3, data = residuos_pgls_combined)

# Mostrar los resultados del MANOVA
cat("Resultados del análisis MANOVA:\n")
print(summary(manova_result))


### el resultado de manova es muy significativo, es decir que las familias varian en el morfoespacio multivariado
### para saber especificamente que variables son las que se diferencian hacemos un test "post-hoc", en este caso Tukey

# -------------------------------------------------------------------
# --- Pruebas Post-hoc para el MANOVA usando TukeyHSD ---
# -------------------------------------------------------------------

# Realizar ANOVA univariado y usar TukeyHSD para cada variable morfométrica como post-hoc
posthoc_results <- lapply(1:ncol(Y), function(i) {
  model <- aov(Y[, i] ~ residuos_pgls_combined$Family3)
  TukeyHSD(model)
})

# Mostrar los resultados de las pruebas post-hoc
cat("\nResultados de las pruebas post-hoc (TukeyHSD):\n")
print(posthoc_results)

# -------------------------------------------------------------------
# --- Cálculo de intervalos de confianza ---
# -------------------------------------------------------------------

# Calcular intervalos de confianza usando TukeyHSD para cada variable
conf_intervals <- lapply(1:ncol(Y), function(i) {
  model <- aov(Y[, i] ~ residuos_pgls_combined$Family3)
  as.data.frame(TukeyHSD(model)$`residuos_pgls_combined$Family3`)
})

# Unir los intervalos de confianza en un solo data frame
ci_data <- do.call(rbind, lapply(seq_along(conf_intervals), function(i) {
  data.frame(Variable = colnames(Y)[i],
             Estimate = conf_intervals[[i]]$diff,
             LowerCI = conf_intervals[[i]]$lwr,
             UpperCI = conf_intervals[[i]]$upr,
             Comparison = rownames(conf_intervals[[i]]))
}))

# -------------------------------------------------------------------
# --- Crear figura para los intervalos de confianza ---
# -------------------------------------------------------------------

# Cargar librería adicional para paletas de colores
library(RColorBrewer)

# Crear un conjunto de colores únicos para las variables
color_palette <- brewer.pal(n = length(unique(ci_data$Variable)), name = "Set1")

# Crear un gráfico de intervalos de confianza para cada variable, con colores diferentes
ci_plots <- lapply(seq_along(unique(ci_data$Variable)), function(i) {
  variable <- unique(ci_data$Variable)[i]
  ggplot(ci_data[ci_data$Variable == variable, ], aes(x = Comparison, y = Estimate, color = Variable)) +
    geom_pointrange(aes(ymin = LowerCI, ymax = UpperCI), position = position_dodge(width = 0.5)) +
    coord_flip() +
    labs(title = paste("Intervalos de confianza de", variable),
         x = "Comparación", y = "Estimación") +
    scale_color_manual(values = color_palette[i]) +  # Asignar color específico a cada variable
    theme_minimal()
})

# Combinar los 8 gráficos en un solo panel, 2 columnas por 4 filas
ci_panel <- ggarrange(plotlist = ci_plots, ncol = 2, nrow = 4)

# Guardar el panel de intervalos de confianza en un archivo PDF en formato A3 vertical
# A3 en pulgadas es aproximadamente 11.7 x 16.5 (29.7 cm x 42 cm)
pdf("intervalos_confianza_panel_A3.pdf", width = 11.7, height = 16.5)  # Ajuste a A3 vertical
print(ci_panel)
dev.off()

# Guardar el panel en formato SVG
svg("intervalos_confianza_panel_A3.svg", width = 11.7, height = 16.5)
print(ci_panel)
dev.off()  # Cerrar el archivo SVG

# Guardar el panel en formato EPS
postscript("intervalos_confianza_panel_A3.eps", width = 11.7, height = 16.5, paper = "special", horizontal = FALSE)
print(ci_panel)
dev.off()  # Cerrar el archivo EPS




#### en el plot se ve que ninguna variable por si sola discrimina todos los grupos (todas pasan por el cero en alguna combinacion)
#### sino que mas bien algunas variables discriminan unos grupos y otras discriminan otros.


### ahora hacemos un analisis discriminante para ver como se organizan las familias en el morfoespacio, dada la diferencia morfometrica

# -------------------------------------------------------------------
# --- Discriminant Function Analysis (DFA) ---
# -------------------------------------------------------------------

# Ejecutar el análisis DFA utilizando Linear Discriminant Analysis (LDA)
dfa_result <- lda(Family3 ~ ., data = residuos_pgls_combined[, -c(1, 3:5)])  # Usar Family3 para discriminar

# Mostrar los resultados del DFA
cat("Resultados del análisis DFA:\n")
print(dfa_result)

# Predecir las clases con el DFA
predictions <- predict(dfa_result)

# -------------------------------------------------------------------
# --- Crear gráficos para visualizar los resultados ---
# -------------------------------------------------------------------

# Gráfico de los scores DFA
dfa_scores <- as.data.frame(predictions$x)
dfa_scores$Family3 <- residuos_pgls_combined$Family3

# Plot DFA Scores para DL1 vs DL2
p1 <- ggplot(dfa_scores, aes(x = LD1, y = LD2, color = Family3)) +
  geom_point(size = 2) +  # Tamaño reducido a 2
  labs(title = "Análisis de Función Discriminante (DFA) - DL1 vs DL2",
       x = "Discriminante Lineal 1", y = "Discriminante Lineal 2") +
  theme_minimal()

# Plot DFA Scores para DL1 vs DL3
p2 <- ggplot(dfa_scores, aes(x = LD1, y = LD3, color = Family3)) +
  geom_point(size = 2) +  # Tamaño reducido a 2
  labs(title = "Análisis de Función Discriminante (DFA) - DL1 vs DL3",
       x = "Discriminante Lineal 1", y = "Discriminante Lineal 3") +
  theme_minimal()

# Plot DFA Scores para DL1 vs DL4
p3 <- ggplot(dfa_scores, aes(x = LD1, y = LD4, color = Family3)) +
  geom_point(size = 2) +  # Tamaño reducido a 2
  labs(title = "Análisis de Función Discriminante (DFA) - DL1 vs DL4",
       x = "Discriminante Lineal 1", y = "Discriminante Lineal 4") +
  theme_minimal()

# Plot DFA Scores para DL2 vs DL3
p4 <- ggplot(dfa_scores, aes(x = LD2, y = LD3, color = Family3)) +
  geom_point(size = 2) +  # Tamaño reducido a 2
  labs(title = "Análisis de Función Discriminante (DFA) - DL2 vs DL3",
       x = "Discriminante Lineal 2", y = "Discriminante Lineal 3") +
  theme_minimal()

# Combinar todas las figuras de DFA en un solo panel
dfa_panel <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Mostrar el panel de gráficos de DFA
print(dfa_panel)



##### ahora un plot bonito de eso en 3D #####

# Cargar la librería plotly
library(plotly)

# Crear un dataframe con los scores de DFA (LD1, LD2, LD3)
dfa_scores <- as.data.frame(predictions$x)
dfa_scores$Family3 <- residuos_pgls_combined$Family3

# Crear el gráfico 3D para LD1, LD2 y LD3
plot_3d <- plot_ly(dfa_scores, x = ~LD1, y = ~LD2, z = ~LD3, color = ~Family3, colors = "Set1", type = "scatter3d", mode = "markers") %>%
  layout(
    title = "Análisis de Función Discriminante (DFA) - LD1 vs LD2 vs LD3",
    scene = list(
      xaxis = list(title = "Discriminante Lineal 1 (LD1)"),
      yaxis = list(title = "Discriminante Lineal 2 (LD2)"),
      zaxis = list(title = "Discriminante Lineal 3 (LD3)")
    )
  )

# Mostrar el gráfico 3D
plot_3d




##### plot de resumen de DF scores

library(ggplot2)
library(MASS)  # Para la función lda

# -------------------------------------------------------------------
# --- Discriminant Function Analysis (DFA) ---
# -------------------------------------------------------------------

# Ejecutar el análisis DFA utilizando Linear Discriminant Analysis (LDA)
dfa_result <- lda(Family3 ~ ., data = residuos_pgls_combined[, -c(1, 3:5)])  # Usar Family3 para discriminar

# Mostrar los resultados del DFA
cat("Resultados del análisis DFA:\n")
print(dfa_result)

# Predecir las clases con el DFA
predictions <- predict(dfa_result)

# -------------------------------------------------------------------
# --- Crear gráfico estilo el mostrado ---
# -------------------------------------------------------------------

# Crear un dataframe con los scores DFA
dfa_scores <- as.data.frame(predictions$x)
dfa_scores$Family3 <- residuos_pgls_combined$Family3

# Crear el gráfico para DF scores
p <- ggplot(dfa_scores, aes(x = Family3, y = LD1, color = Family3)) +
  geom_jitter(size = 3, width = 0.2) +  # Añadir los puntos con un poco de dispersión horizontal para evitar solapamientos
  theme_minimal() +
  labs(title = "Discriminant Function Analysis (DFA) - LD1",
       x = "Familias", y = "DF scores") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar las etiquetas del eje x
  theme(legend.position = "none")  # Ocultar la leyenda

# Mostrar el gráfico
print(p)



### si bien, se diferencian mas unos clados de otros, se ve que colibries, una vez corregido e tamaño, sigue teniendo un espacio
### unico significativo, aunque una parte del espacio parece estar compartido por Apodidae.



## ahora vamos a graficar un scatterplot que muestre la influencia de la masa sobre cada variable, usando como input los datos transformados a log.
## esto es lo que llamamos las "trayectorias alométricas", que deperminan el peso del tamaño sobre la modificacion de cada variable morfometrica

# Cargar las librerías necesarias
library(ggplot2)
library(gridExtra)

# Seleccionar todas las variables morfométricas que desees graficar
morpho_variables <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", 
                      "Beak.Depth", "Tarsus.Length", "Wing.Length", 
                      "Kipps.Distance", "Tail.Length")

# Crear una lista para almacenar los gráficos
plots_list <- list()

# Crear los gráficos de dispersión (scatterplots) para cada variable
for (variable_morfometrica in morpho_variables) {
  
  # Crear un dataframe con los datos transformados y la covariable Mass
  data_for_plot <- data.frame(
    Variable = log_transformed_data[[variable_morfometrica]],
    Mass = log_transformed_data$Mass,
    Family3 = log_transformed_data$Family3
  )
  
  # Verificar si las columnas existen
  if (any(is.na(data_for_plot$Variable)) || any(is.na(data_for_plot$Mass))) {
    stop(paste("Algunos valores de la variable", variable_morfometrica, "o de Mass no están disponibles."))
  }
  
  # Crear el gráfico de dispersión (scatterplot)
  scatter_plot <- ggplot(data_for_plot, aes(x = Mass, y = Variable, color = Family3)) +
    geom_point(size = 1.3) +  # Tamaño de los puntos reducido a la mitad (1)
    geom_smooth(method = "lm", se = FALSE) +  # Línea de tendencia lineal para cada familia
    labs(title = paste(variable_morfometrica, "vs Mass"),
         x = "Mass (log)", y = paste("Log de", variable_morfometrica),
         color = "Familia") +  # Etiquetas de los ejes y leyenda
    theme_minimal() +  # Tema minimalista
    theme(legend.position = "right")  # Colocar la leyenda a la derecha
  
  # Agregar el gráfico a la lista de gráficos
  plots_list[[variable_morfometrica]] <- scatter_plot
}

# Crear los paneles 2x2 con los gráficos
panel_1 <- grid.arrange(plots_list[[1]], plots_list[[2]], plots_list[[3]], plots_list[[4]], ncol = 2, nrow = 2)
panel_2 <- grid.arrange(plots_list[[5]], plots_list[[6]], plots_list[[7]], plots_list[[8]], ncol = 2, nrow = 2)

# Mostrar los paneles
print(panel_1)
print(panel_2)










































































#################  phylomorphospace ####################

# Cargar las librerías necesarias
library(ape)
library(phytools)
library(RColorBrewer)

# Cargar los residuos obtenidos del PGLS y la base de datos taxonómica
residuos_pgls <- read.csv("residuos_pgls.csv", row.names = 1)
morpho_data <- read.csv("BD_morpho/Humm_Strisores_all_AVONET.csv", header = TRUE, sep = ";", fileEncoding = "UTF-8")

# Asegurarnos de que los nombres de las especies en el archivo morfológico coincidan con el formato del árbol y residuos
morpho_data$Species3 <- gsub(" ", "_", morpho_data$Species3)  # Reemplazar espacios por guiones bajos si es necesario
rownames(morpho_data) <- morpho_data$Species3  # Usar Species3 como rownames

# Seleccionar las columnas de interés
taxonomic_info <- morpho_data[, c("Species3", "Family3", "Order3", "Clade")]

# Combinar los residuos con la información taxonómica
residuos_pgls_combined <- cbind(taxonomic_info, residuos_pgls)

# Asegurarse de que las especies del árbol y de los residuos coincidan
common_species <- intersect(rownames(residuos_pgls_combined), sum_tree_all$tip.label)
if (length(common_species) == 0) {
  stop("No hay coincidencias entre las especies en el árbol y en los residuos.")
}

# Filtrar el árbol y los datos para que contengan solo las especies comunes
tree_subset <- drop.tip(sum_tree_all, setdiff(sum_tree_all$tip.label, common_species))
residuos_pgls_filtered <- residuos_pgls_combined[rownames(residuos_pgls_combined) %in% common_species, ]

# Asegurar que el orden de las especies en los residuos coincida con el del árbol
residuos_pgls_filtered <- residuos_pgls_filtered[match(tree_subset$tip.label, rownames(residuos_pgls_filtered)), ]

# Calcular el PCA sobre los residuos del PGLS
residuals_matrix <- as.matrix(residuos_pgls_filtered[, !colnames(residuos_pgls_filtered) %in% c("Family3", "Species3", "Order3", "Clade")])
pca_res <- prcomp(residuals_matrix, center = TRUE, scale. = TRUE)

# Obtener los scores de las PCs
pca_scores <- as.data.frame(pca_res$x[, 1:3])  # PC1, PC2 y PC3
pca_scores$Family3 <- residuos_pgls_filtered$Family3

# Asignar colores a las familias
family_list <- unique(pca_scores$Family3)
family_colors <- setNames(rainbow(length(family_list)), family_list)

# Cambiar el color de Trochilidae a fucsia
family_colors["Trochilidae"] <- "magenta"

tip_colors <- family_colors[pca_scores$Family3]

# Asignar colores a las ramas del árbol basadas en las familias
tree_colored <- tree_subset  # Copia del árbol para colorear los linajes
edge_colors <- rep("black", nrow(tree_colored$edge))  # Color inicial (negro) para las ramas

for (i in seq_along(family_list)) {
  family_species <- rownames(pca_scores[pca_scores$Family3 == family_list[i], ])
  
  # Asegurarse de que haya al menos dos especies en la familia para calcular el MRCA
  if (length(family_species) > 1) {
    mrca_node <- findMRCA(tree_subset, tips = family_species)
    
    # Verificar que el MRCA sea un valor numérico válido y pintar el subárbol
    if (!is.na(mrca_node) && is.numeric(mrca_node)) {
      branches <- paintSubTree(tree_colored, mrca_node, state = family_list[i])
      
      # Colorear las ramas asociadas a la familia
      edge_indices <- which(tree_colored$edge[, 1] == mrca_node | tree_colored$edge[, 2] %in% which(tree_colored$tip.label %in% family_species))
      edge_colors[edge_indices] <- family_colors[family_list[i]]
    }
  }
}

# Para asegurarnos de que los colores de las ramas sean manejados correctamente
edge_colors <- setNames(edge_colors, tree_colored$edge[, 2])

# Phylomorphospace plot de PC1 vs PC2
par(mfrow = c(1, 2))  # Dos gráficos en un mismo panel
phylomorphospace(tree_colored, pca_res$x[, 1:2], label = "off", node.size = 0,
                 control = list(col.node = "black", col.edge = edge_colors, node.by.map = TRUE))

# Añadir los puntos coloreados por familia
points(pca_res$x[, 1], pca_res$x[, 2], col = tip_colors, pch = 19, cex = 1.5)

# Añadir leyenda
legend("topright", legend = family_list, col = family_colors, pch = 19, title = "Familia")

# Phylomorphospace plot de PC1 vs PC3
phylomorphospace(tree_colored, pca_res$x[, c(1, 3)], label = "off", node.size = 0,
                 control = list(col.node = "black", col.edge = edge_colors, node.by.map = TRUE))

# Añadir los puntos coloreados por familia
points(pca_res$x[, 1], pca_res$x[, 3], col = tip_colors, pch = 19, cex = 1.5)

# Añadir leyenda
legend("topright", legend = family_list, col = family_colors, pch = 19, title = "Familia")



### ahora en 3d ##

library(phytools)

# Asegurar que el orden de las especies en los colores y el árbol coincidan
tip_colors_named <- setNames(tip_colors, tree_subset$tip.label)

# Crear el phylomorphospace 3D sin nombres y con los tips coloreados
phylomorphospace3d(tree_subset, pca_res$x[, 1:3], 
                   colors = tip_colors_named, 
                   label = "off", 
                   node.size = 0)  # Sin etiquetas

# Añadir leyenda para las familias
legend("topright", legend = family_list, pch = 19, col = family_colors, title = "Familia")



