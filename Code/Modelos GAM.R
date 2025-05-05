
library(caret)
library(MASS)
library(readxl)
library(dplyr)

############################################################### 1) Preparación de datos

datos2 <- read_excel("C:/Users/PC//Desktop/TESIS/reflectancia_todos_poligonos.xlsx")
datos3 <- datos2[, !(names(datos2) %in% c('Name', 'Empresa o DueÃ±o', 'Faena','Fecha'))]
datos3 <- datos3 %>%  mutate(across(all_of(setdiff(names(datos3), "Estado")), as.numeric))

############################################################### 1.1) Crear ínidces espectrales

B1 <- datos3$blue  # Azul
B2 <- datos3$green # Verde
B3 <- datos3$red # Rojo
B4 <- datos3$nir08 # NIR

# Cálculo de índices espectrales
datos3$ndvi <- (B4 - B3) / (B4 + B3)
datos3$dvi <- B4 - B3
datos3$rvi <- B4 / B3
datos3$sbi <- sqrt(B4^2 + B3^2)
datos3$evi <- 2.5 * (B4 - B3) / (B4 + 6 * B3 - 7.2 * B1 + 1)
datos3$ndwi <- (B2 - B4) / (B2 + B4)

median(datos3$As)

# Tratar la negatividad de las variables predictoras

variables_predictoras <- c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", 
                           "green", "red", "nir08", "coastal", "swir16", "swir22")

# Sumar 1 a todas las columnas de las variables predictoras
datos3[variables_predictoras] <- datos3[variables_predictoras] + 1


############################# Corr de Spearman vs Pearson

# Seleccionar variables predictoras y metales
variables_predictoras <- c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22", "ndvi", "dvi", "rvi", "sbi", "evi", "ndwi")
metales <- c("Ni(g/t)", "Cu(g/t)", "Zn(g/t)", "Co(g/t)", "As(g/t)")

# Crear un dataframe con las variables relevantes
datos_relevantes <- datos3[, c(variables_predictoras, metales)]

# Calcular correlaciones de Spearman
correlaciones_spearman <- cor(datos_relevantes, method = "spearman", use = "complete.obs")

# Calcular correlaciones de Pearson
correlaciones_pearson <- cor(datos_relevantes, method = "pearson", use = "complete.obs")

# Extraer correlaciones para cada metal con las variables predictoras
resultados_comparativos <- lapply(metales, function(metal) {
  spearman <- correlaciones_spearman[metal, variables_predictoras]
  pearson <- correlaciones_pearson[metal, variables_predictoras]
  data.frame(Variable = variables_predictoras, Spearman = spearman, Pearson = pearson)
})

# Mostrar resultados para cada metal
names(resultados_comparativos) <- metales

# Transformar datos para graficar
library(ggplot2)
library(tidyr)

for (metal in metales) {
  # Extraer los resultados del metal específico
  df_metal <- resultados_comparativos[[metal]]
  
  # Transformar a formato largo para graficar
  df_largo <- pivot_longer(df_metal, cols = c("Spearman", "Pearson"), 
                           names_to = "Correlación", values_to = "Valor")
  
  # Crear el gráfico
  p <- ggplot(df_largo, aes(x = Variable, y = Valor, fill = Correlación)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(Valor, 2)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3.5) +
    labs(
      title = paste("Correlación de Spearman vs Pearson -", metal),
      x = "Variables Predictoras",
      y = "Valor de Correlación",
      fill = "Tipo de Correlación"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  
  print(p)
}



####################### GAM

library(mgcv)


####### Co

# Renombrar la columna
datos3 <- datos3 %>%
  rename(Co = `Co(g/t)`)

set.seed(123)

# Configuración de variables, fórmulas y familias
configuracion_ni_co <- list(
  "Co" = list(
    formula = as.formula("Co ~ s(nir08) + s(swir16) + s(blue) + s(green) + s(coastal) + s(swir22)"),
    family = Gamma(link = "log")
  )
)

# Lista para almacenar resultados
resultados_ni_co <- list()

# Número de folds para k-fold cross-validation
k <- 10
mediana_co <- 15  # Mediana de la variable de respuesta "Co"

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- if (sum(!is.na(obs) & !is.na(pred)) > 0) {
    cor(pred, obs, use = "complete.obs")^2
  } else {
    NA
  }
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  rmse_estandarizado <- rmse / mediana_co  # RMSE estandarizado por la mediana
  return(list(R2 = r2, RMSE = rmse, Bias = bias, RMSE_Estandarizado = rmse_estandarizado))
}

# Iterar por cada metal
for (metal in names(configuracion_ni_co)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  config <- configuracion_ni_co[[metal]]
  formula <- config$formula
  family <- config$family
  
  # Filtrar datos sin NA
  variables <- all.vars(formula)
  datos_filtrados <- na.omit(datos3[, variables])
  
  # Verificar las columnas disponibles en el dataset
  print("Columnas disponibles en el dataset:")
  print(colnames(datos3))
  
  # Aplicar Brightness Normalization
  brightness_normalization <- function(df, vars) {
    brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
    brightness[brightness == 0] <- 1  # Evitar divisiones por cero
    df[, vars] <- df[, vars] / brightness
    return(df)
  }
  variables_predictoras <- setdiff(variables, "Co")
  datos_filtrados <- brightness_normalization(datos_filtrados, variables_predictoras)
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados$Co, k = k, list = TRUE)
  
  # Lista para almacenar métricas
  metricas_folds <- list(R2 = c(), RMSE = c(), Bias = c(), RMSE_Estandarizado = c())
  
  for (i in 1:k) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    # Ajustar el modelo GAM
    modelo <- tryCatch({
      gam(formula, data = TRAIN, family = family)
    }, error = function(e) {
      message(sprintf("Error en el fold %d: %s", i, e$message))
      return(NULL)
    })
    
    if (is.null(modelo)) next
    
    # Predecir en el conjunto de validación
    pred <- tryCatch({
      predict(modelo, newdata = VALIDATION, type = "response")
    }, error = function(e) {
      message("Error en la predicción: ", e$message)
      return(rep(NA, nrow(VALIDATION)))
    })
    
    # Calcular métricas
    metricas <- calcular_metricas(VALIDATION$Co, pred)
    metricas_folds$R2 <- c(metricas_folds$R2, metricas$R2)
    metricas_folds$RMSE <- c(metricas_folds$RMSE, metricas$RMSE)
    metricas_folds$Bias <- c(metricas_folds$Bias, metricas$Bias)
    metricas_folds$RMSE_Estandarizado <- c(metricas_folds$RMSE_Estandarizado, metricas$RMSE_Estandarizado)
  }
  
  # Promedio de métricas
  metricas_promedio <- lapply(metricas_folds, mean, na.rm = TRUE)
  
  # Guardar resultados
  resultados_ni_co[[metal]] <- list(modelo = modelo, metricas_promedio = metricas_promedio)
  
  # Imprimir métricas promedio
  cat("Métricas promedio para", metal, ":\n")
  print(metricas_promedio)
}

# Imprimir resultados globales
for (metal in names(resultados_ni_co)) {
  cat("\nResultados del modelo para:", metal, "\n")
  print(resultados_ni_co[[metal]])
}




########################## Zn y As

datos3 <- datos3 %>%
  rename(Zn = `Zn(g/t)`)

datos3 <- datos3 %>%
  rename(As = `As(g/t)`)


set.seed(123)

# Configuración de variables y fórmulas
configuracion_cu_zn_as <- list(
  "Zn" = list(
    formula = as.formula("log1p(Zn) ~ s(log1p(blue)) + s(log1p(green)) + s(log1p(red)) + s(log1p(nir08)) + s(log1p(coastal)) + s(log1p(swir16)) + s(log1p(swir22)) + s(log1p(rvi)) + s(log1p(sbi)) + s(log1p(ndwi))"),
    family = nb,
    mediana = 120.5
  ),
  "As" = list(
    formula = as.formula("log1p(As) ~ nir08 + s(swir16) + s(ndvi) + dvi + rvi + sbi"),
    family = nb,
    mediana = 20
  )
)

# Lista para almacenar resultados
resultados_cu_zn_as <- list()

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_r2 <- function(obs, pred) {
  if (sum(!is.na(obs) & !is.na(pred)) > 0) {
    return(cor(pred, obs, use = "complete.obs")^2)
  } else {
    return(NA)
  }
}

calcular_rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

calcular_bias <- function(obs, pred) {
  mean(pred - obs, na.rm = TRUE)
}

# Iterar por cada metal
for (metal in names(configuracion_cu_zn_as)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  config <- configuracion_cu_zn_as[[metal]]
  formula <- config$formula
  family <- config$family
  mediana <- config$mediana
  
  # Filtrar datos sin NA
  variables <- all.vars(formula)
  datos_filtrados <- na.omit(datos3[, variables])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variables[1]]], k = k, list = TRUE)
  
  # Listas para almacenar métricas
  r2_folds <- c()
  rmse_folds <- c()
  bias_folds <- c()
  
  for (i in 1:k) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    # Ajustar el modelo GAM
    modelo <- tryCatch({
      gam(formula, data = TRAIN, family = family())
    }, error = function(e) {
      message(sprintf("Error en el fold %d: %s", i, e$message))
      return(NULL)
    })
    
    if (is.null(modelo)) next
    
    # Predecir en el conjunto de validación
    pred <- tryCatch({
      predict(modelo, newdata = VALIDATION, type = "response")
    }, error = function(e) {
      message("Error en la predicción: ", e$message)
      return(rep(NA, nrow(VALIDATION)))
    })
    
    # Si hay transformación logarítmica en la variable de respuesta, deshacer la transformación
    if (grepl("log1p", as.character(formula)[2])) {
      pred <- expm1(pred)
    }
    
    # Calcular métricas
    r2 <- calcular_r2(VALIDATION[[variables[1]]], pred)
    rmse <- calcular_rmse(VALIDATION[[variables[1]]], pred)
    bias <- calcular_bias(VALIDATION[[variables[1]]], pred)
    
    # Estandarizar el RMSE por la mediana
    rmse_estandarizado <- rmse / mediana
    
    # Guardar métricas
    r2_folds <- c(r2_folds, r2)
    rmse_folds <- c(rmse_folds, rmse_estandarizado)
    bias_folds <- c(bias_folds, bias)
  }
  
  # Promedio de métricas
  r2_promedio <- mean(r2_folds, na.rm = TRUE)
  rmse_promedio <- mean(rmse_folds, na.rm = TRUE)
  bias_promedio <- mean(bias_folds, na.rm = TRUE)
  
  # Guardar resultados
  resultados_cu_zn_as[[metal]] <- list(
    modelo = modelo,
    R2 = r2_promedio,
    RMSE_estandarizado = rmse_promedio,
    Bias = bias_promedio
  )
  
  # Imprimir resultados para el metal actual
  cat("R2 promedio para", metal, ":", r2_promedio, "\n")
  cat("RMSE estandarizado promedio para", metal, ":", rmse_promedio, "\n")
  cat("Bias promedio para", metal, ":", bias_promedio, "\n")
}

# Imprimir resultados globales
for (metal in names(resultados_cu_zn_as)) {
  cat("\nResultados del modelo para:", metal, "\n")
  print(resultados_cu_zn_as[[metal]])
}



################################# Gráfico de residuos estandarizados

# Calcular residuos estandarizados del modelo final
for (metal in names(resultados_ni_co)) {
  cat("\nGráfico de residuos estandarizados para:", metal, "\n")
  
  # Obtener el modelo final ajustado
  modelo_final <- resultados_ni_co[[metal]]$modelo
  
  # Variables predictoras
  variables_predictoras <- setdiff(all.vars(configuracion_ni_co[[metal]]$formula), "Co")
  
  # Aplicar Brightness Normalization al conjunto completo
  datos_filtrados_normalizados <- brightness_normalization(datos_filtrados, variables_predictoras)
  
  # Calcular predicciones en el conjunto completo
  predicciones <- predict(modelo_final, newdata = datos_filtrados_normalizados, type = "response")
  
  # Calcular residuos estandarizados
  residuos <- datos_filtrados_normalizados$Co - predicciones
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Graficar residuos estandarizados
  plot(predicciones, residuos_estandarizados,
       main = paste("Cobalto"),
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       cex.lab = 1.5,
       cex.axis = 1.5,  
       cex.main = 2.5)
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(predicciones, residuos_estandarizados), col = "blue", lwd = 2)
}


################################################################################# Mejores modelos GLM y GAM


######################## Ni

set.seed(123)

# Configuración de variables y fórmula para Ni
configuracion_ni <- list(
  "Ni(g/t)" = list(
    formula = as.formula("`Ni(g/t)` ~ green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22"),
    family = gaussian(link = "identity")
  )
)

# Lista para almacenar resultados
resultados_ni <- list()

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- if (sum(!is.na(obs) & !is.na(pred)) > 0) cor(pred, obs, use = "complete.obs")^2 else NA
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  list(R2 = r2, RMSE = rmse, Bias = bias)
}

# Entrenamiento del modelo y cálculo de métricas
for (metal in names(configuracion_ni)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  config <- configuracion_ni[[metal]]
  formula <- config$formula
  family <- config$family
  
  # Filtrar datos sin NA
  variables <- all.vars(formula)
  datos_filtrados <- na.omit(datos3[, variables])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variables[1]]], k = k, list = TRUE)
  
  # Lista para almacenar métricas por fold
  metricas_folds <- list(R2 = c(), RMSE = c(), Bias = c())
  
  for (i in seq_along(folds)) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    # Ajustar el modelo
    modelo <- tryCatch({
      glm(formula, data = TRAIN, family = family)
    }, error = function(e) {
      message(sprintf("Error en el fold %d: %s", i, e$message))
      return(NULL)
    })
    
    if (is.null(modelo)) next
    
    # Predecir en el conjunto de validación
    pred <- predict(modelo, newdata = VALIDATION, type = "response")
    
    # Calcular métricas
    metricas <- calcular_metricas(VALIDATION[[variables[1]]], pred)
    metricas_folds$R2 <- c(metricas_folds$R2, metricas$R2)
    metricas_folds$RMSE <- c(metricas_folds$RMSE, metricas$RMSE)
    metricas_folds$Bias <- c(metricas_folds$Bias, metricas$Bias)
  }
  
  # Promedio de métricas
  metricas_promedio <- lapply(metricas_folds, mean, na.rm = TRUE)
  
  # Guardar resultados
  resultados_ni[[metal]] <- list(modelo = modelo, metricas = metricas_promedio)
  
  # Imprimir métricas promedio
  cat("Métricas promedio para", metal, ":\n")
  print(metricas_promedio)
}

# Gráfico de valores predichos vs observados con cuadro ajustado
graficar_pred_vs_real <- function(modelo, datos, metal, texto_personalizado) {
  pred <- predict(modelo, newdata = datos, type = "response")
  
  data_plot <- data.frame(
    Observado = datos[[metal]],
    Predicho = pred
  )
  
  max_val <- max(c(data_plot$Observado, data_plot$Predicho), na.rm = TRUE) * 1.1
  
  # Calcular posición del cuadro
  cuadro_x <- 0.05 * max_val
  cuadro_y <- 0.95 * max_val
  
  ggplot(data_plot, aes(x = Observado, y = Predicho)) +
    geom_point(alpha = 0.7, color = "#1f77b4", size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred", size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "green", linetype = "dashed", size = 1.2) +
    labs(
      x = "Valores Observados",
      y = "Valores Predichos"
    ) +
    coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 16, color = "black"),
      panel.grid.major = element_line(color = "gray", linetype = "dashed"),
      panel.grid.minor = element_blank()
    ) +
    # Agregar cuadro y texto personalizado
    annotate("rect", xmin = cuadro_x, xmax = cuadro_x + max_val * 0.36,
             ymin = cuadro_y - max_val * 0.26, ymax = cuadro_y*1.06,  
             fill = "white", alpha = 0.8, color = "black") +
    annotate("text", x = cuadro_x + max_val * 0.02, y = cuadro_y - max_val * 0.1,
             label = texto_personalizado, hjust = 0, size = 6.2, color = "black") +
    # Agregar simbología de líneas dentro del cuadro
    annotate("segment", x = cuadro_x + max_val * 0.02, xend = cuadro_x + max_val * 0.1,
             y = cuadro_y - max_val * 0.29, yend = cuadro_y - max_val * 0.29, 
             color = "darkred", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x + max_val * 0.12, y = cuadro_y - max_val * 0.29,
             label = "Ajuste ideal", hjust = 0, size = 5, color = "black") +
    annotate("segment", x = cuadro_x + max_val * 0.02, xend = cuadro_x + max_val * 0.1,
             y = cuadro_y - max_val * 0.34, yend = cuadro_y - max_val * 0.34, 
             color = "green", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x + max_val * 0.12, y = cuadro_y - max_val * 0.34,
             label = "Tendencia lineal", hjust = 0, size = 5, color = "black")
}

# Texto personalizado para Ni
texto_personalizado_ni <- "R2: 0.386\nRMSE: 0.399\nBias: -0.077"

# Generar gráfico para Ni
modelo_ni <- resultados_ni[["Ni(g/t)"]]$modelo
print(graficar_pred_vs_real(modelo_ni, datos_filtrados, "Ni(g/t)", texto_personalizado_ni))


########################## Cu, Zn y As

set.seed(123)

# Configuración de variables y fórmulas
configuracion_cu_zn_as <- list(
  "Cu(g/t)" = list(
    formula = as.formula("log1p(`Cu(g/t)`) ~ nir08 + swir16 + swir22 + evi + ndwi")
  ),
  "Zn(g/t)" = list(
    formula = as.formula("log1p(`Zn(g/t)`) ~ log1p(blue) + log1p(green) + log1p(red) + log1p(nir08) + log1p(coastal) + log1p(swir16) + log1p(swir22) + log1p(rvi) + log1p(sbi) + log1p(ndwi)")
  ),
  "As(g/t)" = list(
    formula = as.formula("log1p(`As(g/t)`) ~ nir08 + swir16 + ndvi + dvi + rvi + sbi")
  )
)

# Lista para almacenar resultados
resultados_cu_zn_as <- list()

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- if (sum(!is.na(obs) & !is.na(pred)) > 0) {
    cor(pred, obs, use = "complete.obs")^2
  } else {
    NA
  }
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Entrenamiento del modelo y cálculo de métricas
for (metal in names(configuracion_cu_zn_as)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  formula <- configuracion_cu_zn_as[[metal]]$formula
  
  # Filtrar datos sin NA
  variables <- all.vars(formula)
  datos_filtrados <- na.omit(datos3[, variables])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variables[1]]], k = k, list = TRUE)
  
  # Lista para almacenar métricas
  metricas_folds <- list(R2 = c(), RMSE = c(), Bias = c())
  
  for (i in seq_along(folds)) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    # Ajustar el modelo
    modelo <- tryCatch({
      glm.nb(formula, data = TRAIN)
    }, error = function(e) {
      message(sprintf("Error en el fold %d: %s", i, e$message))
      return(NULL)
    })
    
    if (is.null(modelo)) next
    
    # Predecir en el conjunto de validación
    pred <- tryCatch({
      predict(modelo, newdata = VALIDATION, type = "response")
    }, error = function(e) {
      message("Error en la predicción: ", e$message)
      return(rep(NA, nrow(VALIDATION)))
    })
    
    # Si hay transformación logarítmica en la variable de respuesta, deshacer la transformación
    if (grepl("log1p", as.character(formula)[2])) {
      pred <- expm1(pred)
    }
    
    # Calcular métricas
    metricas <- calcular_metricas(VALIDATION[[variables[1]]], pred)
    metricas_folds$R2 <- c(metricas_folds$R2, metricas$R2)
    metricas_folds$RMSE <- c(metricas_folds$RMSE, metricas$RMSE)
    metricas_folds$Bias <- c(metricas_folds$Bias, metricas$Bias)
  }
  
  # Promedio de métricas
  metricas_promedio <- lapply(metricas_folds, mean, na.rm = TRUE)
  
  # Guardar resultados
  resultados_cu_zn_as[[metal]] <- list(modelo = modelo, metricas_promedio = metricas_promedio)
  
  # Imprimir métricas promedio
  cat("Métricas promedio para", metal, ":\n")
  print(metricas_promedio)
}

# Gráfico con cuadro ajustado y simbología gráfica en la esquina superior derecha
graficar_pred_vs_real <- function(modelo, datos, metal, texto_personalizado) {
  # Predicciones
  pred <- predict(modelo, newdata = datos, type = "response")
  
  # Crear DataFrame para graficar
  data_plot <- data.frame(
    Real = datos[[metal]],
    Predicho = pred
  )
  
  # Calcular límites
  max_x <- max(data_plot$Real, na.rm = TRUE) * 1.1
  max_y <- max(data_plot$Predicho, na.rm = TRUE) * 1.1
  
  # Calcular posición del cuadro (esquina superior derecha)
  cuadro_x_min <- max_x * 0.6
  cuadro_x_max <- max_x * 0.97
  cuadro_y_min <- max_y * 0.01
  cuadro_y_max <- max_y * 0.33
  
  # Graficar
  ggplot(data_plot, aes(x = Real, y = Predicho)) +
    geom_point(alpha = 0.7, color = "#1f77b4", size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred", size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "green", linetype = "dashed", size = 1.2) +
    labs(
      x = "Valores Observados",
      y = "Valores Predichos"
    ) +
    coord_cartesian(xlim = c(0, max_x), ylim = c(0, max_y)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 16, color = "black"),
      panel.grid.major = element_line(color = "gray", linetype = "dashed"),
      panel.grid.minor = element_blank()
    ) +
    # Agregar cuadro y texto personalizado
    annotate("rect", xmin = cuadro_x_min * 1, xmax = cuadro_x_max * 1.03,
             ymin = cuadro_y_min * 14, ymax = cuadro_y_max* 1.28,
             fill = "white", alpha = 0.8, color = "black") +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             y = cuadro_y_max - (cuadro_y_max - cuadro_y_min) * 0.15,
             label = texto_personalizado, hjust = 0, size = 6.2, color = "black") +
    # Agregar simbología de líneas dentro del cuadro
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * 0.3, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * 0.3,
             color = "darkred", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * 0.3,
             label = "Ajuste ideal", hjust = 0, size = 5, color = "black") +
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * 0.1, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * 0.1,
             color = "green", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * 0.1,
             label = "Tendencia lineal", hjust = 0, size = 5, color = "black")
}

# Texto personalizado para cada metal
textos_personalizados <- list(
  "Zn(g/t)" = "R2: 0.075\nRMSE: 16.694\nBias: -521.227",
  "Cu(g/t)" = "R2: 0.076\nRMSE: 1.844\nBias: -1062.217",
  "As(g/t)" = "R2: 0.222\nRMSE: 17.280\nBias: -105.260"
)

# Generar gráficos para cada metal
for (metal in names(resultados_cu_zn_as)) {
  modelo <- resultados_cu_zn_as[[metal]]$modelo
  datos_filtrados <- na.omit(datos3[, all.vars(configuracion_cu_zn_as[[metal]]$formula)])
  
  # Generar texto personalizado
  texto_personalizado <- textos_personalizados[[metal]]
  
  
  print(graficar_pred_vs_real(modelo, datos_filtrados, metal, texto_personalizado))
}




############################################## Co con GAM

# Renombrar la columna
datos3 <- datos3 %>% rename(Co = `Co(g/t)`)

# Configuración del modelo
formula_co <- as.formula("Co ~ s(nir08) + s(swir16) + s(blue) + s(green) + s(coastal) + s(swir22)")
family_co <- Gamma(link = "log")

# Aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, vars] <- df[, vars] / brightness
  return(df)
}

# Filtrar y normalizar los datos
variables_predictoras <- c("nir08", "swir16", "blue", "green", "coastal", "swir22")
datos_filtrados_co <- na.omit(datos3[, c("Co", variables_predictoras)])
datos_filtrados_co <- brightness_normalization(datos_filtrados_co, variables_predictoras)

# Ajustar el modelo GAM
modelo_co <- gam(formula_co, data = datos_filtrados_co, family = family_co)

# Calcular métricas promedio (opcional)
r2 <- cor(predict(modelo_co, type = "response"), datos_filtrados_co$Co)^2
rmse <- sqrt(mean((predict(modelo_co, type = "response") - datos_filtrados_co$Co)^2))
bias <- mean(predict(modelo_co, type = "response") - datos_filtrados_co$Co)
texto_personalizado_co <- "R2: 0.214\nRMSE: 1.050\nBias: -0.258"

# Función para graficar
graficar_pred_vs_real <- function(modelo, datos, metal, texto_personalizado) {
  pred <- predict(modelo, newdata = datos, type = "response")
  
  data_plot <- data.frame(
    Real = datos[[metal]],
    Predicho = pred
  )
  
  max_val <- max(c(data_plot$Real, data_plot$Predicho), na.rm = TRUE) * 1.1
  
  # Calcular posición del cuadro
  cuadro_x_min <- max_val * 0.05
  cuadro_x_max <- max_val * 0.3
  cuadro_y_min <- max_val * 0.8
  cuadro_y_max <- max_val * 0.9
  
  ggplot(data_plot, aes(x = Real, y = Predicho)) +
    geom_point(alpha = 0.7, color = "#1f77b4", size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred", size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "green", linetype = "dashed", size = 1.2) +
    labs(
      
      x = "Valores Observados",
      y = "Valores Predichos"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 16, color = "black"),
      panel.grid.major = element_line(color = "gray", linetype = "dashed"),
      panel.grid.minor = element_blank()
    ) +
    # Agregar cuadro y texto personalizado
    annotate("rect", xmin = cuadro_x_min, xmax = cuadro_x_max *1.3,
             ymin = cuadro_y_min * 0.91, ymax = cuadro_y_max * 1.15,
             fill = "white", alpha = 0.8, color = "black") +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.1, 
             y = cuadro_y_max - (cuadro_y_max - cuadro_y_min) * 0.2,
             label = texto_personalizado, hjust = 0, size = 6.2, color = "black") +
    # Agregar simbología de líneas dentro del cuadro
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1,
             color = "darkred", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1,
             label = "Ajuste ideal", hjust = 0, size = 5, color = "black") +
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1.6, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1.6,
             color = "green", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1.6,
             label = "Tendencia lineal", hjust = 0, size = 5, color = "black")
}

# Generar gráfico para Co
print(graficar_pred_vs_real(modelo_co, datos_filtrados_co, "Co", texto_personalizado_co))

########################## Análisis de Sensibilidad

# Análisis de sensibilidad para modelos GLM y GAM
analizar_sensibilidad <- function(modelo, metal) {
  # Extraer los coeficientes del modelo
  coeficientes <- tryCatch({
    coef(modelo)
  }, error = function(e) {
    cat("\nNo se pudieron extraer los coeficientes para:", metal, "\n")
    return(NULL)
  })
  
  if (!is.null(coeficientes)) {
    # Convertir coeficientes en un DataFrame
    coeficientes_df <- as.data.frame(coeficientes)
    coeficientes_df$Variable <- rownames(coeficientes_df)
    colnames(coeficientes_df) <- c("Importancia", "Variable")
    
    # Excluir el intercepto (si está presente)
    coeficientes_df <- coeficientes_df[coeficientes_df$Variable != "(Intercept)", ]
    
    # Ordenar por la importancia absoluta de los coeficientes
    coeficientes_df <- coeficientes_df[order(abs(coeficientes_df$Importancia), decreasing = TRUE), ]
    
    # Graficar
    ggplot(coeficientes_df, aes(x = reorder(Variable, abs(Importancia)), y = Importancia)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        x = "Variable Predictora",
        y = "Coeficiente (Importancia)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16, color = "black"),  
        axis.text = element_text(size = 14, color = "black")  
      )
  }
}


# Crear tabla de análisis de sensibilidad
crear_tabla_sensibilidad <- function(modelo, metal) {
  # Extraer los coeficientes del modelo
  coeficientes <- tryCatch({
    coef(modelo)
  }, error = function(e) {
    cat("\nNo se pudieron extraer los coeficientes para:", metal, "\n")
    return(NULL)
  })
  
  if (!is.null(coeficientes)) {
    # Convertir coeficientes en un DataFrame
    coeficientes_df <- as.data.frame(coeficientes)
    coeficientes_df$Variable <- rownames(coeficientes_df)
    colnames(coeficientes_df) <- c("Coeficiente", "Variable")
    
    # Excluir el intercepto (si está presente)
    coeficientes_df <- coeficientes_df[coeficientes_df$Variable != "(Intercept)", ]
    
    # Ordenar por la importancia absoluta de los coeficientes
    coeficientes_df <- coeficientes_df[order(abs(coeficientes_df$Coeficiente), decreasing = TRUE), ]
    
    # Imprimir la tabla
    cat("\nTabla de Análisis de Sensibilidad para:", metal, "\n")
    print(coeficientes_df)
    
    # Retornar la tabla
    return(coeficientes_df)
  }
  
  return(NULL)
}


# Modelo y datos para Cu
modelo_cu <- resultados_cu_zn_as[["Cu(g/t)"]]$modelo
analizar_sensibilidad(modelo_cu, "Cu")
tabla_cu <- crear_tabla_sensibilidad(modelo_cu, "Cu")


# Modelo y datos para Zn
modelo_zn <- resultados_cu_zn_as[["Zn(g/t)"]]$modelo
analizar_sensibilidad(modelo_zn, "Zn")
tabla_zn <- crear_tabla_sensibilidad(modelo_zn, "Zn")


# Modelo y datos para As
modelo_as <- resultados_cu_zn_as[["As(g/t)"]]$modelo
analizar_sensibilidad(modelo_as, "As")
tabla_as <- crear_tabla_sensibilidad(modelo_as, "As")





