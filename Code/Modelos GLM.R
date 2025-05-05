######################## Librerías generales
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

median(datos3$`As(g/t)`)
summary(datos3$`Ni(g/t)`)
# Tratar la negatividad de las variables predictoras

variables_predictoras <- c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", 
                           "green", "red", "nir08", "coastal", "swir16", "swir22")

# Sumar 1 a todas las columnas de las variables predictoras
datos3[variables_predictoras] <- datos3[variables_predictoras] + 1

############################################################### 2.2.1) Análisis de Heterocedasticidad 

library(lmtest)

# Para cada uno de los metales se aplican las combinaciones de variables de a continuación para el análisis:
# Ni : green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22
# Cu: nir08 + swir16 + swir22 + evi + ndwi
# Zn: blue + green + red + nir08 + coastal + swir16 + swir22 + rvi + sbi + ndwi
# Co: nir08 + swir16 + blue + green + coastal + swir22
# As: nir08 + swir16 + ndvi + dvi + rvi + sbi

############################################ Gaussian(identity)

glm_gaussian <- gam(`Ni(g/t)` ~ s(blue) + green + red + nir08 + coastal + swir16 + swir22, data = datos3, 
                    family = gaussian(link = "log"))
bptest(glm_gaussian)
residuos_glm <- rstandard(glm_gaussian)


plot(glm_gaussian$fitted.values, residuos_glm,
     main = "Gaussian(identity)",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.main = 2.5,
     cex.lab = 1.5)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(glm_gaussian$fitted.values, residuos_glm), col = "blue", lwd = 2)

############################################ Gamma(link=log)

glm_gamma <- glm(`Co(g/t)` ~ blue + green + red + nir08 + coastal + swir16 + swir22, data = datos3, family = Gamma(link = "log"))
bptest(glm_gamma)
residuos_glm <- rstandard(glm_gamma)


plot(glm_gamma$fitted.values, residuos_glm,
     main = "Gamma(link = log)",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(glm_gamma$fitted.values, residuos_glm), col = "blue", lwd = 2)

############################################ Gaussian(link=log)

glm_gaussian2 <- glm(`Ni(g/t)` ~ blue + green + red + nir08 + coastal + swir16 + swir22, data = datos3, 
                     family = gaussian(link = "identity"))
bptest(glm_gaussian2)
residuos_glm <- rstandard(glm_gaussian2)


plot(glm_gaussian2$fitted.values, residuos_glm,
     main = "Gaussian(link = identity)",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(glm_gaussian2$fitted.values, residuos_glm), col = "blue", lwd = 2)


############################################ Binomial negativa

glm_binomial_neg <- glm.nb(`As(g/t)` ~ blue + green + red + nir08 + coastal + swir16 + swir22, data = datos3)

# Breusch-Pagan Test para evaluar homocedasticidad
bptest(glm_binomial_neg)

# Calcular residuos estandarizados
residuos_glm <- rstandard(glm_binomial_neg)

# Gráfico de residuos vs valores ajustados
plot(glm_binomial_neg$fitted.values, residuos_glm,
     main = "Binomial Negativa",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(glm_binomial_neg$fitted.values, residuos_glm), col = "blue", lwd = 2)

qqnorm(residuos_glm)
qqline(residuos_glm, col = "blue")


# AIC/BIC

modelos <- c("Gaussian (identity)", "Gamma (log)", "Gaussian (log)", "Negative Binomial")
aic <- c(AIC(glm_gaussian), AIC(glm_gamma), AIC(glm_gaussian2), AIC(glm_binomial_neg))
bic <- c(BIC(glm_gaussian), BIC(glm_gamma), BIC(glm_gaussian2), BIC(glm_binomial_neg))

comparacion_modelos <- data.frame(Modelo = modelos, AIC = aic, BIC = bic)
print(comparacion_modelos)




############################################################### 2.2) GLM con k-fold = 10

library(caret)

################################# Ni y Co

set.seed(123)

k <- 10

# Función para generar todas las combinaciones de las variables adicionales
generar_combinaciones <- function(variables_adicionales) {
  combinaciones <- list()
  for (n in 1:length(variables_adicionales)) {
    combinaciones <- c(combinaciones, combn(variables_adicionales, n, simplify = FALSE))
  }
  return(combinaciones)
}

# Variables, combinaciones base y familias por cada variable de respuesta
variables_por_respuesta <- list(
  "Ni(g/t)" = list(
    base = c("green", "nir08"),
    adicionales = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "red", "coastal", "swir16", "swir22"),
    extra_comb = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22"),
    family = gaussian(link = "identity")
  ),
  "Co(g/t)" = list(
    base = c("nir08", "swir16"),
    adicionales = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "green", "red", "coastal", "swir22"),
    extra_comb = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22"),
    family = Gamma(link = "log")
    )
)


# Almacenar resultados para cada variable
resultados_globales <- list()

# Iterar sobre cada variable de respuesta
for (variable_respuesta in names(variables_por_respuesta)) {
  
  datos_var <- variables_por_respuesta[[variable_respuesta]]
  variables_base <- datos_var$base
  variables_adicionales <- datos_var$adicionales
  extra_comb <- datos_var$extra_comb
  family_glm <- datos_var$family
  
  # Generar combinaciones dinámicamente
  combinaciones_adicionales <- generar_combinaciones(variables_adicionales)
  combinaciones <- lapply(combinaciones_adicionales, function(x) c(variables_base, x))
  
  # Añadir la combinación fija adicional
  combinaciones <- c(combinaciones, list(extra_comb))
  
  # Cantidad de combinaciones probadas
  cantidad_combinaciones <- length(combinaciones)
  
  # Validación cruzada
  folds <- createFolds(datos3[[variable_respuesta]], k = k)
  
  # Almacenar métricas para la mejor combinación
  mejor_combinacion <- NULL
  mejor_r2 <- -Inf  # Inicializar con un valor bajo
  mejor_rmse <- NA
  mejor_bias <- NA
  
  for (variables in combinaciones) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Crear fórmula dinámica para el GLM
      formula <- as.formula(paste0("`", variable_respuesta, "` ~ ", paste(variables, collapse = " + ")))
      
      # Ajustar el modelo GLM con manejo de errores
      modelo_glm <- tryCatch(
        glm(formula, data = TRAIN, family = family_glm),
        error = function(e) {
          message(sprintf("Error en el ajuste GLM para la combinación: %s. Error: %s", 
                          paste(variables, collapse = ", "), e$message))
          return(NULL)
        }
      )
      
      if (is.null(modelo_glm)) next  # Saltar en caso de error
      
      # Predecir en el conjunto de validación
      pred <- tryCatch(
        predict(modelo_glm, newdata = VALIDATION, type = "response"),
        error = function(e) {
          message("Error en la predicción: ", e$message)
          return(rep(NA, nrow(VALIDATION)))  # Llenar con NA en caso de error
        }
      )
      obs <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas
      if (!all(is.na(pred))) {
        r2[i] <- cor(pred, obs, use = "complete.obs")^2
        rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        bias[i] <- mean(pred - obs, na.rm = TRUE)
      } else {
        r2[i] <- NA
        rmse[i] <- NA
        bias[i] <- NA
      }
    }
    
    # Promedio de métricas para esta combinación
    r2_promedio <- mean(r2, na.rm = TRUE)
    rmse_promedio <- mean(rmse, na.rm = TRUE)
    bias_promedio <- mean(bias, na.rm = TRUE)
    
    # Actualizar la mejor combinación
    if (!is.na(r2_promedio) && r2_promedio > mejor_r2) {
      mejor_r2 <- r2_promedio
      mejor_rmse <- rmse_promedio
      mejor_bias <- bias_promedio
      mejor_combinacion <- variables
    }
  }
  
  # Guardar resultados
  resultados_globales[[variable_respuesta]] <- list(
    Mejor_Combinacion = mejor_combinacion,
    R2 = mejor_r2,
    RMSE = mejor_rmse,
    Bias = mejor_bias,
    Cantidad_Combinaciones = cantidad_combinaciones
  )
}

# Imprimir resultados para cada variable
for (variable_respuesta in names(resultados_globales)) {
  cat("\nResultados para:", variable_respuesta, "\n")
  print(resultados_globales[[variable_respuesta]])
}

summary(datos3)



########################## Cu, Zn, As  (el del As se demora 20 min)


set.seed(123)
k <- 10

# Configuración de variables base y adicionales por metal
variables_por_respuesta <- list(
  "Cu(g/t)" = list(
    base = c("nir08", "swir16", "swir22"),
    adicionales = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "green", "red", "coastal")
  ),
  "Zn(g/t)" = list(
    base = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22"),
    adicionales = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi")
  ),
  "As(g/t)" = list(
    base = c("nir08"),
    adicionales = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "green", "red", "coastal", "swir16", "swir22")
  )
)



# Función para generar combinaciones de variables adicionales
generar_combinaciones <- function(variables_adicionales) {
  combinaciones <- list()
  for (n in 1:length(variables_adicionales)) {
    combinaciones <- c(combinaciones, combn(variables_adicionales, n, simplify = FALSE))
  }
  return(combinaciones)
}

# Almacenar resultados globales
resultados_globales <- list()

# Iterar por cada metal
for (variable_respuesta in names(variables_por_respuesta)) {
  
  datos_var <- variables_por_respuesta[[variable_respuesta]]
  variables_base <- datos_var$base
  variables_adicionales <- datos_var$adicionales
  
  # Generar combinaciones dinámicas
  combinaciones_adicionales <- generar_combinaciones(variables_adicionales)
  combinaciones <- lapply(combinaciones_adicionales, function(x) c(variables_base, x))
  
  # Contar la cantidad de combinaciones probadas
  cantidad_combinaciones <- length(combinaciones)
  
  # Validación cruzada
  folds <- createFolds(datos3[[variable_respuesta]], k = k)
  
  # Variables para almacenar los mejores resultados
  mejor_combinacion <- NULL
  mejor_r2 <- -Inf
  mejor_rmse <- NA
  mejor_bias <- NA
  particiones_validas <- 0  # Contador para particiones válidas
  
  for (variables in combinaciones) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    particiones_validas_combinacion <- 0  # Contador por combinación
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Aplicar transformación logarítmica a la variable de respuesta
      TRAIN[[variable_respuesta]] <- log1p(TRAIN[[variable_respuesta]])  # log(1 + Y)
      
      # Crear fórmula dinámica
      formula <- as.formula(paste0("`", variable_respuesta, "` ~ ", paste(variables, collapse = " + ")))
      
      # Ajustar el modelo glm.nb con manejo de errores
      modelo_glm <- tryCatch(
        glm.nb(formula, data = TRAIN),
        error = function(e) {
          message(sprintf("Error en el ajuste GLM para la combinación: %s. Error: %s", 
                          paste(variables, collapse = ", "), e$message))
          return(NULL)
        }
      )
      
      if (is.null(modelo_glm)) next  # Saltar en caso de error
      
      # Predecir en el conjunto de validación (sin transformación logarítmica)
      pred <- tryCatch(
        expm1(predict(modelo_glm, newdata = VALIDATION, type = "response")),  # Volver a la escala original
        error = function(e) {
          message("Error en la predicción: ", e$message)
          return(rep(NA, nrow(VALIDATION)))
        }
      )
      obs <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas
      if (!all(is.na(pred))) {
        r2[i] <- cor(pred, obs, use = "complete.obs")^2
        rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        bias[i] <- mean(pred - obs, na.rm = TRUE)  # Sesgo promedio
        particiones_validas_combinacion <- particiones_validas_combinacion + 1
      } else {
        r2[i] <- NA
        rmse[i] <- NA
        bias[i] <- NA
      }
    }
    
    # Promedio de métricas para esta combinación
    r2_promedio <- mean(r2, na.rm = TRUE)
    rmse_promedio <- mean(rmse, na.rm = TRUE)
    bias_promedio <- mean(bias, na.rm = TRUE)
    
    # Actualizar la mejor combinación
    if (!is.na(r2_promedio) && r2_promedio > mejor_r2) {
      mejor_r2 <- r2_promedio
      mejor_rmse <- rmse_promedio
      mejor_bias <- bias_promedio
      mejor_combinacion <- variables
      particiones_validas <- particiones_validas_combinacion  # Actualizar el número de particiones válidas
    }
  }
  
  # Guardar los mejores resultados en la lista global
  resultados_globales[[variable_respuesta]] <- list(
    Mejor_Combinacion = mejor_combinacion,
    R2 = mejor_r2,
    RMSE = mejor_rmse,
    Bias = mejor_bias,
    Cantidad_Combinaciones = cantidad_combinaciones,
    Particiones_Validas = particiones_validas
  )
}

# Imprimir resultados para cada variable
for (variable_respuesta in names(resultados_globales)) {
  cat("\nResultados para:", variable_respuesta, "\n")
  print(resultados_globales[[variable_respuesta]])
}



############################################################### 2.3) Gráficos de residuos


# Para cada uno de los metales se aplican las combinaciones de variables de a continuación para graficar:
# Ni : green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22
# Cu: nir08 + swir16 + swir22 + evi + ndwi
# Zn: blue + green + red + nir08 + coastal + swir16 + swir22 + rvi + sbi + ndwi
# Co: nir08 + swir16 + blue + green + coastal + swir22
# As: nir08 + swir16 + ndvi + dvi + rvi + sbi


# Configuración de las fórmulas y familias para Ni y Co
metales_ni_co <- list(
  "Ni(g/t)" = list(  # Usar el nombre exacto del dataset
    formula = as.formula("`Ni(g/t)` ~ green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22"),
    family = gaussian(link = "identity")
  ),
  "Co(g/t)" = list(  # Usar el nombre exacto del dataset
    formula = as.formula("`Co(g/t)` ~ nir08 + swir16 + blue + green + coastal + swir22"),
    family = Gamma(link = "log")
  )
)

# Iterar por cada metal
for (metal in names(metales_ni_co)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(all.vars(metales_ni_co[[metal]]$formula), metal)])
  
  # Ajustar el modelo GLM
  formula <- metales_ni_co[[metal]]$formula
  family <- metales_ni_co[[metal]]$family
  modelo <- glm(formula, data = datos_filtrados, family = family)
  
  # Predicciones y residuos
  predicciones <- predict(modelo, type = "response")
  residuos <- datos_filtrados[[metal]] - predicciones
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Gráfico de residuos estandarizados
  plot(predicciones, residuos_estandarizados,
       main = "Cobalto",
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 2.5)
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(predicciones, residuos_estandarizados), col = "blue", lwd = 2)
}

#################### Cu, Zn y As


# Configuración de las fórmulas para Cu, Zn y As
metales_cu_zn_as <- list(
  "Cu(g/t)" = list(
    formula = as.formula("log1p(`Cu(g/t)`) ~ nir08 + swir16 + swir22 + evi + ndwi")
  ),
  "Zn(g/t)" = list(
    formula = as.formula("log1p(`Zn(g/t)`) ~ blue + green + red + nir08 + coastal + swir16 + swir22 + rvi + sbi + ndwi")
  ),
  "As(g/t)" = list(
    formula = as.formula("log1p(`As(g/t)`) ~ nir08 + swir16 + ndvi + dvi + rvi + sbi")
  )
)

# Iterar por cada metal
for (metal in names(metales_cu_zn_as)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Ajustar el modelo GLM
  formula <- metales_cu_zn_as[[metal]]$formula
  modelo <- glm.nb(formula, data = datos3)
  
  # Predicciones en escala original
  predicciones_log <- predict(modelo, type = "response")  # Predicciones en escala logarítmica
  predicciones <- expm1(predicciones_log)  # Transformar a escala original
  residuos <- datos3[[metal]] - predicciones
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Gráfico de residuos estandarizados
  plot(predicciones, residuos_estandarizados,
       main = "Cobre",
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       cex.lab = 1.5,
       cex.axis = 1.5,  
       cex.main = 2.5)
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(predicciones, residuos_estandarizados), col = "blue", lwd = 2)

}



################################################################## Experimento 2 - GLM con transformaciones en Y

####################### Ni y Co


set.seed(123)

# Configuración de variables y familias definidas
configuracion_metales <- list(
  "Ni(g/t)" = list(
    formula = as.formula("`Ni(g/t)` ~ green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22"),
    family = gaussian(link = "identity")
  ),
  "Co(g/t)" = list(
    formula = as.formula("`Co(g/t)` ~ nir08 + swir16 + blue + green + coastal + swir22"),
    family = Gamma(link = "log")
  )
)

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- cor(pred, obs, use = "complete.obs")^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Lista para almacenar resultados
resultados <- list()

# Número de folds para cross-validation
k <- 10
epsilon <- 0  # Constante para la transformación logarítmica (log(Y + 1))

# Iterar por cada metal
for (metal in names(configuracion_metales)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  formula <- configuracion_metales[[metal]]$formula
  family <- configuracion_metales[[metal]]$family
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(all.vars(formula), metal)])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[metal]], k = k, list = TRUE)
  
  # Lista para almacenar métricas de transformaciones
  metricas_transformaciones <- list(Logarítmica = list())
  
  ## ---- Transformación Logarítmica (log(Y + epsilon)) ----
  r2_log <- rmse_log <- bias_log <- c()
  particiones_validas_log <- 0
  
  for (i in 1:k) {
    train <- datos_filtrados[-folds[[i]], ]
    test <- datos_filtrados[folds[[i]], ]
    
    # Transformación logarítmica
    train$Y_log <- log1p(train[[metal]])  # log(1 + Y)
    
    modelo_log <- tryCatch(
      glm(formula, data = train[, c(all.vars(formula), "Y_log")], family = family),
      error = function(e) NULL
    )
    
    if (!is.null(modelo_log)) {
      # Volver a la escala original con seguridad para valores no negativos
      pred_log <- tryCatch(
        pmax(expm1(predict(modelo_log, newdata = test, type = "response")), 0),
        error = function(e) rep(NA, nrow(test))
      )
      obs_log <- test[[metal]]
      
      if (!all(is.na(pred_log))) {
        metricas <- calcular_metricas(obs_log, pred_log)
        r2_log <- c(r2_log, metricas$R2)
        rmse_log <- c(rmse_log, metricas$RMSE)
        bias_log <- c(bias_log, metricas$Bias)
        particiones_validas_log <- particiones_validas_log + 1
      }
    }
  }
  
  metricas_transformaciones$Logarítmica <- list(
    R2 = mean(r2_log, na.rm = TRUE),
    RMSE = mean(rmse_log, na.rm = TRUE),
    Bias = mean(bias_log, na.rm = TRUE),
    Particiones_Validas = particiones_validas_log
  )
  
  # Guardar resultados
  resultados[[metal]] <- metricas_transformaciones
}

# Imprimir resultados
for (metal in names(resultados)) {
  cat("\nResultados para:", metal, "\n")
  print(resultados[[metal]])
}


####################### Cu, Zn y As

set.seed(123)

# Configuración de variables y combinación definida
variables_por_respuesta <- list(
  "Cu(g/t)" = c("nir08", "swir16", "swir22", "evi", "ndwi"),
  "Zn(g/t)" = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi", "ndwi"),
  "As(g/t)" = c("nir08", "swir16", "ndvi", "dvi", "rvi", "sbi")
)

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- cor(pred, obs, use = "complete.obs")^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Lista para almacenar resultados
resultados_globales <- list()

# Iterar por cada metal con una combinación definida
for (variable_respuesta in names(variables_por_respuesta)) {
  cat("\nAnálisis para:", variable_respuesta, "\n")
  
  variables <- variables_por_respuesta[[variable_respuesta]]
  formula <- as.formula(paste0("`", variable_respuesta, "` ~ ", paste(variables, collapse = " + ")))
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(variables, variable_respuesta)])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variable_respuesta]], k = k, list = TRUE)
  
  # Almacenar métricas para cada transformación
  metricas_log <- list(R2 = c(), RMSE = c(), Bias = c())
  metricas_quad <- list(R2 = c(), RMSE = c(), Bias = c())
  particiones_validas_log <- 0
  particiones_validas_quad <- 0
  
  for (i in 1:k) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    ## ---- Transformación Logarítmica ----
    TRAIN$Y_log <- log1p(TRAIN[[variable_respuesta]])  # log(1 + Y)
    
    modelo_glm_log <- tryCatch(
      glm.nb(Y_log ~ ., data = TRAIN[, c(variables, "Y_log")]),
      error = function(e) {
        message(sprintf("Error en el ajuste GLM para la transformación logarítmica en el fold %d: %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(modelo_glm_log)) {
      pred_log <- tryCatch(
        expm1(predict(modelo_glm_log, newdata = VALIDATION[, variables], type = "response")),
        error = function(e) {
          message("Error en la predicción logarítmica: ", e$message)
          return(rep(NA, nrow(VALIDATION)))
        }
      )
      obs_log <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas si las predicciones son válidas
      if (!all(is.na(pred_log))) {
        metricas <- calcular_metricas(obs_log, pred_log)
        metricas_log$R2 <- c(metricas_log$R2, metricas$R2)
        metricas_log$RMSE <- c(metricas_log$RMSE, metricas$RMSE)
        metricas_log$Bias <- c(metricas_log$Bias, metricas$Bias)
        particiones_validas_log <- particiones_validas_log + 1
      }
    }
    
    ## ---- Transformación Cuadrática ----
    TRAIN$Y_quad <- sqrt(TRAIN[[variable_respuesta]])  # sqrt(Y)
    
    modelo_glm_quad <- tryCatch(
      glm.nb(Y_quad ~ ., data = TRAIN[, c(variables, "Y_quad")]),
      error = function(e) {
        message(sprintf("Error en el ajuste GLM para la transformación cuadrática en el fold %d: %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(modelo_glm_quad)) {
      pred_quad <- tryCatch(
        predict(modelo_glm_quad, newdata = VALIDATION[, variables], type = "response")^2,
        error = function(e) {
          message("Error en la predicción cuadrática: ", e$message)
          return(rep(NA, nrow(VALIDATION)))
        }
      )
      obs_quad <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas si las predicciones son válidas
      if (!all(is.na(pred_quad))) {
        metricas <- calcular_metricas(obs_quad, pred_quad)
        metricas_quad$R2 <- c(metricas_quad$R2, metricas$R2)
        metricas_quad$RMSE <- c(metricas_quad$RMSE, metricas$RMSE)
        metricas_quad$Bias <- c(metricas_quad$Bias, metricas$Bias)
        particiones_validas_quad <- particiones_validas_quad + 1
      }
    }
  }
  
  # Promedio de métricas para cada transformación
  resultados_globales[[variable_respuesta]] <- list(
    Logarítmica = list(
      R2 = mean(metricas_log$R2, na.rm = TRUE),
      RMSE = mean(metricas_log$RMSE, na.rm = TRUE),
      Bias = mean(metricas_log$Bias, na.rm = TRUE),
      Particiones_Validas = particiones_validas_log
    ),
    Cuadrática = list(
      R2 = mean(metricas_quad$R2, na.rm = TRUE),
      RMSE = mean(metricas_quad$RMSE, na.rm = TRUE),
      Bias = mean(metricas_quad$Bias, na.rm = TRUE),
      Particiones_Validas = particiones_validas_quad
    )
  )
}

# Imprimir resultados para cada variable de respuesta
for (variable_respuesta in names(resultados_globales)) {
  cat("\nResultados para:", variable_respuesta, "\n")
  print(resultados_globales[[variable_respuesta]])
}

##################################### Grafico de residuos estandarizados


# Configuración de las fórmulas para Cu, Zn y As
metales_cu_zn_as <- list(
  "Zn(g/t)" = list(
    formula = as.formula("sqrt(`Zn(g/t)`) ~ blue + green + red + nir08 + coastal + swir16 + swir22 + rvi + sbi + ndwi")
  )
)

# Iterar por cada metal
for (metal in names(metales_cu_zn_as)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Ajustar el modelo GLM
  formula <- metales_cu_zn_as[[metal]]$formula
  modelo <- glm.nb(formula, data = datos3)
  
  # Predicciones en escala cuadrática
  predicciones_sqrt <- predict(modelo, type = "response")  # Predicciones en escala transformada
  predicciones <- predicciones_sqrt^2  # Transformar de vuelta a la escala original
  residuos <- datos3[[metal]] - predicciones
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Gráfico de residuos estandarizados
  plot(predicciones, residuos_estandarizados,
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       cex.lab = 1.5,
       main = paste("Residuos Estandarizados para", metal))
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(predicciones, residuos_estandarizados), col = "blue", lwd = 2)
  
}



############################################################### 4) Experimento 3

############################################################### 4.1) GLM con k-fold y Transformaciones en X


#################### Ni y Co


set.seed(123)

# Configuración de variables y familias definidas
configuracion_metales <- list(
  "Ni(g/t)" = list(
    formula = as.formula("`Ni(g/t)` ~ green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22"),
    family = gaussian(link = "identity")
  ),
  "Co(g/t)" = list(
    formula = as.formula("`Co(g/t)` ~ nir08 + swir16 + blue + green + coastal + swir22"),
    family = Gamma(link = "log")
  )
)

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- cor(pred, obs, use = "complete.obs")^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Función de Brightness Normalization
brightness_normalization <- function(df, variables) {
  brightness <- sqrt(rowSums(df[, variables]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, variables] <- df[, variables] / brightness
  return(df)
}

# Función de Standard Normalization
standard_normalization <- function(df, variables) {
  df[, variables] <- scale(df[, variables])
  return(df)
}

# Función de Derivadas Espectrales
derivadas_espectrales <- function(df, variables) {
  for (i in 2:length(variables)) {
    df[[paste0("deriv_", variables[i])]] <- df[[variables[i]]] - df[[variables[i - 1]]]
  }
  return(df)
}

# Lista para almacenar resultados
resultados <- list()

# Número de folds para cross-validation
k <- 10

# Iterar por cada metal
for (metal in names(configuracion_metales)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  formula <- configuracion_metales[[metal]]$formula
  family <- configuracion_metales[[metal]]$family
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(all.vars(formula), metal)])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[metal]], k = k, list = TRUE)
  
  # Almacenar métricas para cada transformación
  metricas_transformaciones <- list(
    Logarítmica = list(),
    Cuadrática = list(),
    Brightness_Normalization = list(),
    Standard_Normalization = list(),
    Derivadas_Espectrales = list()
  )
  
  # Variables predictoras
  variables_predictoras <- all.vars(formula)[-1]
  
  ## ---- 1. Transformación Logarítmica ----
  r2_log <- rmse_log <- bias_log <- c()
  for (i in 1:k) {
    train <- datos_filtrados[-folds[[i]], ]
    test <- datos_filtrados[folds[[i]], ]
    
    # Aplicar transformación logarítmica
    train[, variables_predictoras] <- log1p(train[, variables_predictoras])
    test[, variables_predictoras] <- log1p(test[, variables_predictoras])
    
    modelo_log <- tryCatch(glm(formula, data = train, family = family), error = function(e) NULL)
    
    if (!is.null(modelo_log)) {
      pred_log <- predict(modelo_log, newdata = test, type = "response")
      obs_log <- test[[metal]]
      metricas <- calcular_metricas(obs_log, pred_log)
      r2_log <- c(r2_log, metricas$R2)
      rmse_log <- c(rmse_log, metricas$RMSE)
      bias_log <- c(bias_log, metricas$Bias)
    }
  }
  metricas_transformaciones$Logarítmica <- list(
    R2 = mean(r2_log, na.rm = TRUE),
    RMSE = mean(rmse_log, na.rm = TRUE),
    Bias = mean(bias_log, na.rm = TRUE)
  )
  
  ## ---- 2. Transformación Cuadrática ----
  r2_quad <- rmse_quad <- bias_quad <- c()
  for (i in 1:k) {
    train <- datos_filtrados[-folds[[i]], ]
    test <- datos_filtrados[folds[[i]], ]
    
    # Aplicar transformación cuadrática
    train[, variables_predictoras] <- sqrt(train[, variables_predictoras])
    test[, variables_predictoras] <- sqrt(test[, variables_predictoras])
    
    modelo_quad <- tryCatch(glm(formula, data = train, family = family), error = function(e) NULL)
    
    if (!is.null(modelo_quad)) {
      pred_quad <- predict(modelo_quad, newdata = test, type = "response")
      obs_quad <- test[[metal]]
      metricas <- calcular_metricas(obs_quad, pred_quad)
      r2_quad <- c(r2_quad, metricas$R2)
      rmse_quad <- c(rmse_quad, metricas$RMSE)
      bias_quad <- c(bias_quad, metricas$Bias)
    }
  }
  metricas_transformaciones$Cuadrática <- list(
    R2 = mean(r2_quad, na.rm = TRUE),
    RMSE = mean(rmse_quad, na.rm = TRUE),
    Bias = mean(bias_quad, na.rm = TRUE)
  )
  
  ## ---- 3. Brightness Normalization ----
  r2_bright <- rmse_bright <- bias_bright <- c()
  for (i in 1:k) {
    train <- brightness_normalization(datos_filtrados[-folds[[i]], ], variables_predictoras)
    test <- brightness_normalization(datos_filtrados[folds[[i]], ], variables_predictoras)
    
    modelo_bright <- tryCatch(glm(formula, data = train, family = family), error = function(e) NULL)
    
    if (!is.null(modelo_bright)) {
      pred_bright <- predict(modelo_bright, newdata = test, type = "response")
      obs_bright <- test[[metal]]
      metricas <- calcular_metricas(obs_bright, pred_bright)
      r2_bright <- c(r2_bright, metricas$R2)
      rmse_bright <- c(rmse_bright, metricas$RMSE)
      bias_bright <- c(bias_bright, metricas$Bias)
    }
  }
  metricas_transformaciones$Brightness_Normalization <- list(
    R2 = mean(r2_bright, na.rm = TRUE),
    RMSE = mean(rmse_bright, na.rm = TRUE),
    Bias = mean(bias_bright, na.rm = TRUE)
  )
  
  ## ---- 4. Standard Normalization ----
  r2_std <- rmse_std <- bias_std <- c()
  for (i in 1:k) {
    train <- standard_normalization(datos_filtrados[-folds[[i]], ], variables_predictoras)
    test <- standard_normalization(datos_filtrados[folds[[i]], ], variables_predictoras)
    
    modelo_std <- tryCatch(glm(formula, data = train, family = family), error = function(e) NULL)
    
    if (!is.null(modelo_std)) {
      pred_std <- predict(modelo_std, newdata = test, type = "response")
      obs_std <- test[[metal]]
      metricas <- calcular_metricas(obs_std, pred_std)
      r2_std <- c(r2_std, metricas$R2)
      rmse_std <- c(rmse_std, metricas$RMSE)
      bias_std <- c(bias_std, metricas$Bias)
    }
  }
  metricas_transformaciones$Standard_Normalization <- list(
    R2 = mean(r2_std, na.rm = TRUE),
    RMSE = mean(rmse_std, na.rm = TRUE),
    Bias = mean(bias_std, na.rm = TRUE)
  )
  
  ## ---- 5. Derivadas Espectrales ----
  r2_deriv <- rmse_deriv <- bias_deriv <- c()
  for (i in 1:k) {
    train <- derivadas_espectrales(datos_filtrados[-folds[[i]], ], variables_predictoras)
    test <- derivadas_espectrales(datos_filtrados[folds[[i]], ], variables_predictoras)
    
    modelo_deriv <- tryCatch(glm(formula, data = train, family = family), error = function(e) NULL)
    
    if (!is.null(modelo_deriv)) {
      pred_deriv <- predict(modelo_deriv, newdata = test, type = "response")
      obs_deriv <- test[[metal]]
      metricas <- calcular_metricas(obs_deriv, pred_deriv)
      r2_deriv <- c(r2_deriv, metricas$R2)
      rmse_deriv <- c(rmse_deriv, metricas$RMSE)
      bias_deriv <- c(bias_deriv, metricas$Bias)
    }
  }
  metricas_transformaciones$Derivadas_Espectrales <- list(
    R2 = mean(r2_deriv, na.rm = TRUE),
    RMSE = mean(rmse_deriv, na.rm = TRUE),
    Bias = mean(bias_deriv, na.rm = TRUE)
  )
  
  # Guardar resultados
  resultados[[metal]] <- metricas_transformaciones
}

# Imprimir resultados
for (metal in names(resultados)) {
  cat("\nResultados para:", metal, "\n")
  print(resultados[[metal]])
}


##################### Cu, Zn, As

#########################
set.seed(123)

# Configuración de variables y combinación definida
variables_por_respuesta <- list(
  "Cu(g/t)" = c("nir08", "swir16", "swir22", "evi", "ndwi"),
  "Zn(g/t)" = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi", "ndwi"),
  "As(g/t)" = c("nir08", "swir16", "ndvi", "dvi", "rvi", "sbi")
)

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- cor(pred, obs, use = "complete.obs")^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Función de Brightness Normalization
brightness_normalization <- function(df, variables) {
  brightness <- sqrt(rowSums(df[, variables]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, variables] <- df[, variables] / brightness
  return(df)
}

# Función de Standard Normalization
standard_normalization <- function(df, variables) {
  df[, variables] <- scale(df[, variables])
  return(df)
}

# Función de Derivadas Espectrales
derivadas_espectrales <- function(df, variables) {
  for (i in 2:length(variables)) {
    df[[paste0("deriv_", variables[i])]] <- df[[variables[i]]] - df[[variables[i - 1]]]
  }
  return(df)
}

# Lista para almacenar resultados
resultados_globales <- list()

# Iterar por cada metal con una combinación definida
for (variable_respuesta in names(variables_por_respuesta)) {
  cat("\nAnálisis para:", variable_respuesta, "\n")
  
  variables <- variables_por_respuesta[[variable_respuesta]]
  formula <- as.formula(paste0("log1p(`", variable_respuesta, "`) ~ ", paste(variables, collapse = " + ")))
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(variables, variable_respuesta)])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variable_respuesta]], k = k, list = TRUE)
  
  # Almacenar métricas por transformación a las variables predictoras
  metricas_transformaciones <- list(
    Logarítmica = list(R2 = c(), RMSE = c(), Bias = c()),
    Cuadrática = list(R2 = c(), RMSE = c(), Bias = c()),
    Brightness_Normalization = list(R2 = c(), RMSE = c(), Bias = c()),
    Standard_Normalization = list(R2 = c(), RMSE = c(), Bias = c()),
    Derivadas_Espectrales = list(R2 = c(), RMSE = c(), Bias = c())
  )
  
  for (i in 1:k) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    ## ---- Transformación Logarítmica ----
    train_log <- TRAIN
    valid_log <- VALIDATION
    train_log[, variables] <- log1p(train_log[, variables])
    valid_log[, variables] <- log1p(valid_log[, variables])
    
    modelo_glm_log <- tryCatch(glm.nb(formula, data = train_log), error = function(e) NULL)
    
    if (!is.null(modelo_glm_log)) {
      pred_log <- expm1(predict(modelo_glm_log, newdata = valid_log, type = "response"))
      obs_log <- VALIDATION[[variable_respuesta]]
      metricas <- calcular_metricas(obs_log, pred_log)
      metricas_transformaciones$Logarítmica$R2 <- c(metricas_transformaciones$Logarítmica$R2, metricas$R2)
      metricas_transformaciones$Logarítmica$RMSE <- c(metricas_transformaciones$Logarítmica$RMSE, metricas$RMSE)
      metricas_transformaciones$Logarítmica$Bias <- c(metricas_transformaciones$Logarítmica$Bias, metricas$Bias)
    }
    
    ## ---- Brightness Normalization ----
    train_bright <- brightness_normalization(TRAIN, variables)
    valid_bright <- brightness_normalization(VALIDATION, variables)
    
    modelo_bright <- tryCatch(glm.nb(formula, data = train_bright), error = function(e) NULL)
    
    if (!is.null(modelo_bright)) {
      pred_bright <- expm1(predict(modelo_bright, newdata = valid_bright, type = "response"))
      obs_bright <- VALIDATION[[variable_respuesta]]
      metricas <- calcular_metricas(obs_bright, pred_bright)
      metricas_transformaciones$Brightness_Normalization$R2 <- c(metricas_transformaciones$Brightness_Normalization$R2, metricas$R2)
      metricas_transformaciones$Brightness_Normalization$RMSE <- c(metricas_transformaciones$Brightness_Normalization$RMSE, metricas$RMSE)
      metricas_transformaciones$Brightness_Normalization$Bias <- c(metricas_transformaciones$Brightness_Normalization$Bias, metricas$Bias)
    }
    
    ## ---- Standard Normalization ----
    train_std <- standard_normalization(TRAIN, variables)
    valid_std <- standard_normalization(VALIDATION, variables)
    
    modelo_std <- tryCatch(glm.nb(formula, data = train_std), error = function(e) NULL)
    
    if (!is.null(modelo_std)) {
      pred_std <- expm1(predict(modelo_std, newdata = valid_std, type = "response"))
      obs_std <- VALIDATION[[variable_respuesta]]
      metricas <- calcular_metricas(obs_std, pred_std)
      metricas_transformaciones$Standard_Normalization$R2 <- c(metricas_transformaciones$Standard_Normalization$R2, metricas$R2)
      metricas_transformaciones$Standard_Normalization$RMSE <- c(metricas_transformaciones$Standard_Normalization$RMSE, metricas$RMSE)
      metricas_transformaciones$Standard_Normalization$Bias <- c(metricas_transformaciones$Standard_Normalization$Bias, metricas$Bias)
    }
    
    ## ---- Derivadas Espectrales ----
    train_deriv <- derivadas_espectrales(TRAIN, variables)
    valid_deriv <- derivadas_espectrales(VALIDATION, variables)
    
    modelo_deriv <- tryCatch(glm.nb(formula, data = train_deriv), error = function(e) NULL)
    
    if (!is.null(modelo_deriv)) {
      pred_deriv <- expm1(predict(modelo_deriv, newdata = valid_deriv, type = "response"))
      obs_deriv <- VALIDATION[[variable_respuesta]]
      metricas <- calcular_metricas(obs_deriv, pred_deriv)
      metricas_transformaciones$Derivadas_Espectrales$R2 <- c(metricas_transformaciones$Derivadas_Espectrales$R2, metricas$R2)
      metricas_transformaciones$Derivadas_Espectrales$RMSE <- c(metricas_transformaciones$Derivadas_Espectrales$RMSE, metricas$RMSE)
      metricas_transformaciones$Derivadas_Espectrales$Bias <- c(metricas_transformaciones$Derivadas_Espectrales$Bias, metricas$Bias)
    }
  }
  
  # Guardar resultados globales
  resultados_globales[[variable_respuesta]] <- lapply(metricas_transformaciones, function(metrica) {
    list(
      R2 = mean(metrica$R2, na.rm = TRUE),
      RMSE = mean(metrica$RMSE, na.rm = TRUE),
      Bias = mean(metrica$Bias, na.rm = TRUE)
    )
  })
}

# Imprimir resultados para cada variable de respuesta
for (variable_respuesta in names(resultados_globales)) {
  cat("\nResultados para:", variable_respuesta, "\n")
  print(resultados_globales[[variable_respuesta]])
}



################## cuadrática

set.seed(123)

# Configuración de variables y combinación definida
variables_por_respuesta <- list(
  "Cu(g/t)" = c("nir08", "swir16", "swir22", "evi", "ndwi"),
  "Zn(g/t)" = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi", "ndwi"),
  "As(g/t)" = c("nir08", "swir16", "ndvi", "dvi", "rvi", "sbi")
)

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- cor(pred, obs, use = "complete.obs")^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Lista para almacenar resultados
resultados_globales <- list()

# Iterar por cada metal con una combinación definida
for (variable_respuesta in names(variables_por_respuesta)) {
  cat("\nAnálisis para:", variable_respuesta, "\n")
  
  variables <- variables_por_respuesta[[variable_respuesta]]
  formula <- as.formula(paste0("log1p(`", variable_respuesta, "`) ~ ", paste(variables, collapse = " + ")))
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(variables, variable_respuesta)])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variable_respuesta]], k = k, list = TRUE)
  
  # Almacenar métricas por transformación a las variables predictoras
  metricas_transformaciones <- list(
    Logarítmica = list(R2 = c(), RMSE = c(), Bias = c()),
    Cuadrática = list(R2 = c(), RMSE = c(), Bias = c())
  )
  
  for (i in 1:k) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    ## ---- Transformación Logarítmica (variables predictoras) ----
    train_log <- TRAIN
    valid_log <- VALIDATION
    train_log[, variables] <- log1p(train_log[, variables])
    valid_log[, variables] <- log1p(valid_log[, variables])
    
    modelo_glm_log <- tryCatch(
      glm.nb(formula, data = train_log),
      error = function(e) {
        message(sprintf("Error en el ajuste GLM para la transformación logarítmica en el fold %d: %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(modelo_glm_log)) {
      pred_log <- tryCatch(
        expm1(predict(modelo_glm_log, newdata = valid_log, type = "response")),
        error = function(e) {
          message("Error en la predicción logarítmica: ", e$message)
          return(rep(NA, nrow(VALIDATION)))
        }
      )
      obs_log <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas si las predicciones son válidas
      if (!all(is.na(pred_log))) {
        metricas <- calcular_metricas(obs_log, pred_log)
        metricas_transformaciones$Logarítmica$R2 <- c(metricas_transformaciones$Logarítmica$R2, metricas$R2)
        metricas_transformaciones$Logarítmica$RMSE <- c(metricas_transformaciones$Logarítmica$RMSE, metricas$RMSE)
        metricas_transformaciones$Logarítmica$Bias <- c(metricas_transformaciones$Logarítmica$Bias, metricas$Bias)
      }
    }
    
    ## ---- Transformación Cuadrática (variables predictoras) ----
    train_quad <- TRAIN
    valid_quad <- VALIDATION
    train_quad[, variables] <- sqrt(train_quad[, variables])
    valid_quad[, variables] <- sqrt(valid_quad[, variables])
    
    modelo_glm_quad <- tryCatch(
      glm.nb(formula, data = train_quad),
      error = function(e) {
        message(sprintf("Error en el ajuste GLM para la transformación cuadrática en el fold %d: %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(modelo_glm_quad)) {
      pred_quad <- tryCatch(
        expm1(predict(modelo_glm_quad, newdata = valid_quad, type = "response")),
        error = function(e) {
          message("Error en la predicción cuadrática: ", e$message)
          return(rep(NA, nrow(VALIDATION)))
        }
      )
      obs_quad <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas si las predicciones son válidas
      if (!all(is.na(pred_quad))) {
        metricas <- calcular_metricas(obs_quad, pred_quad)
        metricas_transformaciones$Cuadrática$R2 <- c(metricas_transformaciones$Cuadrática$R2, metricas$R2)
        metricas_transformaciones$Cuadrática$RMSE <- c(metricas_transformaciones$Cuadrática$RMSE, metricas$RMSE)
        metricas_transformaciones$Cuadrática$Bias <- c(metricas_transformaciones$Cuadrática$Bias, metricas$Bias)
      }
    }
  }
  
  # Guardar resultados globales
  resultados_globales[[variable_respuesta]] <- lapply(metricas_transformaciones, function(metrica) {
    list(
      R2 = mean(metrica$R2, na.rm = TRUE),
      RMSE = mean(metrica$RMSE, na.rm = TRUE),
      Bias = mean(metrica$Bias, na.rm = TRUE)
    )
  })
}

# Imprimir resultados para cada variable de respuesta
for (variable_respuesta in names(resultados_globales)) {
  cat("\nResultados para:", variable_respuesta, "\n")
  print(resultados_globales[[variable_respuesta]])
}

############# Standard Normalization

set.seed(123)

# Configuración de variables y combinación definida
variables_por_respuesta <- list(
  "Cu(g/t)" = c("nir08", "swir16", "swir22", "evi", "ndwi"),
  "Zn(g/t)" = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi", "ndwi"),
  "As(g/t)" = c("nir08", "swir16", "ndvi", "dvi", "rvi", "sbi")
)

# Número de folds para k-fold cross-validation
k <- 10

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- cor(pred, obs, use = "complete.obs")^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Función de Standard Normalization
standard_normalization <- function(df, variables) {
  df[, variables] <- scale(df[, variables])
  return(df)
}

# Lista para almacenar resultados
resultados_globales <- list()

# Iterar por cada metal con una combinación definida
for (variable_respuesta in names(variables_por_respuesta)) {
  cat("\nAnálisis para:", variable_respuesta, "\n")
  
  variables <- variables_por_respuesta[[variable_respuesta]]
  formula <- as.formula(paste0("log1p(`", variable_respuesta, "`) ~ ", paste(variables, collapse = " + ")))
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(variables, variable_respuesta)])
  
  # Aplicar normalización estándar a las variables predictoras
  datos_filtrados <- standard_normalization(datos_filtrados, variables)
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variable_respuesta]], k = k, list = TRUE)
  
  # Almacenar métricas
  metricas_standard <- list(R2 = c(), RMSE = c(), Bias = c())
  
  for (i in 1:k) {
    # Dividir datos en entrenamiento y validación
    TRAIN <- datos_filtrados[-folds[[i]], ]
    VALIDATION <- datos_filtrados[folds[[i]], ]
    
    modelo_glm_standard <- tryCatch(
      glm.nb(formula, data = TRAIN),
      error = function(e) {
        message(sprintf("Error en el ajuste GLM para normalización estándar en el fold %d: %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(modelo_glm_standard)) {
      pred_standard <- tryCatch(
        expm1(predict(modelo_glm_standard, newdata = VALIDATION, type = "response")),
        error = function(e) {
          message("Error en la predicción de normalización estándar: ", e$message)
          return(rep(NA, nrow(VALIDATION)))
        }
      )
      obs_standard <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas si las predicciones son válidas
      if (!all(is.na(pred_standard))) {
        metricas <- calcular_metricas(obs_standard, pred_standard)
        metricas_standard$R2 <- c(metricas_standard$R2, metricas$R2)
        metricas_standard$RMSE <- c(metricas_standard$RMSE, metricas$RMSE)
        metricas_standard$Bias <- c(metricas_standard$Bias, metricas$Bias)
      }
    }
  }
  
  # Guardar resultados globales
  resultados_globales[[variable_respuesta]] <- list(
    Standard_Normalization = list(
      R2 = mean(metricas_standard$R2, na.rm = TRUE),
      RMSE = mean(metricas_standard$RMSE, na.rm = TRUE),
      Bias = mean(metricas_standard$Bias, na.rm = TRUE)
    )
  )
}

# Imprimir resultados para cada variable de respuesta
for (variable_respuesta in names(resultados_globales)) {
  cat("\nResultados para:", variable_respuesta, "\n")
  print(resultados_globales[[variable_respuesta]])
}


##################################### Grafico de residuos estandarizados

metales_cu_zn_as <- list(
  "Zn(g/t)" = list(
    formula = as.formula("log1p(`Zn(g/t)`) ~ blue + green + red + nir08 + coastal + swir16 + swir22 + rvi + sbi + ndwi")
  )
)

# Iterar por cada metal
for (metal in names(metales_cu_zn_as)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Crear una copia del conjunto de datos para evitar modificar el original
  datos_transformados <- datos3
  
  # Aplicar transformación logarítmica a las variables predictoras
  variables_predictoras <- all.vars(metales_cu_zn_as[[metal]]$formula)[-1]  # Excluir la respuesta
  datos_transformados[, variables_predictoras] <- log1p(datos_transformados[, variables_predictoras])
  
  # Ajustar el modelo GLM
  formula <- metales_cu_zn_as[[metal]]$formula
  modelo <- glm.nb(formula, data = datos_transformados)
  
  # Predicciones en escala original
  predicciones_log <- predict(modelo, type = "response")  # Predicciones en escala logarítmica
  predicciones <- expm1(predicciones_log)  # Transformar a escala original
  residuos <- datos3[[metal]] - predicciones
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Gráfico de residuos estandarizados
  plot(predicciones, residuos_estandarizados,
       main = "Zinc",
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       cex.lab = 1.5,
       cex.axis = 1.5,  
       cex.main = 2.5)
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(predicciones, residuos_estandarizados), col = "blue", lwd = 2)
}


######## Para Co con brightness

# Función para aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE)) 
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, vars] <- df[, vars] / brightness 
  return(df)
}

# Configuración de las fórmulas y familias para Co
metales_ni_co <- list(
  "Co(g/t)" = list(  
    formula = as.formula("`Co(g/t)` ~ nir08 + swir16 + blue + green + coastal + swir22"),
    family = Gamma(link = "log")
  )
)

# Iterar por cada metal
for (metal in names(metales_ni_co)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Filtrar datos sin NA
  datos_filtrados <- na.omit(datos3[, c(all.vars(metales_ni_co[[metal]]$formula), metal)])
  
  # Aplicar Brightness Normalization a las variables predictoras
  variables_predictoras <- all.vars(metales_ni_co[[metal]]$formula)[-1]  # Excluir la respuesta
  datos_filtrados <- brightness_normalization(datos_filtrados, variables_predictoras)
  
  # Ajustar el modelo GLM
  formula <- metales_ni_co[[metal]]$formula
  family <- metales_ni_co[[metal]]$family
  modelo <- glm(formula, data = datos_filtrados, family = family)
  
  # Predicciones y residuos
  predicciones <- predict(modelo, type = "response")
  residuos <- datos_filtrados[[metal]] - predicciones
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Gráfico de residuos estandarizados
  plot(predicciones, residuos_estandarizados,
       main = "Cobalto",
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       cex.lab = 1.5,
       cex.axis = 1.5,    # Tamaño de los números en los ejes
       cex.main = 2.5)
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(predicciones, residuos_estandarizados), col = "blue", lwd = 2)
}


######################################################################### 5) Mejor Modelo con GLM, 
########   Gráfico de Valores Predichos vs Reales, el Análisis de Sensibilidad y el Análisis de Contribución de Variables

set.seed(123)

# Configuración de variables, fórmulas y familias
configuracion_ni_co <- list(
  "Ni(g/t)" = list(
    formula = as.formula("`Ni(g/t)` ~ green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22"),
    family = gaussian(link = "identity"),
    transform_predictors = FALSE
  ),
  "Co(g/t)" = list(
    formula = as.formula("`Co(g/t)` ~ nir08 + swir16 + blue + green + coastal + swir22"),
    family = Gamma(link = "log"),
    transform_predictors = TRUE
  )
)

# Lista para almacenar resultados
resultados_ni_co <- list()

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

# Iterar por cada metal
for (metal in names(configuracion_ni_co)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  config <- configuracion_ni_co[[metal]]
  formula <- config$formula
  family <- config$family
  transform_predictors <- config$transform_predictors
  
  # Filtrar datos sin NA
  variables <- all.vars(formula)
  datos_filtrados <- na.omit(datos3[, variables])
  
  # Aplicar transformación logarítmica a las variables predictoras si corresponde
  if (transform_predictors) {
    for (var in variables[-1]) {  # Excluir la variable de respuesta
      datos_filtrados[[var]] <- log1p(datos_filtrados[[var]])
    }
  }
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variables[1]]], k = k, list = TRUE)
  
  # Lista para almacenar R2
  r2_folds <- c()
  
  for (i in 1:k) {
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
    pred <- tryCatch({
      predict(modelo, newdata = VALIDATION, type = "response")
    }, error = function(e) {
      message("Error en la predicción: ", e$message)
      return(rep(NA, nrow(VALIDATION)))
    })
    
    # Calcular R2
    r2 <- calcular_r2(VALIDATION[[variables[1]]], pred)
    r2_folds <- c(r2_folds, r2)
  }
  
  # Promedio de R2
  r2_promedio <- mean(r2_folds, na.rm = TRUE)
  
  # Guardar resultados
  resultados_ni_co[[metal]] <- list(modelo = modelo, R2 = r2_promedio)
  
  # Imprimir R2
  cat("R2 promedio para", metal, ":", r2_promedio, "\n")
}

# Imprimir resultados globales
for (metal in names(resultados_ni_co)) {
  cat("\nResultados del modelo para:", metal, "\n")
  print(resultados_ni_co[[metal]])
}


########################## Cu, Zn y As
set.seed(123)

# Configuración de variables y fórmulas
configuracion_cu_zn_as <- list(
  "Cu(g/t)" = list(
    formula = as.formula("log1p(`Cu(g/t)`) ~ nir08 + swir16 + swir22 + evi + ndwi")
  ),
  "Zn(g/t)" = list(
    formula = as.formula("`Zn(g/t)` ~ sqrt(blue) + sqrt(green) + sqrt(red) + sqrt(nir08) + sqrt(coastal) + sqrt(swir16) + sqrt(swir22) + sqrt(rvi) + sqrt(sbi) + sqrt(ndwi)")
  ),
  "As(g/t)" = list(
    formula = as.formula("log1p(`As(g/t)`) ~ nir08 + swir16 + ndvi + dvi + rvi + sbi")
  )
)

# Lista para almacenar resultados
resultados_cu_zn_as <- list()

# Número de folds para k-fold cross-validation
k <- 10

# Iterar por cada metal
for (metal in names(configuracion_cu_zn_as)) {
  cat("\nAnálisis para:", metal, "\n")
  
  # Extraer configuración
  formula <- configuracion_cu_zn_as[[metal]]$formula
  
  # Filtrar datos sin NA
  variables <- all.vars(formula)
  datos_filtrados <- na.omit(datos3[, variables])
  
  # Crear folds para k-fold cross-validation
  folds <- createFolds(datos_filtrados[[variables[1]]], k = k, list = TRUE)
  
  # Lista para almacenar R2
  r2_folds <- c()
  
  for (i in 1:k) {
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
    
    # Calcular R2
    r2 <- calcular_r2(VALIDATION[[variables[1]]], pred)
    r2_folds <- c(r2_folds, r2)
  }
  
  # Promedio de R2
  r2_promedio <- mean(r2_folds, na.rm = TRUE)
  
  # Guardar resultados
  resultados_cu_zn_as[[metal]] <- list(modelo = modelo, R2 = r2_promedio)
  
  # Imprimir R2
  cat("R2 promedio para", metal, ":", r2_promedio, "\n")
}

# Imprimir resultados globales
for (metal in names(resultados_cu_zn_as)) {
  cat("\nResultados del modelo para:", metal, "\n")
  print(resultados_cu_zn_as[[metal]])
}

################################ Gráfico de valores predichos vs reales

# Gráfico con ajuste proporcional de ambos ejes
graficar_pred_vs_real <- function(modelo, datos, metal, transformacion_respuesta = NULL) {
  # Predicciones con el modelo
  pred <- tryCatch({
    predict(modelo, newdata = datos, type = "response")
  }, error = function(e) {
    message(sprintf("Error en la predicción para %s: %s", metal, e$message))
    return(rep(NA, nrow(datos)))
  })
  
  # Deshacer transformación en la variable de respuesta si corresponde
  if (!is.null(transformacion_respuesta)) {
    if (transformacion_respuesta == "log1p") {
      pred <- expm1(pred)
    }
  }
  
  # Crear DataFrame con valores reales y predichos
  data_plot <- data.frame(
    Real = datos[[metal]],
    Predicho = pred
  )
  
  # Calcular el límite máximo
  max_val <- max(max(data_plot$Real, na.rm = TRUE), max(data_plot$Predicho, na.rm = TRUE)) * 1.1
  
  # Graficar
  ggplot(data_plot, aes(x = Real, y = Predicho)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = paste("Valores Predichos vs Reales para", metal),
      x = "Valores Observados",
      y = "Valores Predichos"
    ) +
    coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    theme_minimal()
}

# Iterar sobre los modelos y generar gráficos
# Modelos para Ni y Co
for (metal in names(resultados_ni_co)) {
  cat("\nGenerando gráfico para:", metal, "\n")
  
  datos_modelo <- na.omit(datos3[, all.vars(configuracion_ni_co[[metal]]$formula)])  # Filtrar datos
  if (configuracion_ni_co[[metal]]$transform_predictors) {
    # Aplicar transformación logarítmica a las variables predictoras si corresponde
    for (var in all.vars(configuracion_ni_co[[metal]]$formula)[-1]) {
      datos_modelo[[var]] <- log1p(datos_modelo[[var]])
    }
  }
  
  print(graficar_pred_vs_real(
    modelo = resultados_ni_co[[metal]]$modelo,
    datos = datos_modelo,
    metal = metal
  ))
}

# Modelos para Cu, Zn y As
for (metal in names(resultados_cu_zn_as)) {
  cat("\nGenerando gráfico para:", metal, "\n")
  
  datos_modelo <- na.omit(datos3[, all.vars(configuracion_cu_zn_as[[metal]]$formula)])  # Filtrar datos
  
  if (metal == "Zn(g/t)") {
    # Aplicar transformación cuadrática a las variables predictoras
    for (var in all.vars(configuracion_cu_zn_as[[metal]]$formula)[-1]) {
      datos_modelo[[var]] <- sqrt(datos_modelo[[var]])
    }
  }
  
  transformacion <- if (grepl("log1p", as.character(configuracion_cu_zn_as[[metal]]$formula)[2])) "log1p" else NULL
  
  print(graficar_pred_vs_real(
    modelo = resultados_cu_zn_as[[metal]]$modelo,
    datos = datos_modelo,
    metal = metal,
    transformacion_respuesta = transformacion
  ))
}




################################ Análisis de Sensibilidad

# Análisis de Sensibilidad basado en los coeficientes del modelo
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
    
    # Imprimir tabla
    cat("\nAnálisis de Sensibilidad para:", metal, "\n")
    print(coeficientes_df)
  }
}

# Análisis para Ni y Co
for (metal in names(resultados_ni_co)) {
  cat("\n### Análisis de Sensibilidad para", metal, "###\n")
  analizar_sensibilidad(resultados_ni_co[[metal]]$modelo, metal)
}

# Análisis para Cu, Zn y As
for (metal in names(resultados_cu_zn_as)) {
  cat("\n### Análisis de Sensibilidad para", metal, "###\n")
  analizar_sensibilidad(resultados_cu_zn_as[[metal]]$modelo, metal)
}




