######################## Librerías generales
library(caret)
library(MASS)
library(readxl)
library(dplyr)

############################################################### 1) Preparación de datos

datos2 <- read_excel("C:/Users/PC//Desktop/TESIS/reflectancia_todos_poligonos.xlsx")
datos3 <- datos2[, !(names(datos2) %in% c('Name', 'Empresa o DueÃ±o', 'Faena','Fecha'))]
datos3 <- datos3 %>%  mutate(across(all_of(setdiff(names(datos3), "Estado")), as.numeric))

############################################################### 1.0) Correlación entre variables dependientes e independientes

library(ggplot2)
library(reshape2)

variables_dependientes <- c("Ni(g/t)", "Cu(g/t)", "Zn(g/t)", "Co(g/t)", "As(g/t)")
variables_independientes <- c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22",
                              "ndvi", "dvi", "rvi", "sbi", "evi", "ndwi")

# Crear una tabla vacía para almacenar las correlaciones
tabla_correlaciones <- data.frame(matrix(nrow = length(variables_dependientes), ncol = length(variables_independientes)))
rownames(tabla_correlaciones) <- variables_dependientes
colnames(tabla_correlaciones) <- variables_independientes

# Calcular correlaciones de Pearson
for (dep in variables_dependientes) {
  for (indep in variables_independientes) {
    correlacion <- cor(datos3[[dep]], datos3[[indep]], use = "complete.obs", method = "pearson")
    tabla_correlaciones[dep, indep] <- correlacion
  }
}

# Convertir los nombres de las filas en una columna explícita
tabla_correlaciones$Dependiente <- rownames(tabla_correlaciones)

# Transformar la tabla de correlaciones a formato largo
tabla_correlaciones_larga <- melt(tabla_correlaciones, id.vars = "Dependiente", variable.name = "Independiente", value.name = "Correlacion")

ggplot(tabla_correlaciones_larga, aes(x = Independiente, y = Dependiente, fill = Correlacion)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab") +
  geom_text(aes(label = round(Correlacion, 2)), color = "black", size = 3) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "Variables independientes", y = "Variables dependientes")


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

############################################################### 2) Experimento 1 

############################################################### 2.1) Selección de variables
library(pls)

# Función para generar combinaciones dinámicas
generar_combinaciones <- function(variables_adicionales) {
  combinaciones <- list()
  for (k in 1:length(variables_adicionales)) {
    combinaciones <- c(combinaciones, combn(variables_adicionales, k, simplify = FALSE))
  }
  return(combinaciones)
}

# Configuración de variables base y adicionales por variable de respuesta
variables_base_por_respuesta <- list(
  `Ni(g/t)` = c("green", "nir08"),
  `Cu(g/t)` = c("nir08", "swir16", "swir22"),
  `Zn(g/t)` = c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22"),
  `Co(g/t)` = c("nir08", "swir16"),
  `As(g/t)` = c("nir08")
)

variables_adicionales_por_respuesta <- list(
  `Ni(g/t)` = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "red", "coastal", "swir16", "swir22"),
  `Cu(g/t)` = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "green", "red", "coastal"),
  `Zn(g/t)` = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi"),
  `Co(g/t)` = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "green", "red", "coastal", "swir22"),
  `As(g/t)` = c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", "green", "red", "coastal", "swir16", "swir22")
)
conteo_combinaciones <- list()

# Generar combinaciones para cada metal y contar
for (variable_respuesta in names(variables_base_por_respuesta)) {
  # Definir variables base y adicionales específicas
  variables_base <- variables_base_por_respuesta[[variable_respuesta]]
  variables_adicionales <- variables_adicionales_por_respuesta[[variable_respuesta]]
  
  # Generar combinaciones
  combinaciones_adicionales <- generar_combinaciones(variables_adicionales)
  combinaciones <- lapply(combinaciones_adicionales, function(x) c(variables_base, x))
  
  # Contar combinaciones
  conteo_combinaciones[[variable_respuesta]] <- length(combinaciones)
}

# Mostrar el conteo de combinaciones por metal
cat("Conteo de combinaciones por metal:\n")
print(conteo_combinaciones)

# Combinación fija a incluir
combinacion_fija <- c("blue", "green", "red", "nir08", "coastal", "swir16", "swir22")

# DataFrame para almacenar todos los resultados
resultados_totales <- data.frame(VariableRespuesta = character(), 
                                 Combinacion = character(), 
                                 R2 = numeric(), 
                                 RMSE = numeric(), 
                                 BIAS = numeric(), 
                                 stringsAsFactors = FALSE)

# Iterar por cada variable de respuesta
for (variable_respuesta in names(variables_base_por_respuesta)) {
  # Definir variables base y adicionales específicas
  variables_base <- variables_base_por_respuesta[[variable_respuesta]]
  variables_adicionales <- variables_adicionales_por_respuesta[[variable_respuesta]]
  
  # Generar combinaciones
  combinaciones_adicionales <- generar_combinaciones(variables_adicionales)
  combinaciones <- lapply(combinaciones_adicionales, function(x) c(variables_base, x))
  
  # Agregar la combinación fija al conjunto de combinaciones
  combinaciones <- c(combinaciones, list(combinacion_fija))
  
  # Configuración para validación cruzada específica de la variable
  set.seed(123)
  folds <- createFolds(datos3[[variable_respuesta]], k = 10)
  
  # DataFrame para resultados de esta variable de respuesta
  resultados_variable <- data.frame(Combinacion = character(), R2 = numeric(), 
                                    RMSE = numeric(), BIAS = numeric(), 
                                    stringsAsFactors = FALSE)
  
  for (variables in combinaciones) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:10) { # 10 particiones
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Crear fórmula dinámica para PLS
      formula <- as.formula(paste0("`", variable_respuesta, "` ~ ", paste(variables, collapse = " + ")))
      
      # Ajustar el modelo PLS con validación cruzada
      modelo_pls <- plsr(formula, data = TRAIN, validation = "CV")
      ncomp_opt <- selectNcomp(modelo_pls, method = "onesigma", plot = FALSE)
      
      # Asegurar un número válido de componentes
      ncomp_opt <- max(1, min(ncomp_opt, dim(TRAIN)[2] - 1))
      
      # Predicción en conjunto de validación
      pred <- predict(modelo_pls, VALIDATION, ncomp = ncomp_opt)
      
      # Verificar formato de predicciones y adaptar
      pred <- if (is.vector(pred)) {
        as.numeric(pred)
      } else if (is.matrix(pred)) {
        pred[, 1]
      } else if (is.array(pred) && length(dim(pred)) > 2) {
        pred[, 1, 1]
      } else {
        stop("Error: formato de predicción desconocido")
      }
      
      obs <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas
      r2[i] <- cor(pred, obs, use = "complete.obs")^2
      rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
      bias[i] <- 1 - coef(lm(pred ~ obs - 1))
    }
    
    # Promedio de métricas para la combinación actual
    resultados_variable <- rbind(resultados_variable, data.frame(
      Combinacion = paste(variables, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
  
  # Ordenar resultados y seleccionar top 3 combinaciones dinámicas
  resultados_top1 <- resultados_variable[order(-resultados_variable$R2), ][1:1, ]
  
  # Agregar la combinación fija
  resultado_fijo <- resultados_variable[resultados_variable$Combinacion == paste(combinacion_fija, collapse = ", "), ]
  
  # Consolidar resultados finales para esta variable
  resultados_top1 <- rbind(resultados_top1, resultado_fijo)
  resultados_top1$VariableRespuesta <- variable_respuesta
  
  # Almacenar en el DataFrame total
  resultados_totales <- rbind(resultados_totales, resultados_top1)
}

print(resultados_totales)


############################################################### 2.2) Gráficos de residuos y Q-Q

# Para cada uno de los metales se aplican las combinaciones de variables de a continuación para graficar:
# Ni : green + nir08 + swir16+ ndvi + ndwi
# Cu: green + blue + red + nir08 + coastal + swir16 + swir22 + ndvi + rvi + sbi + ndwi
# Zn: green + blue + red + nir08 + coastal + swir16 + swir22 + sbi + ndwi
# Co: blue + nir08 + coastal + swir16 + ndvi + dvi + evi + ndwi
# As: red + nir08 + coastal + swir16 + swir22 + rvi + sbi

formula <- as.formula("`Zn(g/t)` ~ red + nir08 + coastal + swir16 + swir22 + rvi + sbi")

modelo_pls <- plsr(formula, data = datos3, validation = "CV", segments = 10,  segment.type = "random")

ncomp_opt <- selectNcomp(modelo_pls, method = "onesigma", plot = FALSE)
print(ncomp_opt)
# Predicciones con el número óptimo de componentes
predicciones <- predict(modelo_pls, datos3, ncomp = max(1, ncomp_opt))
residuos <- datos3$`As(g/t)` - predicciones
residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))

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


# Gráfico Q-Q
qqnorm(residuos,
       main = "Arsénico",
       xlab = "Cuantiles Teóricos",
       ylab = "Cuantiles Observados",
       cex.main = 2.5,
       cex.lab = 1.5,
       cex.axis = 1.5)
qqline(residuos, col = "blue")



############################################################### 3) Experimento 2

############################################################### 3.1) PLS con k-fold y Transformaciones en Y


set.seed(123)

k <- 10

# Función para aplicar las transformaciones a los datos
aplicar_transformacion <- function(data, transform_func) {
  transform_func(data)
}

# Función para generar resultados con PLS y transformaciones
evaluar_pls_transformaciones <- function(datos, variable_respuesta, combinaciones, transformaciones, folds) {
  resultados <- data.frame(Transformacion = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())
  
  for (nombre_transformacion in names(transformaciones)) {
    trans <- transformaciones[[nombre_transformacion]]$transform
    inv_trans <- transformaciones[[nombre_transformacion]]$inverse
    # Crear columna transformada sin comillas invertidas en el nombre
    columna_transformada <- paste0(variable_respuesta, "_trans")
    datos[[columna_transformada]] <- trans(datos[[variable_respuesta]])
    
    for (variables in combinaciones) {
      r2 <- c()
      rmse <- c()
      bias <- c()
      
      for (i in 1:k) {
        TRAIN <- datos[-folds[[i]], ]
        VALIDATION <- datos[folds[[i]], ]
        
        # Crear fórmula dinámica para PLS
        formula <- as.formula(paste0("`", columna_transformada, "` ~ ", paste(variables, collapse = " + ")))
        modelo_pls <- plsr(formula, data = TRAIN, validation = "CV")
        ncomp_opt <- max(1, selectNcomp(modelo_pls, method = "onesigma", plot = FALSE)) # Forzar al menos 1 componente
        
        # Predicción en el espacio transformado
        pred_trans <- predict(modelo_pls, VALIDATION, ncomp = ncomp_opt)
        pred_trans <- as.vector(pred_trans)
        
        # Transformación inversa al espacio original
        pred <- inv_trans(pred_trans)
        obs <- VALIDATION[[variable_respuesta]]
        
        # Calcular métricas en el espacio original
        r2[i] <- cor(pred, obs, use = "complete.obs")^2
        rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        bias[i] <- 1 - coef(lm(pred ~ obs - 1))
      }
      
      # Promedio de métricas
      resultados <- rbind(resultados, data.frame(
        Transformacion = nombre_transformacion,
        Combinacion = paste(variables, collapse = ", "),
        R2 = mean(r2, na.rm = TRUE),
        RMSE = mean(rmse, na.rm = TRUE),
        BIAS = mean(bias, na.rm = TRUE)
      ))
    }
  }
  
  return(resultados)
}

# Combinaciones y transformaciones por variable de respuesta
variables_por_respuesta <- list(
  "Ni(g/t)" = list(
    combinaciones = list(
      c("green", "nir08", "swir16", "ndvi", "ndwi")
    )
  ),
  "Cu(g/t)" = list(
    combinaciones = list(
      c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "ndvi", "rvi", "sbi", "ndwi")
    )
  ),
  "Zn(g/t)" = list(
    combinaciones = list(
      c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi")
    )
  ),
  "Co(g/t)" = list(
    combinaciones = list(
      c("blue", "nir08", "coastal", "swir16", "ndvi", "dvi", "evi", "ndwi")
    )
  ),
  "As(g/t)" = list(
    combinaciones = list(
      c("red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi")
    )
  )
)

# Definir transformaciones
transformaciones <- list(
  log = list(transform = function(x) log(x + 1), inverse = function(x) exp(x) - 1),
  sqrt = list(transform = sqrt, inverse = function(x) x^2),
  none = list(transform = identity, inverse = identity)
)

# Evaluar cada variable de respuesta
resultados_globales <- list()

for (variable in names(variables_por_respuesta)) {
  combinaciones <- variables_por_respuesta[[variable]]$combinaciones
  folds <- createFolds(datos3[[variable]], k = k)
  
  resultados <- evaluar_pls_transformaciones(datos3, variable, combinaciones, transformaciones, folds)
  
  resultados_globales[[variable]] <- resultados
}

# Imprimir resultados
for (variable in names(resultados_globales)) {
  cat("\nResultados para:", variable, "\n")
  print(resultados_globales[[variable]])
}


############################################################### 3.2) Gráficos de residuos y Q-Q

# Para cada uno de los metales se aplican las combinaciones de variables de a continuación para graficar:
# Ni : green + nir08 + swir16+ ndvi + ndwi
# Cu: green + blue + red + nir08 + coastal + swir16 + swir22 + ndvi + rvi + sbi + ndwi
# Zn: green + blue + red + nir08 + coastal + swir16 + swir22 + sbi + ndwi
# Co: blue + nir08 + coastal + swir16 + ndvi + dvi + evi + ndwi
# As: red + nir08 + coastal + swir16 + swir22 + rvi + sbi

# Funciones de transformación logarítmica e inversa
transformacion_log <- function(x) log(x + 1)
inversa_log <- function(x) exp(x) - 1

# Aplicar la transformación a la variable de respuesta
datos3$`Zn_trans` <- transformacion_log(datos3$`Zn(g/t)`)

# Crear fórmula con la variable transformada
formula <- as.formula("Zn_trans ~ green + blue + red + nir08 + coastal + swir16 + swir22 + sbi + ndwi")

# Modelo PLSR con la variable transformada
modelo_pls <- plsr(formula, data = datos3, validation = "CV", segments = 10, segment.type = "random")

# Seleccionar número óptimo de componentes
ncomp_opt <- selectNcomp(modelo_pls, method = "onesigma", plot = FALSE)
print(ncomp_opt)

# Predicciones en el espacio transformado
pred_trans <- predict(modelo_pls, datos3, ncomp = max(1, ncomp_opt))
pred <- inversa_log(pred_trans) # Transformación inversa para el espacio original

# Calcular residuos en el espacio original
residuos <- datos3$`Zn(g/t)` - pred
residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))

# Gráfico de residuos estandarizados
plot(pred, residuos_estandarizados,
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 2)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(pred, residuos_estandarizados), col = "blue", lwd = 2)

# Gráfico Q-Q
qqnorm(residuos,
       cex.main = 2.5,
       cex.lab = 2)
qqline(residuos, col = "blue")





############################################################### 4) Experimento 3

############################################################### 4.1) PLS con k-fold y Transformaciones en X


# Tratar la negatividad de las variables predictoras

variables_predictoras <- c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", 
                           "green", "red", "nir08", "coastal", "swir16", "swir22")

# Sumar 1 a todas las columnas de las variables predictoras
datos3[variables_predictoras] <- datos3[variables_predictoras] + 1


set.seed(123)
k <- 10

# Crear combinaciones dinámicamente para cada variable de respuesta
variables_por_respuesta <- list(
  "Ni(g/t)" = list(combinaciones = list(c("green", "nir08", "swir16", "ndvi", "ndwi"))),
  "Cu(g/t)" = list(combinaciones = list(c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "ndvi", "rvi", "sbi", "ndwi"))),
  "Zn(g/t)" = list(combinaciones = list(c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi"))),
  "Co(g/t)" = list(combinaciones = list(c("blue", "nir08", "coastal", "swir16", "ndvi", "dvi", "evi", "ndwi"))),
  "As(g/t)" = list(combinaciones = list(c("red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi")))
)

# Funciones de transformación
transformaciones <- list(
  log = function(x) log(x),
  sqrt = function(x) sqrt(x),
  brightness = function(df, vars) {
    brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
    brightness[brightness == 0] <- 1
    df[, vars] <- df[, vars] / brightness
    return(df)
  },
  standard = function(df, vars) {
    df[, vars] <- scale(df[, vars])
    return(df)
  }
)

# DataFrame para resultados globales
resultados_globales <- data.frame(VariableRespuesta = character(), Transformacion = character(), 
                                  Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por variables de respuesta
for (variable_respuesta in names(variables_por_respuesta)) {
  combinaciones <- variables_por_respuesta[[variable_respuesta]]$combinaciones
  folds <- createFolds(datos3[[variable_respuesta]], k = k)
  
  # Loop por transformaciones
  for (nombre_trans in names(transformaciones)) {
    transform_func <- transformaciones[[nombre_trans]]
    
    for (variables in combinaciones) {
      r2 <- c()
      rmse <- c()
      bias <- c()
      
      for (i in 1:k) {
        TRAIN <- datos3[-folds[[i]], ]
        VALIDATION <- datos3[folds[[i]], ]
        
        # Aplicar transformación
        if (nombre_trans == "brightness" || nombre_trans == "standard") {
          TRAIN <- transform_func(TRAIN, variables)
          VALIDATION <- transform_func(VALIDATION, variables)
        } else {
          TRAIN[variables] <- transform_func(TRAIN[variables])
          VALIDATION[variables] <- transform_func(VALIDATION[variables])
        }
        
        # Crear fórmula dinámica para PLS
        formula <- as.formula(paste0("`", variable_respuesta, "` ~ ", paste(variables, collapse = " + ")))
        modelo_pls <- plsr(formula, data = TRAIN, validation = "CV")
        ncomp_opt <- max(1, selectNcomp(modelo_pls, method = "onesigma", plot = FALSE))
        
        # Predicción en el conjunto de validación
        pred <- predict(modelo_pls, VALIDATION, ncomp = ncomp_opt)
        pred <- if (is.vector(pred)) {
          as.numeric(pred)
        } else if (is.matrix(pred)) {
          pred[, 1]
        } else if (is.array(pred) && length(dim(pred)) > 2) {
          pred[, 1, 1]
        } else {
          stop("Error: formato de predicción desconocido")
        }
        
        obs <- VALIDATION[[variable_respuesta]]
        
        # Calcular métricas
        r2[i] <- cor(pred, obs, use = "complete.obs")^2
        rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        bias[i] <- 1 - coef(lm(pred ~ obs - 1))
      }
      
      # Guardar resultados
      resultados_globales <- rbind(resultados_globales, data.frame(
        VariableRespuesta = variable_respuesta,
        Transformacion = nombre_trans,
        Combinacion = paste(variables, collapse = ", "),
        R2 = mean(r2, na.rm = TRUE),
        RMSE = mean(rmse, na.rm = TRUE),
        BIAS = mean(bias, na.rm = TRUE)
      ))
    }
  }
}

print(resultados_globales)




################### derivadas espectrales

set.seed(123)
k <- 10

# Crear combinaciones dinámicamente para cada variable de respuesta
variables_por_respuesta <- list(
  "Ni(g/t)" = list(combinaciones = list(c("green", "nir08", "swir16", "ndvi", "ndwi"))),
  "Cu(g/t)" = list(combinaciones = list(c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "ndvi", "rvi", "sbi", "ndwi"))),
  "Zn(g/t)" = list(combinaciones = list(c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi"))),
  "Co(g/t)" = list(combinaciones = list(c("blue", "nir08", "coastal", "swir16", "ndvi", "dvi", "evi", "ndwi"))),
  "As(g/t)" = list(combinaciones = list(c("red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi")))
)

# Función para aplicar Derivadas Espectrales
spectral_derivative <- function(df, variables) {
  # Convertir a matriz para calcular las diferencias entre columnas consecutivas
  deriv_matrix <- t(apply(df[, variables], 1, diff))
  deriv_df <- as.data.frame(deriv_matrix)
  colnames(deriv_df) <- paste0("deriv_", 1:ncol(deriv_df))  # Renombrar columnas derivadas
  return(deriv_df)
}

# DataFrame para resultados globales
resultados_globales <- data.frame(VariableRespuesta = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por variables de respuesta
for (variable_respuesta in names(variables_por_respuesta)) {
  combinaciones <- variables_por_respuesta[[variable_respuesta]]$combinaciones
  folds <- createFolds(datos3[[variable_respuesta]], k = k)
  
  for (variables in combinaciones) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Aplicar Derivadas Espectrales a las variables independientes
      TRAIN_deriv <- spectral_derivative(TRAIN, variables)
      VALIDATION_deriv <- spectral_derivative(VALIDATION, variables)
      
      # Crear fórmula dinámica para PLS
      formula <- as.formula(paste0("`", variable_respuesta, "` ~ ", paste(colnames(TRAIN_deriv), collapse = " + ")))
      
      # Combinar TRAIN y VALIDATION con las derivadas
      TRAIN_combined <- cbind(TRAIN[, variable_respuesta, drop = FALSE], TRAIN_deriv)
      VALIDATION_combined <- cbind(VALIDATION[, variable_respuesta, drop = FALSE], VALIDATION_deriv)
      
      # Ajustar el modelo PLS con validación cruzada para encontrar el número óptimo de componentes
      modelo_pls <- plsr(formula, data = TRAIN_combined, validation = "CV")
      
      # Asegurar que el número de componentes no exceda el número de variables disponibles
      max_components <- min(ncol(TRAIN_combined) - 1, dim(modelo_pls$scores)[2])  # Excluye la respuesta
      ncomp_opt <- max(1, min(selectNcomp(modelo_pls, method = "onesigma", plot = FALSE), max_components))
      
      # Predicción en el conjunto de validación
      pred <- predict(modelo_pls, VALIDATION_combined, ncomp = ncomp_opt)
      pred <- if (is.vector(pred)) {
        as.numeric(pred)
      } else if (is.matrix(pred)) {
        pred[, 1]
      } else if (is.array(pred) && length(dim(pred)) > 2) {
        pred[, 1, 1]
      } else {
        stop("Error: formato de predicción desconocido")
      }
      
      obs <- VALIDATION[[variable_respuesta]]
      
      # Calcular métricas
      r2[i] <- cor(pred, obs, use = "complete.obs")^2
      rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
      bias[i] <- 1 - coef(lm(pred ~ obs - 1))
    }
    
    # Guardar resultados para esta combinación
    resultados_globales <- rbind(resultados_globales, data.frame(
      VariableRespuesta = variable_respuesta,
      Combinacion = paste(variables, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
}

print(resultados_globales)

############################################################### 4.2) Gráficos de residuos y Q-Q

# Función para aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE)) 
  brightness[brightness == 0] <- 1 
  df[, vars] <- df[, vars] / brightness 
  return(df)
}

# Aplicar Brightness Normalization a las variables predictoras
variables_predictoras <- c("blue", "nir08", "coastal", "swir16", "ndvi", "dvi", "evi","ndwi")
datos3_brightness <- brightness_normalization(datos3, variables_predictoras)

# Crear fórmula con las variables normalizadas
formula <- as.formula("`Co(g/t)` ~ blue + nir08 + coastal + swir16 + ndvi + dvi + evi + ndwi")

# Modelo PLSR con Brightness Normalization
modelo_pls <- plsr(formula, data = datos3_brightness, validation = "CV", segments = 10, segment.type = "random")

# Seleccionar número óptimo de componentes
ncomp_opt <- selectNcomp(modelo_pls, method = "onesigma", plot = FALSE)
print(ncomp_opt)

# Predicciones en el espacio original
pred <- predict(modelo_pls, datos3_brightness, ncomp = max(1, ncomp_opt))
pred <- as.numeric(pred)

# Calcular residuos en el espacio original
residuos <- datos3$`Co(g/t)` - pred
residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))

# Gráfico de residuos estandarizados
plot(pred, residuos_estandarizados,
     main = "Cobalto",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2.5)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(pred, residuos_estandarizados), col = "blue", lwd = 2)

# Gráfico Q-Q
qqnorm(residuos,
       main = "Cobalto",
       xlab = "Cuantiles Teóricos",
       ylab = "Cuantiles Observados",
       cex.main = 2.5,
       cex.lab = 1.5,
       cex.axis = 1.5)
qqline(residuos, col = "blue")



############ Para el Zn

# Función para aplicar Standard Normalization
standard_normalization <- function(df, vars) {
  df[, vars] <- scale(df[, vars], center = TRUE, scale = TRUE) # Centrado y escalado
  return(df)
}

# Aplicar Standard Normalization a las variables predictoras
variables_predictoras <- c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi")
datos3_standard <- standard_normalization(datos3, variables_predictoras)

# Crear fórmula con las variables normalizadas
formula <- as.formula("`Zn(g/t)` ~ green + blue + red + nir08 + coastal + swir16 + swir22 + sbi + ndwi")

# Modelo PLSR con Standard Normalization
modelo_pls <- plsr(formula, data = datos3_standard, validation = "CV", segments = 10, segment.type = "random")

# Seleccionar número óptimo de componentes
ncomp_opt <- selectNcomp(modelo_pls, method = "onesigma", plot = FALSE)
print(ncomp_opt)

# Predicciones en el espacio original
pred <- predict(modelo_pls, datos3_standard, ncomp = max(1, ncomp_opt))
pred <- as.numeric(pred)

# Calcular residuos en el espacio original
residuos <- datos3$`Zn(g/t)` - pred
residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))

# Gráfico de residuos estandarizados
plot(pred, residuos_estandarizados,
     main = "Zinc",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2.5)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(pred, residuos_estandarizados), col = "blue", lwd = 2)

# Gráfico Q-Q
qqnorm(residuos,
       main = "Zinc",
       xlab = "Cuantiles Teóricos",
       ylab = "Cuantiles Observados",
       cex.main = 2.5,
       cex.lab = 1.5,
       cex.axis = 1.5)
qqline(residuos, col = "blue")


######################################################################### 5) Mejor Modelo con PLS, 
########   Gráfico de Valores Predichos vs Reales, el Análisis de Sensibilidad y el Análisis de Contribución de Variables

###################################### Modelos: Ni, Cu, As

set.seed(123)
k <- 10

# Configuración de metales y variables predictoras
metales <- list(
  "Ni(g/t)" = c("green", "nir08", "swir16", "ndvi", "ndwi"),
  "Cu(g/t)" = c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "ndvi", "rvi", "sbi", "ndwi"),
  "As(g/t)" = c("red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi")
)

# DataFrame para almacenar resultados de los modelos
resultados_modelos <- list()

# Loop por cada metal
for (metal in names(metales)) {
  cat("\nGenerando modelo para:", metal, "\n")
  
  # Variables predictoras
  vars <- metales[[metal]]
  
  # Crear fórmula dinámica
  formula <- as.formula(paste0("`", metal, "` ~ ", paste(vars, collapse = " + ")))
  
  # Validación cruzada
  folds <- createFolds(datos3[[metal]], k = k)
  
  # Crear modelo PLSR
  modelo_pls <- plsr(formula, data = datos3, validation = "CV", segments = k)
  
  # Selección del número óptimo de componentes
  ncomp_opt <- max(1, selectNcomp(modelo_pls, method = "onesigma", plot = FALSE))
  cat("Número óptimo de componentes para", metal, ":", ncomp_opt, "\n")
  
  # Guardar modelo y componentes óptimos
  resultados_modelos[[metal]] <- list(modelo = modelo_pls, ncomp_opt = ncomp_opt)
}

resultados_modelos


###################################### Modelo Co

set.seed(123)
k <- 10

# Variables predictoras para Co
variables_co <- c("blue", "nir08", "coastal", "swir16", "ndvi", "dvi", "evi", "ndwi")

# Función de Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por 0
  df[, vars] <- df[, vars] / brightness
  return(df)
}

# Aplicar Brightness Normalization al dataset
datos_co <- brightness_normalization(datos3, variables_co)

# Crear fórmula dinámica para Co
formula_co <- as.formula(paste0("`Co(g/t)` ~ ", paste(variables_co, collapse = " + ")))

# Configuración de validación cruzada
folds <- createFolds(datos_co$`Co(g/t)`, k = k)

# Crear modelo PLSR
modelo_co <- plsr(formula_co, data = datos_co, validation = "CV", segments = k)

# Selección del número óptimo de componentes
ncomp_opt_co <- max(1, selectNcomp(modelo_co, method = "onesigma", plot = FALSE))
cat("Número óptimo de componentes para Co:", ncomp_opt_co, "\n")

# Resultados del modelo
list(modelo = modelo_co, ncomp_opt = ncomp_opt_co)



################################################ Modelo Zn

set.seed(123)
k <- 10

# Variables predictoras para Zn
variables_zn <- c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi")

# Función de Standard Normalization
standard_normalization <- function(df, vars) {
  df[, vars] <- scale(df[, vars])  # Centrar y escalar
  return(df)
}

# Aplicar Standard Normalization al dataset
datos_zn <- standard_normalization(datos3, variables_zn)

# Crear fórmula dinámica para Zn
formula_zn <- as.formula(paste0("`Zn(g/t)` ~ ", paste(variables_zn, collapse = " + ")))

# Configuración de validación cruzada
folds <- createFolds(datos_zn$`Zn(g/t)`, k = k)

# Crear modelo PLSR
modelo_zn <- plsr(formula_zn, data = datos_zn, validation = "CV", segments = k)

# Selección del número óptimo de componentes
ncomp_opt_zn <- max(1, selectNcomp(modelo_zn, method = "onesigma", plot = FALSE))
cat("Número óptimo de componentes para Zn:", ncomp_opt_zn, "\n")

# Resultados del modelo
list(modelo = modelo_zn, ncomp_opt = ncomp_opt_zn)


################################ Gráfico de valores predichos vs reales

# Gráfico con cuadro ajustado y simbología gráfica en lugar de texto
graficar_pred_vs_real <- function(modelo, ncomp_opt, datos, metal, texto_personalizado) {
  # Predicciones con el modelo
  pred <- predict(modelo, ncomp = ncomp_opt)
  
  # Crear DataFrame con valores reales y predichos
  data_plot <- data.frame(
    Real = datos[[metal]],
    Predicho = pred[, 1, 1]
  )
  
  # Calcular el límite máximo
  max_val <- max(max(data_plot$Real, na.rm = TRUE), max(data_plot$Predicho, na.rm = TRUE)) * 1.1
  
  # Calcular posición del cuadro
  cuadro_x <- 0.05 * max_val
  cuadro_y <- 0.95 * max_val
  
  # Graficar
  ggplot(data_plot, aes(x = Real, y = Predicho)) +
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
    annotate("rect", xmin = cuadro_x, xmax = cuadro_x + max_val * 0.4,
             ymin = cuadro_y - max_val * 0.25, ymax = cuadro_y * 1.05,  
             fill = "white", alpha = 0.8, color = "black") +
    annotate("text", x = cuadro_x + max_val * 0.02, y = cuadro_y - max_val * 0.1,
             label = texto_personalizado, hjust = 0, size = 6.2, color = "black") +
    # Agregar simbología de líneas dentro del cuadro
    annotate("segment", x = cuadro_x + max_val * 0.02, xend = cuadro_x + max_val * 0.1,
             y = cuadro_y - max_val * 0.3, yend = cuadro_y - max_val * 0.3, 
             color = "darkred", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x + max_val * 0.12, y = cuadro_y - max_val * 0.3,
             label = "Ajuste ideal", hjust = 0, size = 5, color = "black") +
    annotate("segment", x = cuadro_x + max_val * 0.02, xend = cuadro_x + max_val * 0.1,
             y = cuadro_y - max_val * 0.36, yend = cuadro_y - max_val * 0.36, 
             color = "green", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x + max_val * 0.12, y = cuadro_y - max_val * 0.36,
             label = "Tendencia lineal", hjust = 0, size = 5, color = "black")
}

# Texto personalizado para cada metal
textos_personalizados <- list(
  "Co(g/t)" = "R2: 0.106\nRMSE: 1.07\nBias: 0.326",
  "Zn(g/t)" = "R2: 0.109\nRMSE: 16.47\nBias: 0.576",
  "Ni(g/t)" = "R2: 0.349\nRMSE: 0.416\nBias: 0.199",
  "Cu(g/t)" = "R2: 0.069\nRMSE: 1.701\nBias: 0.517",
  "As(g/t)" = "R2: 0.164\nRMSE: 17.97\nBias: 0.640"
)



print(graficar_pred_vs_real(
  modelo_zn,
  ncomp_opt_zn,
  datos_zn,
  "Zn(g/t)",
  textos_personalizados[["Zn(g/t)"]]
))



# Generar gráficos para cada metal en la lista general
for (metal in names(metales)) {
  print(graficar_pred_vs_real(
    resultados_modelos[[metal]]$modelo,
    resultados_modelos[[metal]]$ncomp_opt,
    datos3,
    metal,
    textos_personalizados[[metal]]
  ))
}

# Gráficos específicos para Co y Zn
print(graficar_pred_vs_real(
  modelo_co,
  ncomp_opt_co,
  datos_co,
  "Co(g/t)",
  textos_personalizados[["Co(g/t)"]]
))




################################ Análisis de Sensibilidad

# Análisis de Sensibilidad
analizar_sensibilidad <- function(modelo, metal) {
  # Verificar si el modelo tiene loading.weights válidos
  if (!is.null(modelo$loading.weights)) {
    # Extraer los pesos del primer componente
    pesos <- as.data.frame(modelo$loading.weights[, 1, drop = FALSE])  # Solo el primer componente
    colnames(pesos) <- c("Importancia")
    pesos$Variable <- rownames(pesos)
    
    # Ordenar por la importancia de las variables
    pesos <- pesos[order(abs(pesos$Importancia), decreasing = TRUE), ]
    
    # Imprimir tabla
    cat("\nAnálisis de Sensibilidad para:", metal, "\n")
    print(pesos)
  } else {
    cat("\nNo se encontraron pesos de carga válidos para:", metal, "\n")
  }
}

# Generar análisis para cada metal
for (metal in names(metales)) {
  cat("\n### Análisis de Sensibilidad para", metal, "###\n")
  analizar_sensibilidad(resultados_modelos[[metal]]$modelo, metal)
}

# Para Co
cat("\n### Análisis de Sensibilidad para Co(g/t) ###\n")
analizar_sensibilidad(modelo_co, "Co(g/t)")

# Para Zn
cat("\n### Análisis de Sensibilidad para Zn(g/t) ###\n")
analizar_sensibilidad(modelo_zn, "Zn(g/t)")





################################ Análisis de Contribución de Variables

# Análisis de Contribución de Variables
graficar_contribucion <- function(modelo, metal) {
  # Verificar si el modelo tiene loading.weights válidos
  if (!is.null(modelo$loading.weights)) {
    # Extraer los pesos del primer componente
    pesos <- as.data.frame(modelo$loading.weights[, 1, drop = FALSE])  # Primer componente
    colnames(pesos) <- c("Contribucion")  # Renombrar columna
    
    # Agregar nombres de las variables
    pesos$Variable <- rownames(pesos)
    
    # Calcular contribución absoluta y ordenar
    pesos$Contribucion <- abs(pesos$Contribucion)
    pesos <- pesos[order(pesos$Contribucion, decreasing = TRUE), ]
    
    # Graficar contribución
    ggplot(pesos, aes(x = reorder(Variable, Contribucion), y = Contribucion)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      coord_flip() +
      labs(
        title = paste("Contribución de Variables para", metal),
        x = "Variable",
        y = "Contribución Absoluta"
      ) +
      theme_minimal()
  } else {
    cat("\nNo se encontraron pesos de carga válidos para:", metal, "\n")
  }
}

# Generar gráficos para cada metal
for (metal in names(metales)) {
  print(graficar_contribucion(resultados_modelos[[metal]]$modelo, metal))
}

# Para Co
print(graficar_contribucion(modelo_co, "Co(g/t)"))

# Para Zn
print(graficar_contribucion(modelo_zn, "Zn(g/t)"))


####################################################### Modelos 

# Función para extraer y mostrar los coeficientes beta del modelo
mostrar_coeficientes <- function(modelo, ncomp_opt, metal) {
  # Extraer los coeficientes beta del modelo para el número óptimo de componentes
  betas <- coef(modelo, ncomp = ncomp_opt)
  
  # Convertir los coeficientes a un data frame
  coeficientes <- as.data.frame(betas)
  coeficientes$Variable <- rownames(coeficientes)  # Añadir nombres de las variables
  
  # Mostrar resultados
  cat("\nCoeficientes del Modelo PLSR para:", metal, "\n")
  print(coeficientes)
  
  # Mostrar la fórmula del modelo en formato simple
  cat("\nFórmula del Modelo:\n")
  cat("Y =", round(coeficientes[1, 1], 4))  # Intercepto
  for (i in 2:nrow(coeficientes)) {
    signo <- ifelse(coeficientes[i, 1] >= 0, "+", "-")
    cat(signo, abs(round(coeficientes[i, 1], 4)), "*", coeficientes$Variable[i], "")
  }
  cat("\n")
}

# Aplicar la función a cada modelo
for (metal in names(metales)) {
  mostrar_coeficientes(
    resultados_modelos[[metal]]$modelo,
    resultados_modelos[[metal]]$ncomp_opt,
    metal
  )
}

# Para Co
mostrar_coeficientes(modelo_co, ncomp_opt_co, "Co(g/t)")

# Para Zn
mostrar_coeficientes(modelo_zn, ncomp_opt_zn, "Zn(g/t)")




