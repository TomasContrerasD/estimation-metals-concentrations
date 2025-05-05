
library(caret)
library(MASS)
library(readxl)
library(dplyr)
library(pls)
library(mgcv)
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

library(moments)
skewness(datos3$`Ni(g/t)`)
skewness(datos3$`Cu(g/t)`)
skewness(datos3$`Zn(g/t)`)
skewness(datos3$`Co(g/t)`)
skewness(datos3$`As(g/t)`)




# Tratar la negatividad de las variables predictoras

variables_predictoras <- c("ndvi", "dvi", "rvi", "sbi", "evi", "ndwi", "blue", 
                           "green", "red", "nir08", "coastal", "swir16", "swir22")

# Sumar 1 a todas las columnas de las variables predictoras
datos3[variables_predictoras] <- datos3[variables_predictoras] + 1


############################################################## 3) Creación de Variables Dummies con GLM


############################################################## 3.1) Por tamaño

# Calcular los cuartiles
cuartiles <- quantile(datos3$Area_m2, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# Crear las dummies por cuartil y la variable Cuartil_Tamaño
datos3 <- datos3 %>%
  mutate(
    Q1 = ifelse(Area_m2 <= cuartiles[2], 1, 0),
    Q2 = ifelse(Area_m2 > cuartiles[2] & Area_m2 <= cuartiles[3], 1, 0),
    Q3 = ifelse(Area_m2 > cuartiles[3] & Area_m2 <= cuartiles[4], 1, 0),
    Q4 = ifelse(Area_m2 > cuartiles[4], 1, 0),
    Cuartil_Tamaño = case_when(
      Area_m2 <= cuartiles[2] ~ "Q1",
      Area_m2 > cuartiles[2] & Area_m2 <= cuartiles[3] ~ "Q2",
      Area_m2 > cuartiles[3] & Area_m2 <= cuartiles[4] ~ "Q3",
      Area_m2 > cuartiles[4] ~ "Q4"
    )
  )


############################################################## 3.2) Por estado

# Crear dummies para la variable 'Estado'
dummies_estado <- model.matrix(~ Estado - 1, data = datos3) %>%
  as.data.frame()
colnames(dummies_estado) <- gsub("Estado", "", colnames(dummies_estado))

# Añadir las dummies al dataset
datos3 <- cbind(datos3, dummies_estado)

########################################################### 3.3) Normalizar Area_m2

datos3$Area_m2 <- (datos3$Area_m2 - min(datos3$Area_m2, na.rm = TRUE)) / 
  (max(datos3$Area_m2, na.rm = TRUE) - min(datos3$Area_m2, na.rm = TRUE))


########################################################### Analizar nuevas variables para PLSR


set.seed(123)
k <- 10

# Definir combinaciones de variables predictoras
combinaciones <- list(
  list(base = c("green", "nir08", "swir16", "ndvi", "ndwi"), metal = "Ni(g/t)"),
  list(base = c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "ndvi", "rvi", "sbi", "ndwi"), metal = "Cu(g/t)"),
  list(base = c("red", "nir08", "coastal", "swir16", "swir22", "rvi", "sbi"), metal = "As(g/t)")
)

# Crear combinaciones adicionales
definir_combinaciones <- function(base_vars, extra_vars) {
  list(
    base_vars,
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO"),
    c(base_vars, "Q1","Q2", "Q3", "Q4"),
    c(base_vars, "Area_m2"),
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO", "Area_m2")
  )
}

# DataFrame para almacenar resultados
resultados_globales <- data.frame(Metal = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por cada metal y combinación
for (config in combinaciones) {
  metal <- config$metal
  base_vars <- config$base
  
  # Generar combinaciones de variables
  extra_vars <- definir_combinaciones(base_vars, c("ABANDONADO", "INACTIVO", "ACTIVO", "Q1","Q2", "Q3", "Q4", "Area_m2"))
  
  cat("\n### Análisis para:", metal, "###\n")
  
  # Validación cruzada
  folds <- createFolds(datos3[[metal]], k = k)
  
  # Evaluar cada combinación de variables
  for (vars in extra_vars) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Crear fórmula dinámica
      formula <- as.formula(paste0("`", metal, "` ~ ", paste(vars, collapse = " + ")))
      
      # Ajustar el modelo PLSR
      modelo_pls <- tryCatch({
        plsr(formula, data = TRAIN, validation = "CV", segments = k)
      }, error = function(e) {
        message(sprintf("Error en el ajuste PLSR para la combinación: %s. Error: %s", paste(vars, collapse = ", "), e$message))
        return(NULL)
      })
      
      if (is.null(modelo_pls)) next
      
      # Selección del número óptimo de componentes
      ncomp_opt <- tryCatch({
        max(1, selectNcomp(modelo_pls, method = "onesigma", plot = FALSE))
      }, error = function(e) {
        message("Error al seleccionar el número óptimo de componentes: ", e$message)
        return(NA)
      })
      
      if (is.na(ncomp_opt)) next
      
      # Predicciones
      pred <- tryCatch({
        predict(modelo_pls, newdata = VALIDATION, ncomp = ncomp_opt)
      }, error = function(e) {
        message("Error en la predicción: ", e$message)
        return(rep(NA, nrow(VALIDATION)))
      })
      
      obs <- VALIDATION[[metal]]
      
      # Calcular métricas
      if (!all(is.na(pred))) {
        pred <- pred[, 1, 1]  # Extraer la dimensión de predicción
        r2[i] <- cor(pred, obs, use = "complete.obs")^2
        rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        bias[i] <- mean(pred - obs, na.rm = TRUE)
      } else {
        r2[i] <- NA
        rmse[i] <- NA
        bias[i] <- NA
      }
    }
    
    # Promedio de métricas para la combinación
    resultados_globales <- rbind(resultados_globales, data.frame(
      Metal = metal,
      Combinacion = paste(vars, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
}

# Eliminar combinaciones con resultados no válidos
resultados_globales <- resultados_globales[!is.na(resultados_globales$R2), ]

print(resultados_globales)



################################# Co y Zn

set.seed(123)
k <- 10

# Configuración de metales con normalización y combinaciones de variables
metales <- list(
  "Co(g/t)" = list(
    base_vars = c("blue", "nir08", "coastal", "swir16", "ndvi", "dvi", "evi", "ndwi"),
    normalization = function(df, vars) {
      brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
      brightness[brightness == 0] <- 1  # Evitar divisiones por 0
      df[, vars] <- df[, vars] / brightness
      return(df)
    }
  ),
  "Zn(g/t)" = list(
    base_vars = c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi"),
    normalization = function(df, vars) {
      df[, vars] <- scale(df[, vars])  # Centrar y escalar
      return(df)
    }
  )
)

# Función para definir combinaciones adicionales de variables
definir_combinaciones <- function(base_vars) {
  list(
    base_vars,
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO"),
    c(base_vars, "Q2", "Q3", "Q4"),
    c(base_vars, "Area_m2"),
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO", "Area_m2")
  )
}

# DataFrame para almacenar resultados
resultados_globales <- data.frame(Metal = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por cada metal
for (metal in names(metales)) {
  cat("\n### Análisis para:", metal, "###\n")
  
  config <- metales[[metal]]
  base_vars <- config$base_vars
  normalization <- config$normalization
  
  # Aplicar la normalización al dataset
  datos_filtrados <- normalization(datos3, base_vars)
  
  # Generar combinaciones de variables
  extra_vars <- definir_combinaciones(base_vars)
  
  # Validación cruzada
  folds <- createFolds(datos_filtrados[[metal]], k = k)
  
  # Evaluar cada combinación de variables
  for (vars in extra_vars) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos_filtrados[-folds[[i]], ]
      VALIDATION <- datos_filtrados[folds[[i]], ]
      
      # Crear fórmula dinámica
      formula <- as.formula(paste0("`", metal, "` ~ ", paste(vars, collapse = " + ")))
      
      # Ajustar el modelo PLSR
      modelo_pls <- tryCatch({
        plsr(formula, data = TRAIN, validation = "CV", segments = k)
      }, error = function(e) {
        message(sprintf("Error en el ajuste PLSR para la combinación: %s. Error: %s", paste(vars, collapse = ", "), e$message))
        return(NULL)
      })
      
      if (is.null(modelo_pls)) next
      
      # Selección del número óptimo de componentes
      ncomp_opt <- tryCatch({
        max(1, selectNcomp(modelo_pls, method = "onesigma", plot = FALSE))
      }, error = function(e) {
        message("Error al seleccionar el número óptimo de componentes: ", e$message)
        return(NA)
      })
      
      if (is.na(ncomp_opt)) next
      
      # Predicciones
      pred <- tryCatch({
        predict(modelo_pls, newdata = VALIDATION, ncomp = ncomp_opt)
      }, error = function(e) {
        message("Error en la predicción: ", e$message)
        return(rep(NA, nrow(VALIDATION)))
      })
      
      obs <- VALIDATION[[metal]]
      
      # Calcular métricas
      if (!all(is.na(pred))) {
        pred <- pred[, 1, 1]  # Extraer la dimensión de predicción
        r2[i] <- cor(pred, obs, use = "complete.obs")^2
        rmse[i] <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        bias[i] <- mean(pred - obs, na.rm = TRUE)
      } else {
        r2[i] <- NA
        rmse[i] <- NA
        bias[i] <- NA
      }
    }
    
    # Promedio de métricas para la combinación
    resultados_globales <- rbind(resultados_globales, data.frame(
      Metal = metal,
      Combinacion = paste(vars, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
}

# Eliminar combinaciones con resultados no válidos
resultados_globales <- resultados_globales[!is.na(resultados_globales$R2), ]

print(resultados_globales)


########################################################### Analizar nuevas variables para GLM

set.seed(123)
k <- 10

# Definir combinaciones de variables predictoras
combinaciones <- list(
  list(base = c("green", "nir08", "rvi", "sbi", "evi", "blue", "coastal", "swir16", "swir22"), 
       metal = "Ni(g/t)", 
       family = gaussian(link = "identity"),
       transform_predictors = FALSE),
  list(base = c("nir08", "swir16", "blue", "green", "coastal", "swir22"), 
       metal = "Co(g/t)", 
       family = Gamma(link = "log"),
       transform_predictors = TRUE)
)

# Crear combinaciones adicionales
definir_combinaciones <- function(base_vars) {
  list(
    base_vars,
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO"),
    c(base_vars, "Q1", "Q2", "Q3", "Q4"),
    c(base_vars, "Area_m2"),
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO", "Area_m2")
  )
}

# DataFrame para almacenar resultados
resultados_globales <- data.frame(Metal = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por cada metal y combinación
for (config in combinaciones) {
  metal <- config$metal
  base_vars <- config$base
  family <- config$family
  transform_predictors <- config$transform_predictors
  
  # Generar combinaciones de variables
  extra_vars <- definir_combinaciones(base_vars)
  
  cat("\n### Análisis para:", metal, "###\n")
  
  # Validación cruzada
  folds <- createFolds(datos3[[metal]], k = k)
  
  # Evaluar cada combinación de variables
  for (vars in extra_vars) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Aplicar transformación logarítmica a las variables predictoras si corresponde
      if (transform_predictors) {
        for (var in vars) {
          if (var %in% names(TRAIN)) {
            TRAIN[[var]] <- log1p(TRAIN[[var]])
            VALIDATION[[var]] <- log1p(VALIDATION[[var]])
          }
        }
      }
      
      # Crear fórmula dinámica
      formula <- as.formula(paste0("`", metal, "` ~ ", paste(vars, collapse = " + ")))
      
      # Ajustar el modelo GLM
      modelo_glm <- tryCatch({
        glm(formula, data = TRAIN, family = family)
      }, error = function(e) {
        message(sprintf("Error en el ajuste GLM para la combinación: %s. Error: %s", paste(vars, collapse = ", "), e$message))
        return(NULL)
      })
      
      if (is.null(modelo_glm)) next
      
      # Predicciones
      pred <- tryCatch({
        predict(modelo_glm, newdata = VALIDATION, type = "response")
      }, error = function(e) {
        message("Error en la predicción: ", e$message)
        return(rep(NA, nrow(VALIDATION)))
      })
      
      obs <- VALIDATION[[metal]]
      
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
    
    # Promedio de métricas para la combinación
    resultados_globales <- rbind(resultados_globales, data.frame(
      Metal = metal,
      Combinacion = paste(vars, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
}

# Eliminar combinaciones con resultados no válidos
resultados_globales <- resultados_globales[!is.na(resultados_globales$R2), ]

print(resultados_globales)



##################### Co con GAM

library(mgcv)

datos3 <- datos3 %>%
  rename(Co = `Co(g/t)`)

set.seed(123)
k <- 10

# Definir combinaciones de variables predictoras
combinaciones <- list(
  list(base = c("nir08", "swir16", "blue", "green", "coastal", "swir22"), 
       metal = "Co", 
       family = Gamma(link = "log"),
       transform_predictors = TRUE)
)

# Crear combinaciones adicionales
definir_combinaciones <- function(base_vars) {
  list(
    base_vars,
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO"),
    c(base_vars, "Q1", "Q2", "Q3", "Q4"),
    c(base_vars, "Area_m2"),
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO", "Area_m2")
  )
}

# Función para aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, vars] <- df[, vars] / brightness
  return(df)
}

# DataFrame para almacenar resultados
resultados_globales <- data.frame(Metal = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por cada metal y combinación
for (config in combinaciones) {
  metal <- config$metal
  base_vars <- config$base
  family <- config$family
  transform_predictors <- config$transform_predictors
  
  # Generar combinaciones de variables
  extra_vars <- definir_combinaciones(base_vars)
  
  cat("\n### Análisis para:", metal, "###\n")
  
  # Validación cruzada
  folds <- createFolds(datos3[[metal]], k = k)
  
  # Evaluar cada combinación de variables
  for (vars in extra_vars) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Aplicar Brightness Normalization a las variables predictoras base si corresponde
      if (transform_predictors) {
        TRAIN <- brightness_normalization(TRAIN, base_vars)
        VALIDATION <- brightness_normalization(VALIDATION, base_vars)
      }
      
      # Crear fórmula dinámica
      base_vars_s <- paste(paste0("s(", base_vars, ")"), collapse = " + ")  # Aplicar s() a variables base
      additional_vars <- setdiff(vars, base_vars)  # Variables adicionales sin s()
      additional_vars_formula <- if (length(additional_vars) > 0) {
        paste(additional_vars, collapse = " + ")
      } else {
        ""
      }
      
      # Combinar fórmula
      formula_parts <- c(base_vars_s, additional_vars_formula)
      formula <- as.formula(paste0("`", metal, "` ~ ", paste(formula_parts[formula_parts != ""], collapse = " + ")))
      
      # Ajustar el modelo GAM
      modelo_gam <- tryCatch({
        gam(formula, data = TRAIN, family = family)
      }, error = function(e) {
        message(sprintf("Error en el ajuste GAM para la combinación: %s. Error: %s", paste(vars, collapse = ", "), e$message))
        return(NULL)
      })
      
      if (is.null(modelo_gam)) next
      
      # Predicciones
      pred <- tryCatch({
        predict(modelo_gam, newdata = VALIDATION, type = "response")
      }, error = function(e) {
        message("Error en la predicción: ", e$message)
        return(rep(NA, nrow(VALIDATION)))
      })
      
      obs <- VALIDATION[[metal]]
      
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
    
    # Promedio de métricas para la combinación
    resultados_globales <- rbind(resultados_globales, data.frame(
      Metal = metal,
      Combinacion = paste(vars, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
}

# Eliminar combinaciones con resultados no válidos
resultados_globales <- resultados_globales[!is.na(resultados_globales$R2), ]

print(resultados_globales)



################# Cu, Zn y As


set.seed(123)
k <- 10

# Configuración de variables y fórmulas
configuracion_cu_zn_as <- list(
  "Cu(g/t)" = list(
    formula = as.formula("log1p(`Cu(g/t)`) ~ nir08 + swir16 + swir22 + evi + ndwi"),
    transform_response = TRUE
  ),
  "Zn(g/t)" = list(
    formula = as.formula("log1p(`Zn(g/t)`) ~ log1p(blue) + log1p(green) + log1p(red) + log1p(nir08) + log1p(coastal) + log1p(swir16) + log1p(swir22) + log1p(rvi) + log1p(sbi) + log1p(ndwi)"),
    transform_response = FALSE
  ),
  "As(g/t)" = list(
    formula = as.formula("log1p(`As(g/t)`) ~ nir08 + swir16 + ndvi + dvi + rvi + sbi"),
    transform_response = TRUE
  )
)

# Crear combinaciones de variables predictoras
definir_combinaciones <- function(base_vars) {
  list(
    base_vars,
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO"),
    c(base_vars, "Q1", "Q2", "Q3", "Q4"),
    c(base_vars, "Area_m2"),
    c(base_vars, "ABANDONADO", "INACTIVO", "ACTIVO", "Area_m2")
  )
}

# DataFrame para almacenar resultados
resultados_globales <- data.frame(Metal = character(), Combinacion = character(), R2 = numeric(), RMSE = numeric(), BIAS = numeric())

# Loop por cada metal y combinación
for (metal in names(configuracion_cu_zn_as)) {
  config <- configuracion_cu_zn_as[[metal]]
  formula <- config$formula
  transform_response <- config$transform_response
  
  # Extraer variables base
  base_vars <- all.vars(formula)[-1]  # Excluir la variable de respuesta
  
  # Generar combinaciones de variables
  extra_vars <- definir_combinaciones(base_vars)
  
  cat("\n### Análisis para:", metal, "###\n")
  
  # Validación cruzada
  folds <- createFolds(datos3[[metal]], k = k)
  
  # Evaluar cada combinación de variables
  for (vars in extra_vars) {
    r2 <- c()
    rmse <- c()
    bias <- c()
    
    for (i in 1:k) {
      TRAIN <- datos3[-folds[[i]], ]
      VALIDATION <- datos3[folds[[i]], ]
      
      # Crear fórmula dinámica
      formula_dynamic <- as.formula(paste0("`", metal, "` ~ ", paste(vars, collapse = " + ")))
      
      # Transformar la respuesta si corresponde
      if (transform_response) {
        TRAIN[[metal]] <- log1p(TRAIN[[metal]])
        VALIDATION[[metal]] <- log1p(VALIDATION[[metal]])
      }
      
      # Ajustar el modelo GLM.NB
      modelo_glm_nb <- tryCatch({
        glm.nb(formula_dynamic, data = TRAIN)
      }, error = function(e) {
        message(sprintf("Error en el ajuste GLM.NB para la combinación: %s. Error: %s", paste(vars, collapse = ", "), e$message))
        return(NULL)
      })
      
      if (is.null(modelo_glm_nb)) next
      
      # Predicciones
      pred <- tryCatch({
        predict(modelo_glm_nb, newdata = VALIDATION, type = "response")
      }, error = function(e) {
        message("Error en la predicción: ", e$message)
        return(rep(NA, nrow(VALIDATION)))
      })
      
      # Deshacer transformación de respuesta si corresponde
      if (transform_response) {
        pred <- expm1(pred)
        VALIDATION[[metal]] <- expm1(VALIDATION[[metal]])
      }
      
      obs <- VALIDATION[[metal]]
      
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
    
    # Promedio de métricas para la combinación
    resultados_globales <- rbind(resultados_globales, data.frame(
      Metal = metal,
      Combinacion = paste(vars, collapse = ", "),
      R2 = mean(r2, na.rm = TRUE),
      RMSE = mean(rmse, na.rm = TRUE),
      BIAS = mean(bias, na.rm = TRUE)
    ))
  }
}

# Eliminar combinaciones con resultados no válidos
resultados_globales <- resultados_globales[!is.na(resultados_globales$R2), ]

print(resultados_globales)




############################################ Gráficos

################################################ Modelo Zn PLSR

# Configuración inicial
set.seed(123)
k <- 10

# Variables predictoras para Zn
variables_zn <- c("green", "blue", "red", "nir08", "coastal", "swir16", "swir22", "sbi", "ndwi","Q1", "Q2", "Q3", "Q4")

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


library(ggplot2)

# Gráfico de residuos estandarizados
graficar_residuos_estandarizados <- function(modelo, datos, metal, ncomp_opt) {
  # Predicciones
  pred <- predict(modelo, datos, ncomp = ncomp_opt)[, 1, 1]
  
  # Cálculo de residuos
  residuos <- datos[[metal]] - pred
  residuos_estandarizados <- scale(residuos, center = TRUE, scale = sd(residuos))
  
  # Gráfico de residuos
  plot(pred, residuos_estandarizados,
       xlab = "Valores Ajustados",
       ylab = "Residuos Estandarizados",
       pch = 20,
       main = paste("Residuos Estandarizados para", metal),
       cex.lab = 1.5,
       cex.main = 1.5)
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(pred, residuos_estandarizados), col = "blue", lwd = 2)
}

# Gráfico Q-Q
graficar_qq <- function(modelo, datos, metal, ncomp_opt) {
  # Predicciones
  pred <- predict(modelo, datos, ncomp = ncomp_opt)[, 1, 1]
  
  # Cálculo de residuos
  residuos <- datos[[metal]] - pred
  
  # Gráfico Q-Q
  qqnorm(residuos,
         main = paste("Gráfico Q-Q de Residuos para", metal),
         cex.main = 1.5,
         cex.lab = 1.5)
  qqline(residuos, col = "blue", lwd = 2)
}

# Gráfico de valores predichos vs reales
graficar_pred_vs_real <- function(modelo, datos, metal, ncomp_opt) {
  # Predicciones
  pred <- predict(modelo, datos, ncomp = ncomp_opt)[, 1, 1]
  
  # Crear DataFrame para ggplot
  data_plot <- data.frame(
    Real = datos[[metal]],
    Predicho = pred
  )
  
  # Calcular límite máximo
  max_val <- max(c(data_plot$Real, data_plot$Predicho), na.rm = TRUE) * 1.1
  
  # Gráfico
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

# Análisis de Sensibilidad
analizar_sensibilidad <- function(modelo, metal) {
  if (!is.null(modelo$loading.weights)) {
    # Extraer pesos de los componentes
    pesos <- as.data.frame(modelo$loading.weights[, 1, drop = FALSE])
    colnames(pesos) <- c("Importancia")
    pesos$Variable <- rownames(pesos)
    
    # Ordenar por importancia
    pesos <- pesos[order(abs(pesos$Importancia), decreasing = TRUE), ]
    
    # Imprimir tabla
    cat("\n### Análisis de Sensibilidad para", metal, "###\n")
    print(pesos)
  } else {
    cat("\nNo se encontraron pesos de carga válidos para:", metal, "\n")
  }
}

# Ejecutar los análisis para el modelo de Zinc
cat("\n--- Gráfico de Residuos Estandarizados para Zn ---\n")
graficar_residuos_estandarizados(modelo_zn, datos_zn, "Zn(g/t)", ncomp_opt_zn)

cat("\n--- Gráfico Q-Q para Zn ---\n")
graficar_qq(modelo_zn, datos_zn, "Zn(g/t)", ncomp_opt_zn)

cat("\n--- Gráfico de Valores Predichos vs Reales para Zn ---\n")
print(graficar_pred_vs_real(modelo_zn, datos_zn, "Zn(g/t)", ncomp_opt_zn))

cat("\n--- Análisis de Sensibilidad para Zn ---\n")
analizar_sensibilidad(modelo_zn, "Zn(g/t)")


################################################# Modelos Ni con GLM

set.seed(123)

# Configuración para Ni
configuracion_ni <- list(
  formula = as.formula("`Ni(g/t)` ~ green + nir08 + rvi + sbi + evi + blue + coastal + swir16 + swir22 + Area_m2"),
  family = gaussian(link = "identity")
)

# Lista para almacenar resultados
resultados_ni <- list()

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

# Análisis para Ni
cat("\nAnálisis para: Ni(g/t)\n")

# Extraer configuración
formula <- configuracion_ni$formula
family <- configuracion_ni$family

# Filtrar datos sin NA
variables <- all.vars(formula)
datos_filtrados <- na.omit(datos3[, variables])

# Crear folds para k-fold cross-validation
folds <- createFolds(datos_filtrados[[variables[1]]], k = k, list = TRUE)

# Lista para almacenar R2
r2_folds <- c()
predicciones_completas <- numeric(nrow(datos_filtrados))
residuos_completos <- numeric(nrow(datos_filtrados))

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
  
  # Guardar predicciones y residuos para el fold actual
  indices <- folds[[i]]
  predicciones_completas[indices] <- pred
  residuos_completos[indices] <- datos_filtrados[[variables[1]]][indices] - pred
  
  # Calcular R2
  r2 <- calcular_r2(VALIDATION[[variables[1]]], pred)
  r2_folds <- c(r2_folds, r2)
}

# Promedio de R2
r2_promedio <- mean(r2_folds, na.rm = TRUE)

# Guardar resultados
resultados_ni <- list(modelo = modelo, R2 = r2_promedio)

# Imprimir R2
cat("R2 promedio para Ni(g/t):", r2_promedio, "\n")

# Cálculo de residuos estandarizados
residuos_estandarizados <- scale(residuos_completos, center = TRUE, scale = sd(residuos_completos))

# Gráfico de residuos estandarizados
plot(predicciones_completas, residuos_estandarizados,
     main = "Níquel",
     xlab = "Valores Ajustados",
     ylab = "Residuos Estandarizados",
     pch = 20,
     cex.lab = 1.5,
     cex.axis = 1.5,    # Tamaño de los números en los ejes
     cex.main = 2.5)
abline(h = 0, col = "red", lwd = 2)
lines(lowess(predicciones_completas, residuos_estandarizados), col = "blue", lwd = 2)


# Función para graficar valores predichos vs reales
graficar_pred_vs_real <- function(modelo, datos, metal, texto_personalizado) {
  pred <- predict(modelo, newdata = datos, type = "response")
  
  data_plot <- data.frame(
    Real = datos[[metal]],
    Predicho = pred
  )
  
  # Calcular el límite máximo
  max_val <- max(c(data_plot$Real, data_plot$Predicho), na.rm = TRUE) * 1.1
  
  # Calcular posición del cuadro
  cuadro_x_min <- max_val * 0.05
  cuadro_x_max <- max_val * 0.3
  cuadro_y_min <- max_val * 0.8
  cuadro_y_max <- max_val * 0.9
  
  ggplot(data_plot, aes(x = Real, y = Predicho)) +
    geom_point(alpha = 0.7, color = "#1f77b4", size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred", size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "green", linetype = "dotted", size = 1.2) +
    labs(
      title = paste("Predicciones vs Observaciones para", metal),
      x = "Valores Observados",
      y = "Valores Predichos"
    ) +
    coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray", linetype = "dashed"),
      panel.grid.minor = element_blank()
    ) +
    # Agregar cuadro y texto personalizado
    annotate("rect", xmin = cuadro_x_min, xmax = cuadro_x_max,
             ymin = cuadro_y_min * 0.97, ymax = cuadro_y_max * 1.1,
             fill = "white", alpha = 0.8, color = "black") +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             y = cuadro_y_max - (cuadro_y_max - cuadro_y_min) * 0.15,
             label = texto_personalizado, hjust = 0, size = 5, color = "black") +
    # Agregar simbología de líneas dentro del cuadro
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -0.6, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -0.6,
             color = "darkred", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -0.6,
             label = "Ajuste ideal", hjust = 0, size = 4, color = "black") +
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1,
             color = "green", linetype = "dotted", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1,
             label = "Tendencia lineal", hjust = 0, size = 4, color = "black")
}

# Generar gráfico para Ni
modelo_ni <- resultados_ni$modelo
datos_filtrados_ni <- na.omit(datos3[, all.vars(configuracion_ni$formula)])

# Texto personalizado para el cuadro
texto_personalizado_ni <- "R2: 0.397\nRMSE: 0.397\nBias: 0.087"

# Generar gráfico para Ni
print(graficar_pred_vs_real(modelo_ni, datos_filtrados_ni, "Ni(g/t)", texto_personalizado_ni))



################################## Modelo Co con GAM


# Renombrar la columna
datos3 <- datos3 %>% rename(Co = `Co(g/t)`)

# Configuración del modelo
formula_co <- as.formula("Co ~ s(nir08) + s(swir16) + s(blue) + s(green) + s(coastal) + s(swir22) + ABANDONADO + INACTIVO + ACTIVO")
family_co <- Gamma(link = "log")

# Aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, vars] <- df[, vars] / brightness
  return(df)
}

# Filtrar y normalizar los datos
variables_predictoras <- c("nir08", "swir16", "blue", "green", "coastal", "swir22", "ABANDONADO", "INACTIVO", "ACTIVO")
datos_filtrados_co <- na.omit(datos3[, c("Co", variables_predictoras)])
datos_filtrados_co <- brightness_normalization(datos_filtrados_co, variables_predictoras)

# Ajustar el modelo GAM
modelo_co <- gam(formula_co, data = datos_filtrados_co, family = family_co)

# Calcular métricas promedio (opcional)
r2 <- cor(predict(modelo_co, type = "response"), datos_filtrados_co$Co)^2
rmse <- sqrt(mean((predict(modelo_co, type = "response") - datos_filtrados_co$Co)^2))
bias <- mean(predict(modelo_co, type = "response") - datos_filtrados_co$Co)
texto_personalizado_co <- "R2: 0.218\nRMSE: 1.051\nBias: -0.295"

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
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "green", linetype = "dotted", size = 1.2) +
    labs(
      title = paste("Predicciones vs Observaciones para", metal),
      x = "Valores Observados",
      y = "Valores Predichos"
    ) +
    coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray", linetype = "dashed"),
      panel.grid.minor = element_blank()
    ) +
    # Agregar cuadro y texto personalizado
    annotate("rect", xmin = cuadro_x_min, xmax = cuadro_x_max,
             ymin = cuadro_y_min * 0.97, ymax = cuadro_y_max * 1.1,
             fill = "white", alpha = 0.8, color = "black") +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             y = cuadro_y_max - (cuadro_y_max - cuadro_y_min) * 0.15,
             label = texto_personalizado, hjust = 0, size = 5, color = "black") +
    # Agregar simbología de líneas dentro del cuadro
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -0.6, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -0.6,
             color = "darkred", linetype = "dashed", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -0.6,
             label = "Ajuste ideal", hjust = 0, size = 4, color = "black") +
    annotate("segment", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.05, 
             xend = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.2,
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1, 
             yend = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1,
             color = "green", linetype = "dotted", size = 1.2) +
    annotate("text", x = cuadro_x_min + (cuadro_x_max - cuadro_x_min) * 0.22, 
             y = cuadro_y_min + (cuadro_y_max - cuadro_y_min) * -1,
             label = "Tendencia lineal", hjust = 0, size = 4, color = "black")
}

# Generar gráfico para Co
print(graficar_pred_vs_real(modelo_co, datos_filtrados_co, "Co", texto_personalizado_co))



# Cálculo de residuos estandarizados
predicciones <- predict(modelo_co, type = "response")
residuos <- datos_filtrados_co$Co - predicciones
residuos_estandarizados <- scale(residuos)

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

############################## Análisis de Sensibilidad

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
        axis.title = element_text(size = 16, color = "black"),   # Aumentar tamaño y color negro para títulos
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

# Ajustar el modelo completo para Ni
modelo_ni <- glm(configuracion_ni$formula, data = datos_filtrados, family = configuracion_ni$family)

# Análisis de sensibilidad para Ni
analizar_sensibilidad(modelo_ni, "Ni(g/t)")

# Crear tabla de sensibilidad para Ni
tabla_sensibilidad_ni <- crear_tabla_sensibilidad(modelo_ni, "Ni(g/t)")

# Análisis de sensibilidad para Co
analizar_sensibilidad(modelo_co, "Co")

# Crear tabla de sensibilidad para Co
tabla_sensibilidad_co <- crear_tabla_sensibilidad(modelo_co, "Co")


############## Co

# Función modificada para graficar solo las 10 variables más importantes
analizar_sensibilidad_top10 <- function(modelo, metal) {
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
    
    # Seleccionar solo las 10 primeras variables
    coeficientes_top10 <- head(coeficientes_df, 10)
    
    # Graficar
    ggplot(coeficientes_top10, aes(x = reorder(Variable, abs(Importancia)), y = Importancia)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        x = "Variable Predictora",
        y = "Coeficiente (Importancia)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16, color = "black"),   # Aumentar tamaño y color negro para títulos
        axis.text = element_text(size = 14, color = "black") 
      )
  }
}

# Aplicar la función para el modelo Co
analizar_sensibilidad_top10(modelo_co, "Co")


########################################################## métricas por dummy

# Configuración del modelo
formula_co <- as.formula("Co ~ s(nir08) + s(swir16) + s(blue) + s(green) + s(coastal) + s(swir22) + ABANDONADO + INACTIVO + ACTIVO")
family_co <- Gamma(link = "log")

# Aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, vars] <- df[, vars] / brightness
  return(df)
}

# Filtrar y normalizar los datos
variables_predictoras <- c("nir08", "swir16", "blue", "green", "coastal", "swir22", "ABANDONADO", "INACTIVO", "ACTIVO")
datos_filtrados_co <- na.omit(datos3[, c("Co", "Estado", variables_predictoras)])  # Incluir "Estado"
datos_filtrados_co <- brightness_normalization(datos_filtrados_co, variables_predictoras)

# Ajustar el modelo GAM
modelo_co <- gam(formula_co, data = datos_filtrados_co, family = family_co)

# Función para calcular métricas
calcular_metricas <- function(obs, pred) {
  r2 <- if (sum(!is.na(obs) & !is.na(pred)) > 1) cor(pred, obs, use = "complete.obs")^2 else NA
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  bias <- mean(pred - obs, na.rm = TRUE)
  return(list(R2 = r2, RMSE = rmse, Bias = bias))
}

# Añadir predicciones al dataset
datos_filtrados_co$Predicho <- predict(modelo_co, newdata = datos_filtrados_co, type = "response")

# DataFrame para almacenar resultados por categoría
resultados_por_categoria <- data.frame(
  Categoria = character(),
  R2 = numeric(),
  RMSE = numeric(),
  Bias = numeric(),
  N = numeric(),  # Nueva columna para almacenar la cantidad de datos
  stringsAsFactors = FALSE
)

# Calcular métricas por categoría
for (categoria in unique(datos_filtrados_co$Estado)) {  # Usar "Estado" del dataset filtrado
  # Filtrar datos por categoría
  datos_categoria <- subset(datos_filtrados_co, Estado == categoria)
  
  # Verificar si la categoría tiene datos suficientes
  if (nrow(datos_categoria) > 1) {
    # Calcular métricas
    metricas <- calcular_metricas(datos_categoria$Co, datos_categoria$Predicho)
    resultados_por_categoria <- rbind(resultados_por_categoria, data.frame(
      Categoria = categoria,
      R2 = metricas$R2,
      RMSE = metricas$RMSE,
      Bias = metricas$Bias,
      N = nrow(datos_categoria)  # Almacenar la cantidad de datos
    ))
  }
}

print(resultados_por_categoria)



############################################## ANOVA Co con GAM y dummies

datos3 <- datos3 %>% rename(Co = `Co(g/t)`)
# Configuración del modelo
formula_co <- as.formula("Co ~ s(nir08) + s(swir16) + s(blue) + s(green) + s(coastal) + s(swir22) + ABANDONADO + INACTIVO + ACTIVO")
family_co <- Gamma(link = "log")

# Aplicar Brightness Normalization
brightness_normalization <- function(df, vars) {
  brightness <- sqrt(rowSums(df[, vars]^2, na.rm = TRUE))
  brightness[brightness == 0] <- 1  # Evitar divisiones por cero
  df[, vars] <- df[, vars] / brightness
  return(df)
}

# Filtrar y normalizar los datos
variables_predictoras <- c("nir08", "swir16", "blue", "green", "coastal", "swir22", "ABANDONADO", "INACTIVO", "ACTIVO")
datos_filtrados_co <- na.omit(datos3[, c("Co", variables_predictoras)])
datos_filtrados_co <- brightness_normalization(datos_filtrados_co, variables_predictoras)

# Ajustar el modelo GAM
modelo_co <- gam(formula_co, data = datos_filtrados_co, family = family_co)

# Realizar el análisis ANOVA
anova_resultado <- anova.gam(modelo_co)

print(anova_resultado)




######################################### Random Forest
install.packages("randomForest")


# Cargar las librerías necesarias
library(randomForest)
library(caret) # Para partición y métricas
library(dplyr)

# Eliminar filas con NA (opcional)
datos3 <- na.omit(datos3)

# Seleccionar las variables predictoras y la variable respuesta
predictoras <- datos3 %>% select(blue, green, red, nir08, ndvi, dvi, rvi, sbi, evi, ndwi)
respuesta <- datos3$`Cu(g/t)`

# Dividir los datos en entrenamiento (70%) y prueba (30%)
set.seed(123) # Asegura reproducibilidad
train_index <- createDataPartition(respuesta, p = 0.7, list = FALSE)

train_data <- datos3[train_index, ]
test_data <- datos3[-train_index, ]

# Crear el modelo de Random Forest
modelo_rf <- randomForest(
  `Cu(g/t)` ~ blue + green + red + nir08 + ndvi + dvi + rvi + sbi + evi + ndwi,
  data = train_data,
  importance = TRUE, # Para calcular la importancia de las variables
  ntree = 500, # Número de árboles
  mtry = 3 # Número de variables usadas en cada división
)

# Predecir en los datos de prueba
predicciones <- predict(modelo_rf, newdata = test_data)

# Calcular R2
r2 <- cor(predicciones, test_data$`Cu(g/t)`)^2
cat("El R² del modelo Random Forest es:", round(r2, 3), "\n")

# Importancia de las variables
importancia <- importance(modelo_rf)
importancia_ordenada <- importancia[order(importancia[, 1], decreasing = TRUE), ]
print(importancia_ordenada)

# Graficar la importancia de las variables
varImpPlot(modelo_rf, main = "Importancia de las Variables en el Modelo RF")










