# Estimación de Concentraciones de Metales Pesados en Depósitos de Relaves usando Imágenes Satelitales Multiespectrales

## Descripción
Este repositorio contiene el código y los datos utilizados en la tesis de Magíster en Ciencia de Datos, titulada "Estimación de concentraciones de metales pesados en depósitos de relaves a partir de imágenes satelitales multiespectrales" realizada por Tomás Ignacio Contreras Delporte en la Universidad Adolfo Ibáñez. El trabajo se centra en la estimación de las concentraciones de metales pesados como Cu, Ni, Zn, Co y As en depósitos de relaves en la Región de Coquimbo, Chile, utilizando imágenes satelitales multiespectrales.

## Objetivo
El objetivo principal de este proyecto es evaluar la viabilidad de usar imágenes satelitales para estimar las concentraciones de metales pesados en depósitos de relaves mineros. A través de un análisis de regresión y el uso de diferentes modelos estadísticos y de machine learning, se busca identificar patrones en los datos espectrales de las imágenes para realizar estimaciones precisas.

## Metodología
1) *Área de Estudio:* Región de Coquimbo, Chile, con un enfoque en depósitos de relaves activos, inactivos y abandonados.

2) *Adquisición de Datos:*

  - Datos Satelitales: Imágenes multiespectrales obtenidas de satélites Landsat 8.

  - Datos de Concentraciones de Metales: Información histórica proporcionada por Sernageomin.

3) *Preprocesamiento de Datos:* Se realiza un preprocesamiento de las imágenes satelitales para asegurar la calidad de los datos espectrales, aplicando correcciones geométricas y radiométricas.

4) *Modelado Predictivo:* Uso de modelos estadísticos como Partial Least-Squares Regression (PLSR), Generalized Linear Models (GLM), y Generalized Additive Models (GAM) para estimar las concentraciones de metales pesados.

El modelo consiste en utilizar imágenes del satélite Landsat 8 para estimar las concentraciones de Cu, Ni, Zn, Co y As, en 228 depósitos de relaves en la región de Coquimbo, Chile. Para ello, se utilizan muestras tomadas por Sernageomin entre 2015 y 2017 como variables de respuesta. Luego, para estimar las concentraciones se utilizan 3 modelos: Partial Least Square Regression (PLSR), Generalized Linear Models (GLM) y Generalized Additive Models (GAM).
