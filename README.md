# TFM_CNV

## Revisión de un conjunto de herramientas bioinformáticas para la determinación de CNVs en muestras clínicas de exoma

#### Máster Universitario en Bioinformática y Biología Computacional
###### Autor: de Dios Blázquez, Lucía

### Resumen
La detección y caracterización precisa de las CNVs es esencial para comprender la variación humana y para los estudios de población, genética médica y evolución. A pesar de la importancia evolutiva y relevancia clínica de las CNVs, y de los avances en las técnicas de secuenciación y detección de variantes estructurales genómicas, los datos de caracterización de CNVs disponibles en repositorios públicos, son incompletos y presentan una alta tasa de falsos positivos. Por ello, la caracterización e interpretación de las CNVs continúa siendo un desafío, debido, además, a la dificultad de construir un conjunto de datos de validación de referencia para evaluar con precisión el rendimiento de las herramientas de determinación de variantes. En base a lo expuesto, nos propusimos dos objetivos principales en este trabajo: (i) obtener un gold standard de CNVs de referencia de la muestra control NA12878 (material de referencia del NIST), y (ii) evaluar el rendimiento de cuatro herramientas bioinformáticas independientes, con la finalidad de mejorar el análisis bioinformático de determinación de variantes en muestras clínicas de exoma. Mediante la recopilación de las variantes identificadas en estudios experimentales con diferentes diseños de arrays, generamos un conjunto de CNVs de referencia o gold standard para NA12878. Este conjunto de referencia incompleto pero robusto, incluye 1209 intervalos de CNV identificados en todo el genoma de la muestra de estudio, con una puntuación de confianza asignada para medir de una manear aproximada la incertidumbre asociada al intervalo. Para evaluar el rendimiento de las 4 herramientas de determinación de variantes en datos WES de NA12878, se utilizó el conjunto gold standard formado por 1652 exón-CNV y 14305 no exón-CNV. En este trabajo se analizaron las diferencias en cuanto a número, tipo y tamaños de las CNVs detectadas, así como diferentes métricas de precisión y rendimiento de predicción de cada algoritmo. Las herramientas mostraron poca concordancia y bajos rendimientos generales, que pueden explicarse por las limitaciones de este trabajo. Este estudio demuestra los desafíos a la hora de generar un conjunto de referencia completo y robusto de variaciones que abarquen el exoma humano; y las limitaciones de la predicción precisa de CNVs, a partir de los datos disponibles actualmente. A pesar de que es necesario un estudio más exhaustivo, las herramientas y metodologías utilizadas en este trabajo, así como el conjunto de datos de referencia generado, pueden ser de gran utilidad en la investigación de las CNVs con aplicaciones médicas, lo que permitirá una comprensión más completa del impacto funcional de estas variantes en las enfermedades. 



### Estructura de las carpetas

###### Objetivo 1. Generar un *gold standard* robusto de la muestra NA12878/HG001 del NIST para ser usado en las medidas de rendimiento de las herramientas bioinformáticas de determinación de CNVs en exoma
El conjunto de datos de los estudios iniciales de referencia, regiones repetitivas o de alta homología, datos de otros estudios y conjunto de los datos generados, así como los códigos utilizados en cada sección se encuentran en la carpeta: **gold_standard_NA12878**. Esta se encuentra dividida a su vez en cada apartado.

###### Objetivo 2. Evaluación de las medidas de rendimiento de las herramientas bioinformáticas de determinación de CNVs a partir del conjunto de CNVs de referencia generado para la muestra NA12878/HG001
Los *outputs* generados por cada herramienta de predicción de CNVs (a nivel de exón y a nivel de CNV), el código utilizado para la evaluación comparativa y los datos generados tras este análisis se encuentran en la carpeta **evaluation_cnvtools**.


