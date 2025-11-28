# Computación Paralela: Estimación de Haplotipos 🧬

Este repositorio contiene la implementación y análisis de rendimiento de un algoritmo para la inferencia de haplotipos, comparando una versión secuencial contra una versión paralela optimizada utilizando `multiprocessing` en Python.

## 📋 Descripción

El proyecto aborda el problema combinatorio de determinar la configuración de haplotipos más probable para una población a partir de datos genotípicos. Se analiza la escalabilidad del algoritmo (Speedup y Eficiencia) al aumentar el número de individuos y posiciones genéticas, distribuyendo la carga de trabajo en múltiples núcleos de CPU.

## 📂 Estructura del Proyecto

El repositorio está organizado de la siguiente manera:

* **`00_bed_to_csv.ipynb`**: **Preprocesamiento**. Script que toma archivos crudos en formato `.ped` (PLINK), separa los alelos y genera archivos CSV estructurados listos para el análisis.
* **`01_algoritmo_sin_paralelizar.ipynb`**: **Prototipo Lógico**. Implementación secuencial (`base`) que define las funciones de selección de probabilidad y cálculo de frecuencias. Sirve para validación.
* **`02_script_paralelo.py`**: **Implementación Paralela**. Script principal que contiene la función `haplotipar`. Utiliza `multiprocessing.Pool` para paralelizar el cálculo de frecuencias por individuo.
* **`03_script_paralelo_con_tiempos.py`**: **Benchmark**. Script de automatización que ejecuta el algoritmo iterativamente variando:
    * *Procesos:* [1, 2, 4, 6, 8, 16, 32]
    * *Individuos:* [500, 1000, 1500, 3000]
    * *Posiciones:* [4, 6, 8, 10, 12]
    Los resultados se guardan en `resultados/tiempos.csv`.
* **`04_graficos.ipynb`**: **Visualización**. Notebook para el análisis de datos. Genera gráficos de tiempo de ejecución, eficiencia y escalabilidad utilizando `matplotlib` y `seaborn`.

## ⚙️ Requisitos

* Python 3.8+
* Librerías externas:
    ```bash
    pip install pandas numpy matplotlib seaborn
    ```
## 🚀 Uso

### 1. Ejecución del Algoritmo
Para utilizar la función principal en tu propio código:

```python
from 02_script_paralelo import haplotipar

# Parámetros:
# procesos: Número de núcleos a utilizar
# individuos: Cantidad de muestras a procesar
# posiciones: Número de marcadores genéticos
# resultado: Nombre del archivo de salida

haplotipar(procesos=4, individuos=1000, posiciones=8, resultado='mis_resultados.txt')
