# %%
import numpy as np
import pandas as pd
import itertools
import timeit
import multiprocessing
from functools import partial
from collections import Counter


# %%
NUMERO_DE_REPETICIONES = 10

import time

def monitor_tiempo_acumulado(N=NUMERO_DE_REPETICIONES):
    def decorador(func):
        def wrapper(*args, **kwargs):
            inicio = time.time()
            resultado = func(*args, **kwargs)
            fin = time.time()
            if wrapper.llamadas == N:
                wrapper.tiempo_total = 0
                wrapper.llamadas = 0
                wrapper.tiempo = []

            wrapper.tiempo_total += (fin - inicio)
            wrapper.tiempo.append(fin - inicio)
            wrapper.llamadas += 1

            
            return resultado

        wrapper.tiempo_total = 0
        wrapper.llamadas = 0
        wrapper.tiempo = []
        wrapper.__name__ = func.__name__
        return wrapper
    return decorador

def funcion_seleccion(combinaciones_con_prob):
    """ algoritmo de selección """
    combinaciones_con_prob = np.array(combinaciones_con_prob)
    indice = np.argmax(np.array(combinaciones_con_prob[:, 1], dtype=float) + np.array(combinaciones_con_prob[:, 3], dtype=float))
    return combinaciones_con_prob[indice]

def _worker_task(x, alelos, posiciones, combinaciones):
    """
    Función que procesa un solo individuo (x)
    y devuelve un Counter con sus frecuencias.
    """
    res_local = Counter()
    
    for i in range(len(combinaciones)):
        val = ''.join(alelos[x, np.arange(posiciones), combinaciones[i]])
        res_local[val] += 1 
        
    return res_local


@monitor_tiempo_acumulado()
def calcular_frecuencias(alelos, individuos, posiciones, combinaciones, procesos):
    """
    Calcular frecuencias de cada combinación en paralelo
    usando el número especificado de 'procesos'.
    """
    tarea_parcial = partial(_worker_task,
                            alelos=alelos,
                            posiciones=posiciones,
                            combinaciones=combinaciones)

    indices_individuos = range(individuos)

    with multiprocessing.Pool(processes=procesos) as pool:
        
        lista_de_counters = pool.map(tarea_parcial, indices_individuos)

    res_final = Counter()
    for res_local in lista_de_counters:
        res_final.update(res_local) # .update() suma las cuentas

    return dict(res_final)


def _worker_calcular_hoplotipo(x, posiciones, alelos, combinaciones, res, cantidad_total):
    
    combinaciones_con_prob = []
    for i in range(len(combinaciones)):
        aleloA_ind = combinaciones[i]
        aleloB_ind = np.abs(np.array(aleloA_ind) - 1)
        aleloA = ''.join(alelos[x, np.arange(posiciones), aleloA_ind])
        aleloB = ''.join(alelos[x, np.arange(posiciones), aleloB_ind])
        
        prob_A = res.get(aleloA, 0) / cantidad_total
        prob_B = res.get(aleloB, 0) / cantidad_total
        combinaciones_con_prob.append([aleloA, prob_A, aleloB, prob_B])
    
    return funcion_seleccion(combinaciones_con_prob)

@monitor_tiempo_acumulado()
def calcular_hoplotipo(individuos, posiciones, alelos, combinaciones, res, cantidad_total, procesos):
    """
    Calcular haplotipo más frecuente para cada individuo (en paralelo).
    """
    tarea_parcial = partial(_worker_calcular_hoplotipo,
                            posiciones=posiciones,
                            alelos=alelos,
                            combinaciones=combinaciones,
                            res=res,
                            cantidad_total=cantidad_total)

    indices_individuos = range(individuos)

    with multiprocessing.Pool(processes=procesos) as pool:
        results = pool.map(tarea_parcial, indices_individuos)

    return results


# %%
def haplotipar(procesos=1, individuos=10, posiciones=4, resultado='resultado_script.txt'):

    alelo_uno = pd.read_csv('data_ch6/chr6_allele1_wmeta.csv', nrows=individuos)
    alelo_dos = pd.read_csv('data_ch6/chr6_allele2_wmeta.csv', nrows=individuos)

    alelo_uno = alelo_uno[alelo_uno.columns[2:posiciones+2]].to_numpy()
    alelo_dos = alelo_dos[alelo_dos.columns[2:posiciones+2]].to_numpy()

    alelos = np.concatenate([alelo_uno.reshape(individuos, posiciones, 1), alelo_dos.reshape(individuos, posiciones, 1)], axis=2)


    combinaciones = list(itertools.product([0, 1], repeat=posiciones))

    res = calcular_frecuencias(alelos, individuos, posiciones, combinaciones, procesos)
    cantidad_total = len(combinaciones) * individuos


    results = calcular_hoplotipo(individuos, posiciones, alelos, combinaciones, res, cantidad_total, procesos)
    
    with open(resultado, 'w') as f:
        for i, r in enumerate(results):
            f.write(f"Individuo {i}: Alelo A: {r[0]}, Alelo B: {r[2]}, Prob A: {float(r[1]):.4f}, Prob B: {float(r[3]):.4f}\n")

# %%
resultado_tiempos_csv = []
def main(procesos, individuos=10, posiciones=4, numero_de_repeticiones=10, resultado='resultado.txt'):
    
    tiempos = timeit.repeat(
        lambda: haplotipar(procesos, individuos, posiciones, resultado),
        repeat=numero_de_repeticiones,
        number=1,
        globals=globals()
    )
    tiempos = np.array(tiempos)
    tiempo_promedio = tiempos.mean()
    tiempo_std = tiempos.std()
    
    
    np.array(calcular_hoplotipo.tiempo)
    
    tiempos_calcular_frecuencias = np.array(calcular_frecuencias.tiempo)
    tiempos_calcular_haplotipo = np.array(calcular_hoplotipo.tiempo)
    
    tiempo_promedio_calcular_frecuencias = tiempos_calcular_frecuencias.mean()
    tiempo_promedio_calcular_hoplotipo = tiempos_calcular_haplotipo.mean()
    
    tiempo_std_calcular_hoplotipo = tiempos_calcular_haplotipo.std()
    tiempo_std_calcular_frecuencias = tiempos_calcular_frecuencias.std()
    
    print(f"calcular_frecuencias -> Tiempo promedio por ejecución ({tiempos_calcular_frecuencias.shape[0]} repeticiones):: {tiempo_promedio_calcular_frecuencias:.6f} seg +- {tiempo_std_calcular_frecuencias:4f}")
    print(f"calcular_hoplotipo -> Tiempo promedio por ejecución ({tiempos_calcular_haplotipo.shape[0]} repeticiones): {tiempo_promedio_calcular_hoplotipo:.6f} seg +- {tiempo_std_calcular_hoplotipo:4f}")
    print(f"Tiempo promedio por ejecución ({numero_de_repeticiones} repeticiones): {tiempo_promedio:.4f} seg +- {tiempo_std:4f}")
    
    resultado_tiempos_csv.append([procesos, 
                                  individuos, 
                                  posiciones, 
                                  numero_de_repeticiones, 
                                  tiempo_promedio,
                                  tiempo_std, 
                                  tiempo_promedio_calcular_frecuencias,
                                  tiempo_std_calcular_frecuencias,
                                  tiempo_promedio_calcular_hoplotipo,
                                  tiempo_std_calcular_hoplotipo,])



if __name__ == "__main__":
    
    numeros_de_procesos = [1, 2, 3, 4, 6, 8, 16, 32]
    individuos_list = [500, 1000, 1500, 3000]
    posiciones_list = [4, 6, 8, 10, 12]

    for procesos, individuos, posiciones in itertools.product(numeros_de_procesos, individuos_list, posiciones_list):
        print("*"*20)
        print(f"Cantidad de procesos: {procesos}, cantidad de individuos: {individuos}, cantidad de posiciones: {posiciones}, cantidad de combinaciones posibles: {2 ** posiciones}")
        main(procesos, individuos, posiciones, NUMERO_DE_REPETICIONES,)
    
    print('Guardando tiempos...')
    resultado_tiempos = pd.DataFrame(resultado_tiempos_csv, 
                 columns=['Procesos', 
                          'Individuos', 
                          'Posiciones', 
                          'Repeticiones', 
                          'Tiempo Promedio',
                          'Tiempo Promedio std', 
                          'Tiempo Promedio Calcular Frecuencias', 
                          'Tiempo Promedio Calcular Frecuencias std', 
                          'Tiempo Promedio Calcular Hoplotipo',
                          'Tiempo Promedio Calcular Hoplotipo std'
                          ])
    resultado_tiempos.to_csv('tiempos.csv', index=False)
    print('Listo.')
# %%
