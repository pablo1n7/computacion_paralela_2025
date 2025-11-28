# %%
import numpy as np
import pandas as pd
import itertools
import multiprocessing
from functools import partial
from collections import Counter


# %%

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
def haplotipar(procesos=1, individuos=10, posiciones=4, resultado='resultado.txt'):

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


if __name__ == "__main__":
    haplotipar(1, 10, 4, 'resultados/resultado_paralelo.txt')
    print('Listo.')
# %%
