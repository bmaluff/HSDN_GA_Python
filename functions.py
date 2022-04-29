# -*- coding: utf-8 -*-
'''
Created on 2 jul. 2019

@author: Bader Maluff
'''
from __future__ import division
import random
import sys
from collections import deque, namedtuple
import globales
from types import NoneType

#import statistics


def _generar_poblacion(cantidad_individuos, cantidad_dispositivos, cantidad_solicitudes):
    #Devuelve un cromosoma con valores generados aleatoriamente
    cromo_solicitudes = [lista_random_solicitudes(cantidad_solicitudes) for i in range(cantidad_individuos)]
    cromo_dispositivos = [[random.randint(0,1) for i in range(cantidad_dispositivos)] for j in range(cantidad_individuos)]
    cromo_sdn_flows = [[random.randint(0,1) for i in range(cantidad_solicitudes)] for j in range(cantidad_individuos)]
    return  cromo_solicitudes, cromo_dispositivos, cromo_sdn_flows

def lista_random_solicitudes(cantidad_solicitudes):
	"""Devuelve un vector con la cantidad de elementos del parámetro pasado
	   el vector se rellena con valores no repetidos entre 1 y parámetro pasado
    """
	aux = random.sample([i for i in range(cantidad_solicitudes)], cantidad_solicitudes)
    return aux

def inserta(x, lst, i):
    #Devuelve una nueva lista resultado de insertar x dentro de lst en la posici�n i.
    return lst[:i] + [x] + lst[i:]

def inserta_multiple(x, lst):
    """Devuelve una lista con el resultado de
       insertar x en todas las posiciones de lst.  
    """
    return [inserta(x, lst, i) for i in range(len(lst) + 1)]

def permuta(c):
    """Calcula y devuelve una lista con todas las
       permutaciones posibles que se pueden hacer
       con los elementos contenidos en c.
    """
    if len(c) == 0:
        return [[]]
    return sum([inserta_multiple(c[0], s)
                for s in permuta(c[1:])],
               [])

def seleccion_padres_determinista(aptitudes, cantidad_individuos):
    
    x1 = random.randint(0, cantidad_individuos-1)
    x2 = random.randint(0, cantidad_individuos-1)
    if aptitudes[x1] > aptitudes[x2]:
        padre = x1
    else:
        padre = x2
    
    x1 = random.randint(0, cantidad_individuos-1)
    x2 = random.randint(0, cantidad_individuos-1)
    if aptitudes[x1] > aptitudes[x2]:
        madre = x1
    else:
        madre = x2
        
    return padre, madre

def cantidad_dimensiones(a):
    try:
        if len(a) > 0:
            return 1 + cantidad_dimensiones(a[0])
    except TypeError:
        return 0

def cruza_1_punto(pob, aptitudes, padre, madre,):
    dimensiones = cantidad_dimensiones(pob)
    if dimensiones > 2:
        hijo_1=[]
        hijo_2=[]
        punto_cruza = [random.randint(0, len(pob[i][0])-1) for i in range(len(pob))]
        i = 1
        while i < len(pob):
            hijo_1.append(pob[i][padre][:punto_cruza[i]]+pob[i][madre][punto_cruza[i]:])
            hijo_2.append(pob[i][madre][punto_cruza[i]:]+pob[i][padre][:punto_cruza[i]])
            i+=1
        
    else:
        punto_cruza = random.randint(0, len(pob[0])-1)
        hijo_1=pob[padre][:punto_cruza[i]]+pob[madre][:punto_cruza[i]]
        hijo_2=pob[madre][punto_cruza[i]+1:]+pob[padre][punto_cruza[i]+1:]
    
    return hijo_1, hijo_2

def mutacion_permutacion(individuo, tasa_ruido):
    '''
    Esta función recibe un individuo y realiza la mutación tantas veces como sea necesario
    hasta llegar al porcentaje de mutación deseado.
    '''
    cantidad_solicitudes=len(individuo)
    contador=0.0
    off_set=1/cantidad_solicitudes
    while contador < tasa_ruido:
        origen_insercion=random.randint(0,cantidad_solicitudes-1)
        destino_insercion=random.randint(0,cantidad_solicitudes-1)
        aux=individuo.pop(origen_insercion)
        individuo = inserta(aux, individuo, destino_insercion)
        contador = contador + off_set
    return individuo

def mutacion_flipbit(individuo, tasa_ruido):
#Es una funcion que cambia los bits de forma aleatoria basado en una tasa de ruido mínima    
    cantidad_dispositivos = len(individuo)
    contador=0
    while contador < tasa_ruido:
        dispositivo=random.randint(0,cantidad_dispositivos-1)
        individuo[dispositivo]=flip(individuo[dispositivo])
        contador += (1/cantidad_dispositivos)
    return individuo

def mutacion(poblacion, ruido_permutacion, ruido_flipbit, porcentaje_mutacion):
#Recibe la población completa, elige los individuos al azar hata llegar al porcentaje de mutación

    contador_global=0
    cantidad_poblacion=len(poblacion[0])
    if cantidad_poblacion > 1:
        while contador_global < porcentaje_mutacion:
            individuo=random.randint(1,cantidad_poblacion-1)
            poblacion[0][individuo]=mutacion_permutacion(poblacion[0][individuo], ruido_permutacion)
            poblacion[1][individuo]=mutacion_flipbit(poblacion[1][individuo], ruido_flipbit)
            poblacion[2][individuo]=mutacion_flipbit(poblacion[2][individuo], ruido_flipbit)
            contador_global+=(1/cantidad_poblacion)
    return poblacion

def sel_x_mut(poblacion, aptitudes, porcentaje_mutacion, tasa_ruido_permutacion, tasa_ruido_flipbit):
"""Es la función con los operadore evolutivos, para cada operador llama a una función diferente
	y da como resultado una nueva generación
"""
    cantidad_individuos=len(poblacion[0])
    nueva_poblacion = []
    for i in range(len(poblacion)):
        nueva_poblacion.append([])
        for j in range (len(poblacion[i])):
            nueva_poblacion[i].append([])
            for k in range(len(poblacion[i][j])):
                nueva_poblacion[i][j].append(None)
    #nueva_poblacion = poblacion[:]
    poblados = 0
    posi_mas_apto = obtener_pos_mas_apto(aptitudes)
    nueva_poblacion[0][poblados]=poblacion[0][posi_mas_apto]
    nueva_poblacion[1][poblados]=poblacion[1][posi_mas_apto]
    nueva_poblacion[2][poblados]=poblacion[2][posi_mas_apto]
    poblados +=1
    while poblados < len(nueva_poblacion[0]):
        padre, madre = seleccion_padres_determinista(aptitudes, cantidad_individuos)
        hijo_1, hijo_2 = cruza_1_punto(poblacion, aptitudes, padre, madre)
        hijo_ruta_1, hijo_ruta_2=OX(poblacion, padre, madre)
        nueva_poblacion[0][poblados] = hijo_ruta_1
        nueva_poblacion[1][poblados] = hijo_1[0]
        nueva_poblacion[2][poblados] = hijo_1[1]
        poblados +=1
        if poblados >= len(nueva_poblacion[0]):
            break
        nueva_poblacion[0][poblados] = hijo_ruta_2
        nueva_poblacion[1][poblados] = hijo_2[0]
        nueva_poblacion[2][poblados] = hijo_2[1]
        poblados +=1
    nueva_poblacion= mutacion(nueva_poblacion, tasa_ruido_permutacion, tasa_ruido_flipbit, porcentaje_mutacion)
        
    return nueva_poblacion

def OX(pob, padre, madre):
#esta es la función de cruza para permutaciones, el funcionamiento se describe en el libro Apuntes de Coello
#Recibe las ubicaciones de los padres dentro de la población y devuelve los 2 hijos
    indice_1=random.randint(0, len(pob[0][0])-2)
    indice_2=random.randint(0, len(pob[0][0])-1)
    while indice_1>=indice_2:
        indice_2=random.randint(0, len(pob[0][0])-1)
    h1=[]
    h2=[]
    contador = 0
    for i in range(len(pob[0][0])):
        if i >indice_1 and i<=indice_2:
            h1.append(pob[0][padre][i])
        else:
            while True:
                if pob[0][madre][contador] not in pob[0][padre][indice_1+1:indice_2+1]:
                    h1.append(pob[0][madre][contador])
                    contador += 1
                    break
                else:
                    contador += 1

    indice_1=random.randint(0, len(pob[0][0])-2)
    indice_2=random.randint(0, len(pob[0][0])-1)
    while indice_1>=indice_2:
        indice_2=random.randint(0, len(pob[0][0])-1)
    contador = 0
    
    for i in range(len(pob[0][0])):
        if i >indice_1 and i<=indice_2:
            h2.append(pob[0][madre][i])
        else:
            while True:
                if pob[0][padre][contador] not in pob[0][madre][indice_1+1:indice_2+1]:
                    h2.append(pob[0][padre][contador])
                    contador += 1
                    break
                else:
                    contador += 1
    return h1, h2

def flip(a):
    if a==0:
        return 1
    return 0

def obtener_mas_apto(poblacion, aptitudes):
    aux = aptitudes[0]
    aux_posi = 0
    cantidad_individuos=len(aptitudes)
    for i in range(len(aptitudes)-1):
        if aptitudes[i] > aux:
            aux =aptitudes[i]
            aux_posi = i 
    return poblacion[0][aux_posi], poblacion[1][aux_posi], poblacion[2][aux_posi]

def obtener_pos_mas_apto(aptitudes):
    aux = aptitudes[0]
    aux_posi = 0
    #cantidad_individuos=len(aptitudes)
    for i in range(len(aptitudes)):
        if aptitudes[i] < aux and aptitudes[i] != inf:
            aux =aptitudes[i]
            aux_posi = i
    return aux_posi

def calcular_costos_sfp_old(mat_adya):
    grafo = Graph_old(mat_adya)
    mat_sfp_costos = []
    for i in range(len(mat_adya)):
        mat_sfp_costos.append([])
        costos=grafo.dijkstra(i)
        for j in range (len(mat_adya[i])):
            mat_sfp_costos[i].append(costos[j])
            #mat_sfp_costos[i].append([i,j,costos[j]])
    return mat_sfp_costos

inf = float('inf')
Edge = namedtuple('Edge', 'start, end, cost')


def make_edge(start, end, cost=1):
  return Edge(start, end, cost)


class Graph:
"""	Es la clase que define la topología como un grafo
"""
    def __init__(self, *args):
        # let's check that the data is right
        if len(args[0]) > 0 :
            
            if len(args[0][0]) > 3:
                self.edges=[]
                for i in range(len(args[0])):
                    for j in range(len(args[0])):
                        if args[0][i][j] != 0:
                            self.edges.append(make_edge(i, j, args[0][i][j]))
            else:            
                wrong_edges = [i for i in args[0] if len(i) not in [2, 3]]
                if wrong_edges:
                    raise ValueError('Wrong edges data: {}'.format(wrong_edges))
            
                self.edges = [make_edge(*edge) for edge in args[0]]
        else:
            self.edges = [make_edge(*edge) for edge in args[0]]

    @property
    def vertices(self):
	#Devuelve la lista de vertices del grafo
        return set(
            sum(
                ([edge.start, edge.end] for edge in self.edges), []
            )
        )

    def get_node_pairs(self, n1, n2, both_ends=True):
	#Devuelve los nodos pares, se usa en caso de que el grafo sea indirecto
        if both_ends:
            node_pairs = [[n1, n2], [n2, n1]]
        else:
            node_pairs = [[n1, n2]]
        return node_pairs

    def remove_edge(self, n1, n2, both_ends=True):
	#elimina del grafo el enlace formado por los nodos n1 y n2
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        edges = self.edges[:]
        for edge in edges:
            if [edge.start, edge.end] in node_pairs:
                self.edges.remove(edge)

    def add_edge(self, n1, n2, cost=1, both_ends=True):
	#agrega al grafo el enlace formado por n1 y n2 si este no existe
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        for edge in self.edges:
            if [edge.start, edge.end] in node_pairs:
                return ValueError('Edge {} {} already exists'.format(n1, n2))

        self.edges.append(Edge(start=n1, end=n2, cost=cost))
        if both_ends:
            self.edges.append(Edge(start=n2, end=n1, cost=cost))

    @property
    def neighbours(self):
	#Devuelve el listado de vecinos de un vertice del grafo
        neighbours = {vertex: set() for vertex in self.vertices}
        for edge in self.edges:
            neighbours[edge.start].add((edge.end, edge.cost))

        return neighbours

    def dijkstra(self, source, dest):
	#Devuelve un listado de todos los caminos posibles entre source y dest,
	#también devuelve un listado de costos por cada camino
        
        assert source in self.vertices, 'Such source node doesn\'t exist'
        distances = {vertex: inf for vertex in self.vertices}
        previous_vertices = {
            vertex: None for vertex in self.vertices
        }
        distances[source] = 0
        vertices = self.vertices.copy()
        #print "Vertices: {}".format(vertices)
        while vertices:
            current_vertex = min(
                vertices, key=lambda vertex: distances[vertex])
            vertices.remove(current_vertex)
            if distances[current_vertex] == inf:
                break
            for neighbour, cost in self.neighbours[current_vertex]:
                alternative_route = distances[current_vertex] + cost
                if alternative_route < distances[neighbour]:
                    distances[neighbour] = alternative_route
                    previous_vertices[neighbour] = current_vertex

        path, current_vertex = deque(), dest
        
        while previous_vertices[current_vertex] is not None:
            path.appendleft(current_vertex)
            current_vertex = previous_vertices[current_vertex]
            
        if path:
            path.appendleft(current_vertex)
        return list(path), distances[dest]
    
    def dijkstra_AB (self, source, dest, mat_demandas, mat_estado_capacidades, nodos_usados=None):
	#devuelve el listado de nodos del camino más corto entre source y dest considerando el AB disponibles
	#en la matriz de capacidades de enlaces y la distancia
        
        assert source in self.vertices, 'Such source node doesn\'t exist'
        distances = {vertex: inf for vertex in self.vertices}
        previous_vertices = {
            vertex: None for vertex in self.vertices
        }
        distances[source] = 0
        vertices = self.vertices.copy()
        
        while vertices:
            current_vertex = min(
                vertices, key=lambda vertex: distances[vertex])
            vertices.remove(current_vertex)
            if distances[current_vertex] == inf:
                break
            for neighbour, cost in self.neighbours[current_vertex]:
                alternative_route = distances[current_vertex] + cost
                if alternative_route < distances[neighbour] and \
                mat_estado_capacidades[current_vertex][neighbour] + mat_demandas[current_vertex][neighbour] <= globales.MATRIZ_CAPACIDADES[current_vertex][neighbour]:
                    try:
                        if neighbour not in nodos_usados:
                            distances[neighbour] = alternative_route
                            previous_vertices[neighbour] = current_vertex
                    except TypeError as e:
                        distances[neighbour] = alternative_route
                        previous_vertices[neighbour] = current_vertex
        
        path, current_vertex = deque(), dest
        
        while previous_vertices[current_vertex] is not None:
            path.appendleft(current_vertex)
            current_vertex = previous_vertices[current_vertex]
            
        if path:
            path.appendleft(current_vertex)
        
        return list(path), distances[dest]
    
    def meta_dijkstra (self, source):
	#Devuelve el camino más corto desde el origen a todos los nodos del grafo y sus distancias
        assert source in self.vertices, 'Such source node doesn\'t exist'
        distances = {vertex: inf for vertex in self.vertices}
        previous_vertices = {
            vertex: None for vertex in self.vertices
        }
        distances[source] = 0
        vertices = self.vertices.copy()
        while vertices:
            current_vertex = min(
                vertices, key=lambda vertex: distances[vertex])
            vertices.remove(current_vertex)
            if distances[current_vertex] == inf:
                break
            for neighbour, cost in self.neighbours[current_vertex]:
                alternative_route = distances[current_vertex] + cost
                if alternative_route < distances[neighbour]:
                    distances[neighbour] = alternative_route
                    previous_vertices[neighbour] = current_vertex
        return previous_vertices, distances

def elegir_ruta(destino, posibles_caminos, posibles_costos):
#esta función recibe el resultado de meta_dijkstra y elige el mejor camino para el destino
    camino, actual_vertex = deque(), destino
    
    while posibles_caminos[actual_vertex] is not None:
        camino.appendleft(actual_vertex)
        actual_vertex = posibles_caminos[actual_vertex]
        
    if camino:
        camino.appendleft(actual_vertex)
    
    return camino, posibles_costos[destino]

def calcular_ruteo_costos_sfp():
"""
Calcula la matriz spf_matrix donde se retorna el listado de todos los caminos entre todos los
nodos del grafo y sus costos correspondientes, utiliza las funciones meta_dijkstra y elegir ruta_aux
para ciclar entre todos los nodos del grafo.
"""
    grafo = Graph(globales.MATRIZ_ADYACENCIA_PONDERADA)
    vertices=list(grafo.vertices)
    mat_sfp_costos = []
    mat_sfp_ruteo = []

    for i in range(len(globales.MATRIZ_ADYACENCIA_PONDERADA)):
        mat_sfp_costos.append([])
        mat_sfp_ruteo.append([])
        posibles_caminos, posibles_costos=grafo.meta_dijkstra(vertices[i])
        for j in range (len(globales.MATRIZ_ADYACENCIA_PONDERADA[i])):
            ruta_aux, costo_aux=elegir_ruta(vertices[j],posibles_caminos, posibles_costos)
     
            mat_sfp_costos[i].append(costo_aux)
            mat_sfp_ruteo[i].append(ruta_aux)
    return mat_sfp_ruteo, mat_sfp_costos#, mat_estado_capacidades 

def evaluar_individuo(cromosoma, sfp_costos_gral, sfp_rutas_gral, demandas_gral):
    cromo_dispositivos=cromosoma[1]
    cromo_rutas=cromosoma[0]
    sdn_flow=cromosoma[2]
    matriz_aux=modificar_matriz(globales.MATRIZ_ADYACENCIA_PONDERADA, cromo_dispositivos, globales.TIPOS_DISPOSITIVOS)
    grafo_aux=Graph(matriz_aux)
    estado_capacidades_aux=[[0 for i in range(globales.CANTIDAD_DISPOSITIVOS)] for j in range(globales.CANTIDAD_DISPOSITIVOS)]
    costo_reemplazo=0
    costo_incremento=0
    costo_enlaces=0
    costo_calidad=0
    costo_calidad_incremento=0
    costo_calidad_reemplazo=0
    indice=0
    nodos_verificados=[]
    flujos_sdn=[0 for index in range(globales.CANTIDAD_SOLICITUDES)]
    cant_sdn_cromosoma=0
    for contador in cromo_dispositivos:
                cant_sdn_cromosoma=cant_sdn_cromosoma+contador
    while True:
        recalcular=False
        origen = globales.LISTA_SOLICITUDES[cromo_rutas[indice]][0]
        destino = globales.LISTA_SOLICITUDES[cromo_rutas[indice]][1]
        camino = sfp_rutas_gral[origen][destino]
        last=None
        del nodos_verificados[:]
        if sdn_flow[indice]==0:
            for nodo in camino:
                
                if cromo_dispositivos[nodo] == 0 and globales.TIPOS_DISPOSITIVOS[nodo]==0: #Si en el camino mas corto maestro hay un dispositivo nuevo pero este no es candidato para el upgrade entonces se recalcula
                    recalcular=True
                    restaurar_matriz_capacidades(origen, destino, estado_capacidades_aux,nodos_verificados,demandas_gral) #controlar que el vector sea >0
                    del nodos_verificados[:]
                    break
                
                if last is not None:
                    if estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino] < globales.MATRIZ_CAPACIDADES[nodo][last]:
                        estado_capacidades_aux[nodo][last] = estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino]
                        estado_capacidades_aux[last][nodo] = estado_capacidades_aux[last][nodo]+demandas_gral[origen][destino]
                    else:
                        recalcular = True
                        restaurar_matriz_capacidades(origen, destino, estado_capacidades_aux, nodos_verificados, demandas_gral)
                        del nodos_verificados[:]
                        break

                last = nodo
                nodos_verificados.append(nodo)
            if recalcular:
                sdn_aux=get_cheapest_sdn_aux(grafo_aux, cromo_dispositivos, origen, destino, demandas_gral,
                                                                       estado_capacidades_aux)
                aux_costo_enlaces_sdn= sdn_aux[3]
                nodos_verificados=sdn_aux[4]
                last=None
                aux_costo_enlace_sfp = grafo_aux.dijkstra_AB(origen, destino,demandas_gral,
                                                             estado_capacidades_aux)
                if aux_costo_enlaces_sdn < aux_costo_enlace_sfp[1] and aux_costo_enlaces_sdn > 0 and aux_costo_enlaces_sdn < float('inf'):

                    costo_enlaces = costo_enlaces + aux_costo_enlaces_sdn

                    for nodo in nodos_verificados:
                        if last is not None:
                            if estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino] < globales.MATRIZ_CAPACIDADES[nodo][last]:
                                estado_capacidades_aux[nodo][last] = estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino]
                                estado_capacidades_aux[last][nodo] = estado_capacidades_aux[last][nodo]+demandas_gral[origen][destino]
                            
                        last=nodo

                else:

                    if aux_costo_enlace_sfp[1] < float('inf'):
                        costo_enlaces = costo_enlaces +  aux_costo_enlace_sfp[1]
                        nodos_verificados=aux_costo_enlace_sfp[0]
                        for nodo in nodos_verificados:
                            if last is not None:
                                if estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino] < globales.MATRIZ_CAPACIDADES[nodo][last]:
                                    estado_capacidades_aux[nodo][last] = estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino]
                                    estado_capacidades_aux[last][nodo] = estado_capacidades_aux[last][nodo]+demandas_gral[origen][destino]
                                
                            last=nodo

            else:

                if sfp_costos_gral[origen][destino] < float('inf') and sfp_costos_gral[origen][destino] > 0:
                    costo_enlaces = costo_enlaces + sfp_costos_gral[origen][destino]                                                  
        else:
            cant_sdn_camino=0
            for nodo in camino:
                if cromo_dispositivos[nodo]==0 and globales.TIPOS_DISPOSITIVOS[nodo]==0:
                    recalcular=True
                    break
                    
                if last is not None:
                    if estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino] < globales.MATRIZ_CAPACIDADES[nodo][last]:
                        estado_capacidades_aux[nodo][last] = estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino]
                        estado_capacidades_aux[last][nodo] = estado_capacidades_aux[last][nodo]+demandas_gral[origen][destino]
                    else:
                        recalcular = True
                        break

                cant_sdn_camino=cant_sdn_camino+cromo_dispositivos[nodo]
                last = nodo
                nodos_verificados.append(nodo)
                
            if cant_sdn_cromosoma > 0 and (recalcular or cant_sdn_camino==0):
                restaurar_matriz_capacidades(origen, destino, estado_capacidades_aux, nodos_verificados, demandas_gral)
                del nodos_verificados[:]
                last=None
                sdn_aux=get_cheapest_sdn_aux(grafo_aux, cromo_dispositivos, origen, destino, demandas_gral,
                                                                        estado_capacidades_aux)
                if sdn_aux[3] < float('inf') and sdn_aux[3] > 0:
                    costo_enlaces=costo_enlaces+sdn_aux[3]
                    nodos_verificados=sdn_aux[4]
					
                    for nodo in nodos_verificados:
                        if last is not None:
                            if estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino] < globales.MATRIZ_CAPACIDADES[nodo][last]:
                                estado_capacidades_aux[nodo][last] = estado_capacidades_aux[nodo][last]+demandas_gral[origen][destino]
                                estado_capacidades_aux[last][nodo] = estado_capacidades_aux[last][nodo]+demandas_gral[origen][destino]
                                
                        last=nodo
                
            else:

                if sfp_costos_gral[origen][destino] < float('inf') and sfp_costos_gral[origen][destino] > 0:
                    costo_enlaces = costo_enlaces + sfp_costos_gral[origen][destino]

        
        for nodo in range(len(nodos_verificados)):
            if cromo_dispositivos[nodos_verificados[nodo]] > 0:
                costo_calidad+=1
                for nodo_i in nodos_verificados[nodo:]:
                    if globales.TIPOS_DISPOSITIVOS[nodo_i] == 0:
                        costo_calidad_incremento +=1
                        break
                break

        indice+=1
        if indice >= len(cromo_rutas):
            break
    for i in range(len(cromo_dispositivos)):
        costo_reemplazo = costo_reemplazo + cromo_dispositivos[i]*globales.CUD[i]*globales.TIPOS_DISPOSITIVOS[i]

    costo_calidad_reemplazo = costo_calidad - costo_calidad_incremento
    if costo_calidad_reemplazo < 0:
        print 'Error'

    fitness= globales.ALFA_1*(((costo_reemplazo+costo_incremento)-globales.MIN_CD)/(globales.MAX_CD-globales.MIN_CD)) - globales.ALFA_2*((costo_calidad-globales.MIN_Q)/(globales.MAX_Q-globales.MIN_Q))+globales.ALFA_3*((costo_enlaces-globales.MIN_CL)/(globales.MAX_CL-globales.MIN_CL))

    
         
    return fitness, costo_reemplazo, costo_incremento, costo_calidad, costo_calidad_reemplazo, costo_calidad_incremento, costo_enlaces

def cerar_matriz(matriz):
#Es una función auxiliar para copiar una matriz ya que las copias se hacen por referencias nada más
    dimensiones=cantidad_dimensiones(matriz)
    if dimensiones > 1:
        for indice in matriz:
            cerar_matriz(indice)
    else:
        for indice in matriz:
            indice=0

def restaurar_matriz_capacidades(origen, destino, matriz, camino, mat_demandas):
#Elimina la moficación realizada por la última operación o la última solicitud
    last=None
    for nodo in camino:
        if last is not None:
            if matriz[nodo][last]-mat_demandas[origen][destino] < 0:
                print 'Whaaat???'
            matriz[nodo][last]=matriz[nodo][last]-mat_demandas[origen][destino]
            matriz[last][nodo]=matriz[last][nodo]-mat_demandas[origen][destino]
        last=nodo

def evaluar_poblacion(poblacion, sfp_costos_gral, sfp_rutas_gral, demandas_gral, aptitudes=None , mejor=None, aux=None):
#Engloba a la función de evaluación individual para poder ciclar a través de la población´, devuelve la aptitud de toda la población,
#así como los valores que forman parte de dicha aptitud para su posterior análisis
    fitness=[]
    cr=[]
    ci=[]
    cq=[]
    cqr=[]
    cqi=[]
    cl=[]
    contador=0
    indice=0
    total= len(poblacion[0])
    individuo=[]
    if mejor != None:
        fitness.append(aptitudes[mejor])
        cr.append(aux[1][mejor])
        ci.append(aux[2][mejor])
        cq.append(aux[3][mejor])
        cqr.append(aux[4][mejor])
        cqi.append(aux[5][mejor])
        cl.append(aux[6][mejor])
        contador+=1
        indice+=1
    while contador < total:
        individuo.append(poblacion[0][indice])
        individuo.append(poblacion[1][indice])
        individuo.append(poblacion[2][indice])
        aux_resul=evaluar_individuo(individuo, sfp_costos_gral, sfp_rutas_gral, demandas_gral)
        fitness.append(aux_resul[0])
        cr.append(aux_resul[1])
        ci.append(aux_resul[2])
        cq.append(aux_resul[3])
        cqr.append(aux_resul[4])
        cqi.append(aux_resul[5])
        cl.append(aux_resul[6])
        del individuo[:]
        contador +=1
        indice+=1
            
    return fitness, cr, ci, cq, cqr, cqi, cl

def modificar_matriz(mat_adya, cromo_dispositivos, dispositivos):
"""
Esta función modifica la matriz de adyacencia de acuerdo a lo que represente un cromosoma específico
Esto se hace ya que la spf_matrix posee todos los caminos más cortos posibles con sus respectivos costos,
pero puede ocurrir que un camino inicialmente posible ya no lo sea porque algún nodo para agregación finalmente
no se considere para ese individuo específico
"""
    aux = []
    for i in range(len(dispositivos)):
        aux.append([])
        for j in range(len(dispositivos)):
            if cromo_dispositivos[i] == 0 and dispositivos[i] == 0 or cromo_dispositivos[j] == 0 and dispositivos[j] == 0:
                aux[i].append(0)
            else:
                aux[i].append(mat_adya[i][j])
    return aux


def get_cheapest_sdn_aux(grafo, cromo_dispositivos, origen, destino, demandas_gral, estado_capacidades_aux):
'''
El metodo get_cheapest_sdn_aux busca algún dispositivo que se proponga como sdn para enviar el tráfico de una solicitud
si no encuentra ninguno el resultado es infinito
'''
    cheapest = float('inf')
    costo=float('inf')
    cheapest_dispositivo=0
    cheapest_camino=[]
    auxiliar=grafo.meta_dijkstra(origen)
    auxiliar=ordenar_costos_ruteo(auxiliar[0].items(), auxiliar[1].items())
    index=0
    nodos_verificados=[]
    for i,j in auxiliar[1]:
        del nodos_verificados[:]
        previo=auxiliar[0][index][1]
        if cromo_dispositivos[i]>0 and i != origen and estado_capacidades_aux[i][previo] + demandas_gral[i][previo] <= globales.MATRIZ_CAPACIDADES[i][previo]:
            cheapest=j
            cheapest_dispositivo=i
            cheapest_camino=elegir_ruta_AB(i, auxiliar[0], auxiliar[1]).__getitem__(0)
            for nodo in cheapest_camino:
                nodos_verificados.append(nodo)
            sdn_aux = grafo.dijkstra_AB(i, destino, demandas_gral, estado_capacidades_aux, nodos_verificados)
            if sdn_aux[1]!=float('inf'):
                for nodo in sdn_aux[0][1:]:
                    nodos_verificados.append(nodo)
                break
        index+=1
       
    if cheapest!=float('inf'):
        costo=cheapest+sdn_aux[1]
    
    return cheapest, cheapest_dispositivo, cheapest_camino, costo, nodos_verificados

def ordenar_costos_ruteo(rutas, costos):
#Se utiliza para ordenar las rutas según el costo de ellas. a fin de seleccionar 
#la más cercana
    if len(costos) <= 1:
        return rutas, costos
    medio = len(costos)//2
    R_izq, C_izq= ordenar_costos_ruteo(rutas[:medio], costos[:medio])
    R_der, C_der=ordenar_costos_ruteo(rutas[medio:], costos[medio:])
    copia_rutas=rutas[:]
    copia_costos=costos[:]
    return merge(R_izq, C_izq, R_der, C_der, copia_rutas, copia_costos)

def merge(R_izq, C_izq, R_der, C_der, rutas_mezcladas, costos_mezclados):
    cursor_izq, cursor_der=0, 0
    while cursor_izq < len(C_izq) and cursor_der < len(C_der):
        if C_izq[cursor_izq][1] <= C_der[cursor_der][1]:
            rutas_mezcladas[cursor_izq+cursor_der]=R_izq[cursor_izq]
            costos_mezclados[cursor_izq+cursor_der]=C_izq[cursor_izq]
            cursor_izq+=1
        else:
            rutas_mezcladas[cursor_izq+cursor_der]=R_der[cursor_der]
            costos_mezclados[cursor_izq+cursor_der]=C_der[cursor_der]
            cursor_der+=1
            
    for cursor_izq in range(cursor_izq, len(C_izq)):
        rutas_mezcladas[cursor_izq+cursor_der]= R_izq[cursor_izq]
        costos_mezclados[cursor_izq+cursor_der]= C_izq[cursor_izq]
    
    for cursor_der in range(cursor_der, len(C_der)):
        rutas_mezcladas[cursor_izq+cursor_der]= R_der[cursor_der]
        costos_mezclados[cursor_izq+cursor_der]= C_der[cursor_der]
    
    return rutas_mezcladas, costos_mezclados


def elegir_ruta_AB(destino, posibles_caminos, posibles_costos):
#Elige las rutas de acuerdo a la capacidad de los enlaces disponibles
    camino, actual_vertex = deque(), destino
    diccionario=dict(posibles_caminos)
    diccionario_costos=dict(posibles_costos)
    while diccionario[actual_vertex] is not None:
        camino.appendleft(actual_vertex)
        actual_vertex = diccionario[actual_vertex]
        
    if camino:
        camino.appendleft(actual_vertex)
        
    try:
        if diccionario_costos[destino] ==float('inf'):
            print ''
    except IndexError as e:
        print e
    
    return camino, diccionario_costos[destino]

def verificacion_preliminar():
#Realiza las comprobaciones preliminares del algoritmo, en este caso solo una verificación hasta ahora
#que es la verificacióin de rutas
    if not verificacion_rutas():
        return False
    return True

def verificacion_rutas():
#Verfica que en las solicitudes todos los nodos existan inicialmente
    for i in globales.LISTA_SOLICITUDES:
        if globales.TIPOS_DISPOSITIVOS[i[0]] == 0 or globales.TIPOS_DISPOSITIVOS[i[1]] == 0:
            return False
    return True

def copiar_poblacion(datos):
#copia la poblacion en cada generación para actualizar los datos de acuerdo a los operadores evolutivos.
    aux=[]
    for i in range(len(datos)):
        aux.append([])
        for j in range(len(datos[i])):
            aux[i].append([])
            for k in range(len(datos[i][j])):
                aux[i][j].append(datos[i][j][k])
    return aux

def construir_mat_demandas():
#construye una matriz de demandas de acuerdo al listado de solicitudes inicial
    mat_demandas=[[0 for i in range(globales.CANTIDAD_DISPOSITIVOS)] for j in range(globales.CANTIDAD_DISPOSITIVOS)]
    for i in globales.LISTA_SOLICITUDES:
        mat_demandas[i[0]][i[1]]=i[2]
        mat_demandas[i[1]][i[0]]=i[2]
    return mat_demandas

def duplicar_matriz_2d(matriz):
    duplicado=[]
    for i in range(len(matriz)):
        duplicado.append([])
        for j in range(len(matriz[i])):
            duplicado[i].append(matriz[i][j])
    return duplicado
