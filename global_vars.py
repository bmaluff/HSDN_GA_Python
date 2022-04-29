'''
Created on 6 ago. 2019

@author: Bader Maluff
'''
import json

with open('datos.json') as json_file:
    json_object=json.load(json_file)
CANTIDAD_INDIVIDUOS=json_object['tamano_poblacion']
CANTIDAD_DISPOSITIVOS=json_object['cantidad_dispositivos']
CANTIDAD_SOLICITUDES=json_object['cantidad_solicitudes']
MATRIZ_ADYACENCIA_PONDERADA=json_object['matriz_adyacencia_ponderada']
TIPOS_DISPOSITIVOS=json_object['tipos_dispositivos']
CUD=json_object['costo_unificado_dispositivos']
MIN_CD=json_object['costo_dispositivo_min']
MAX_CD=json_object['costo_dispositivo_max']
MIN_CL=json_object['costo_enlaces_min']
MAX_CL=json_object['costo_enlaces_max']
MIN_Q=json_object['costo_calidad_min']
MAX_Q=json_object['costo_calidad_max']
LISTA_SOLICITUDES=json_object['solicitudes']
ALFA_1=json_object['alfa1']
ALFA_2=json_object['alfa2']
ALFA_3=json_object['alfa3']
PORCENTAJE_MUTACION = json_object['porcentaje_mutacion']
RUIDO_PERMUTACION=json_object['ruido_permutacion']
RUIDO_FLIPBIT=json_object['ruido_flipbit']
MATRIZ_CAPACIDADES = json_object['matriz_capacidades']
