# -*- coding: utf-8 -*-
'''
Created on 2 jul. 2019

@author: Bader Maluff
'''

if __name__ == '__main__':
    pass

import random
import AG_funciones as ag_gene
import globales #esta es un archivo que lee el json con los datos del modelo y se importa como un modulo para poder tener todos los datos de forma global
import sys
import time
tiempo_inicio=time.clock()
if not ag_gene.verificacion_preliminar():
    print "Por favor ingrese los datos correctos"   
    sys.exit()
mat_demandas_gral=ag_gene.construir_mat_demandas()
sfp_rutas_gral, sfp_costos_gral=ag_gene.calcular_ruteo_costos_sfp()
poblacion=list(ag_gene._generar_poblacion(globales.CANTIDAD_INDIVIDUOS, globales.CANTIDAD_DISPOSITIVOS, globales.CANTIDAD_SOLICITUDES))
aux_resul = ag_gene.evaluar_poblacion(poblacion, sfp_costos_gral, sfp_rutas_gral, mat_demandas_gral)
aptitudes=aux_resul[0]
for i in range(globales.CANTIDAD_INDIVIDUOS):
    print '{0};{1};{2};{3};CR: {4};CI:{5};CQ: {6};CQR: {7};CQI: {8};CL: {9}'.format(aptitudes[i], poblacion[0][i], poblacion[1][i], poblacion[2][i], aux_resul[1][i], aux_resul[2][i], aux_resul[3][i], aux_resul[4][i], aux_resul[5][i], aux_resul[6][i])
nueva_generacion=ag_gene.copiar_poblacion(poblacion)
for i in range(50):
    nueva_generacion = ag_gene.sel_x_mut(nueva_generacion, aptitudes, globales.PORCENTAJE_MUTACION, globales.RUIDO_PERMUTACION, globales.RUIDO_FLIPBIT)
    aux_resul=ag_gene.evaluar_poblacion(nueva_generacion, sfp_costos_gral, sfp_rutas_gral, mat_demandas_gral)
    aptitudes = aux_resul[0]
    print "Generacion {}".format(i+1)
    
    for i in range(globales.CANTIDAD_INDIVIDUOS):
        print '{0};{1};{2};{3};CR: {4};CI:{5};CQ: {6};CQR: {7};CQI: {8};CL: {9}'.format(aptitudes[i], nueva_generacion[0][i], nueva_generacion[1][i], nueva_generacion[2][i], aux_resul[1][i], aux_resul[2][i], aux_resul[3][i], aux_resul[4][i], aux_resul[5][i], aux_resul[6][i])
tiempo_fin=time.clock()
print 'Tiempo de ejecución total: {0}'.format(tiempo_fin-tiempo_inicio)
print "El mejor individuo es: "
i=0
print '{0};{1};{2};{3};CR: {4};CI:{5};CQ: {6};CQR: {7};CQI: {8};CL: {9}'.format(aptitudes[i], nueva_generacion[0][i], nueva_generacion[1][i], nueva_generacion[2][i], aux_resul[1][i], aux_resul[2][i], aux_resul[3][i], aux_resul[4][i], aux_resul[5][i], aux_resul[6][i])
file_log=open("resultados_maestro.txt","a")
for i in range(globales.CANTIDAD_INDIVIDUOS):
    file_log.write("Apt: "+str(aptitudes[i])+"; Ord_Flow: "+str(nueva_generacion[0][i])+";SDN_Flow: "+str(nueva_generacion[1][i])+";SDN_Dev: "+str(nueva_generacion[2][i])+"; CR: "+str(aux_resul[1][i])+"; CI: "+str(aux_resul[2][i])+"; CQ: "+str(aux_resul[3][i])+"; CQR: "+str(aux_resul[4][i])+"; CQI: "+str(aux_resul[5][i])+"; CL: "+str(aux_resul[6][i])+"\n")
file_log.write("Tiempo de Ejecución: "+str(tiempo_fin-tiempo_inicio)+'\n')
file_log.close()
