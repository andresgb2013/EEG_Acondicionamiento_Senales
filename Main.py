# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:31:59 2019

@author: Juan Diego
"""
import matplotlib.pyplot as plt;
import numpy as np
import EPOCH ;
#%% Menu de selección de archivo
a=int(input('Seleccione el archivo a cargar : \n 1. P1_RAWEEG_2018-11-15_Electrobisturí1_3min \n 2. P1_RAWEEG_2018-11-15_Electrobisturí2_2min \n 3. P1_RAWEEG_2018-11-15_FinProcedimiento_53min \n 4. P1_RAWEEG_2018-11-15_OjosCerrados_2min] :\n ') ) 
if a==1:
    archivo= 'P1_RAWEEG_2018-11-15_Electrobisturí1_3min.txt'
elif a==2:
    archivo ='P1_RAWEEG_2018-11-15_Electrobisturí2_2min.txt'
elif a==3:
    archivo ='P1_RAWEEG_2018-11-15_FinProcedimiento_53min.txt'
elif a==4:
    archivo ='P1_RAWEEG_2018-11-15_OjosCerrados_2min.txt'
#%% Menu de selección de metodología
senal = EPOCH.datos(archivo);
canales = EPOCH.canales(senal);
s= int(input('Ingrese el numero de segundo por época :'))
c= int(input('Ingrese el canal de interés [0-7] :'))
b= int(input('Escoja una metodología de preprocesamiento : \n 1. Filtrado, División por épocas, Detección y eliminacuión de atípicas por curtosis \n 2. Filtrado, División por épocas, Detección y eliminacuión de atípicas por potencias \n 3. Filtrado, división por épocas, detección de épocas etípicas por métodos de umbral'))
if b==1:
    Data_Filtered = EPOCH.filtrados(canales); #Filtrado
    [soluc,tepc,epocas]= EPOCH.epocas(s,Data_Filtered[c],250);
    soluc= EPOCH.curtosis(epocas,tepc);
    corregidas= EPOCH.eliminar(epocas, soluc)
    plt.title('Metodología 1, aplicación de curtosis ');
    plt.plot(corregidas);
    plt.figure()
        
elif b==2:
    Data_Filtered = EPOCH.filtrados(canales);#Filtrado
    [soluc,tepc,epocas]= EPOCH.epocas(s,Data_Filtered[c],250);
    soluc2= EPOCH.potencia(epocas,tepc);
    corregidas2= EPOCH.eliminar(epocas, soluc2);
    plt.title('Metodología 2, aplicación de Potencias - Welch ');
    plt.plot(corregidas2)
    plt.figure()
elif b==3:
    Data_Filtered = EPOCH.filtrados(canales); #Filtrado
    [soluc3,tepc,epocas]= EPOCH.epocas(s,Data_Filtered[c],250);
    corregidas3= EPOCH.eliminar(epocas, soluc3)
    plt.title('Metodología 1, aplicación de Umbral ');
    plt.plot(corregidas3)
    plt.figure()
