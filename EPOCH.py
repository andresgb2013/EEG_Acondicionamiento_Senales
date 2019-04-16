# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 22:36:22 2019

@author: Andres Gonzalez
"""

import numpy as np;
import matplotlib.pyplot as plt;
import scipy.signal as signal
import numpy as np
import math as mt
from scipy.stats import kurtosis

#%%
def datos(nombre):
    dato=[];
    senal=[]
    archivo = open(nombre,'r')  
    lineas = (archivo.readlines())
    linea = lineas[6:]
    for i in range(0,len(linea)):
        r = linea[i].split(', ')
        dato.append(list(map(float,r[1:9])))
        for j in range(0,8):
            senal.append(dato[i][j]);
    return senal

#%%
def Filtered_Signal(Data,Colums,Fsenal):
    Signal_Fil=[]
    for i in np.arange(0,Colums):
        Signal_Fil =eegfiltnew(Data, Fsenal, 50 , 1 , 0, 0);

    return(Signal_Fil)
#%%
def eegfiltnew(senal, srate, locutoff = 0, hicutoff = 0, revfilt = 0, plot = 0):
    #Constants
    TRANSWIDTHRATIO = 0.25;
    fNyquist = srate/2;   
    
    if hicutoff == 0: #Convert highpass to inverted lowpass
        hicutoff = locutoff
        locutoff = 0
        revfilt = 1 #este valor se cambia para notch y tambien se debe cambiar en este caso
    if locutoff > 0 and hicutoff > 0:
        edgeArray = np.array([locutoff , hicutoff])
    else:
        edgeArray = np.array([hicutoff]);
    
    if np.any(edgeArray<0) or np.any(edgeArray >= fNyquist):
        return False    
    
    # Max stop-band width
    maxTBWArray = edgeArray.copy() # Band-/highpass
    if revfilt == 0: # Band-/lowpass
        maxTBWArray[-1] = fNyquist - edgeArray[-1];
    elif len(edgeArray) == 2: # Bandstop
        maxTBWArray = np.diff(edgeArray) / 2;
    maxDf = np.min(maxTBWArray);
    
    # Default filter order heuristic
    if revfilt == 1: # Highpass and bandstop
        df = np.min([np.max([maxDf * TRANSWIDTHRATIO, 2]) , maxDf]);
    else: # Lowpass and bandpass
        df = np.min([np.max([edgeArray[0] * TRANSWIDTHRATIO, 2]) , maxDf]);
        
    filtorder = 3.3 / (df / srate); # Hamming window
    filtorder = np.ceil(filtorder / 2) * 2; # Filter order must be even.
    
    # Passband edge to cutoff (transition band center; -6 dB)
    dfArray = [[df, [-df, df]] , [-df, [df, -df]]];
    cutoffArray = edgeArray + np.array(dfArray[revfilt][len(edgeArray) - 1]) / 2;
    # Window
    winArray = signal.hamming(int(filtorder) + 1);
    # Filter coefficients
    if revfilt == 1:
        filterTypeArray = ['high', 'stop'];
        b = firws(filtorder, cutoffArray / fNyquist, winArray, filterTypeArray[len(edgeArray) - 1]);
    else:
        b = firws(filtorder, cutoffArray / fNyquist, winArray);
    
    if plot == 1:
        #plot the response of the filter
        mfreqz(b,1,filtorder, fNyquist);
    
    signal_filtered = signal.filtfilt(b, 1, senal);#
    return signal_filtered;    
#%%
def firws(m, f , w , t = None):    
    try:
        m = int(m)
        if  (m%2 != 0) or (m<2):
            return False
    except (ValueError, TypeError):
        return False
    
    if (not type(f) is np.ndarray):
        return False
    
    f = np.squeeze(f)
    if (f.ndim > 1) or (f.size > 2):
        return False
    f = f / 2; 
    
    if np.any(f <= 0) or np.any(f >= 0.5):
        return False
    
    w = np.squeeze(w)
    if (f.ndim == 0):
        b = fkernel(m, f, w)
    else:
        b = fkernel(m, f[0], w)
    if (f.ndim == 0) and (t == 'high'):
        b = fspecinv(b)
    elif (f.size == 2):
        b = b + fspecinv(fkernel(m, f[1], w))
        if t == None or (t != 'stop'):
            b = fspecinv(b)        
    return b    
#%%
def fkernel(m, f, w):
    m = np.arange(-m/2, (m/2)+1)
    b = np.zeros((m.shape[0]))
    b[m==0] = 2*np.pi*f # No division by zero
    b[m!=0] = np.sin(2*np.pi*f*m[m!=0]) / m[m!=0] # Sinc
    b = b * w # Window
    b = b / np.sum(b) # Normalization to unity gain at DC
    return b  
#%%
def fspecinv(b):
    b = -b
    b[int((b.shape[0]-1)/2)] = b[int((b.shape[0]-1)/2)]+1
    return b    
#%%
def canales(senal):
    canal=[];
    for i in range(0,8):
        canal.append(senal[i:len(senal):8])
    return canal

def filtrados(canal):
    filtrado=[];
    for i in range(0,8):
        filtrado.append(Filtered_Signal(canal[i],1,250))
    return filtrado
#%%
def epocas(s,canal,pm):
    # Generacion de epocas
    time=len(canal);
    tepc= int((time/pm)/s);
    n=pm*s;
    umbralpos=np.zeros(n);
    umbralneg=np.zeros(n);
    epocas=np.zeros((n,int(len(canal)/n)));
    for i in range(0,tepc):
        vector = canal[int(i*n):int((i+1)*n)]
        
        medias= np.mean(vector);
        desvestd=np.std(vector);
        umbralpos[i] = (medias+(3.5)*desvestd);
        umbralneg[i] = (medias-(3.5)*desvestd);
        epocas[:,i] = vector;

    #Metodo de Umbral desde acá
    soluc=np.ones(tepc);
    for j in range(tepc):
        atipicx=epocas[:,j];
        comp1=atipicx[atipicx>umbralpos[j]];
        comp2=atipicx[atipicx<umbralneg[j]];
        atipicx=epocas[j,:];
        if comp1.size>0 or comp2.size>0:
            soluc[j]=1;
        elif comp1.size==0 and comp2.size==0:
            soluc[j]=0;
            
    return soluc,tepc,epocas;

#%% Tercer método, aplicación de Curtosis
def curtosis(epocas,tepc):
    soluc=np.ones(tepc);
    a=[]
    for j in range(tepc):
        c=int(kurtosis(epocas[:,j]));
        a.append(c);
        if a[j]>0 :
            soluc[j]= 1;
        else:
            soluc[j]=0          
    return soluc;
#%%
def tendencia(epocas,tepc):
    soluc=np.ones(tepc);
    w=[]
    for j in range(tepc):
        c= signal.detrend(epocas[:,j]);
        w.append(c);
        t=np.arange(0,tepc);     
        m=(c[10]-c[1])/(t[1]-t[10])
        if 1>m:
             soluc[j]=0;
        else:
             soluc[j]=1;
        
    return soluc;
#%%
def potencia(epocas,tepc):
    soluc=np.ones(tepc);
    a=[]
    for j in range(tepc):
        [e1,e2]=signal.welch(epocas[:,j],250)
        a.append(e2)
        if max(e2)>= 15 :
            soluc[j]= 1;
            
        else:
            soluc[j]=0          
    return soluc;
#%%

def eliminar(epocas, atipica):
    new=[];
    for j in range(0,len(atipica)):
        for i in range(0,len(epocas[:,0])):
            if atipica[j]==0:
                new.append(epocas[i,j]);
    return new