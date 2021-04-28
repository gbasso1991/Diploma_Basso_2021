# -*- coding: utf-8 -*-
"""
Turbo Lover - Giuliano Andrés Basso - 23 Feb 2021 

Este script permite obtener los ciclos de magnetizacion, coercitividad, 
remanencia y el SAR (Specfic Absortion Rate) de las muestas medidas con ESAR a distintas
temperaturas.

Ahora tengo que tomar solo 1 archivo "_fondo" y "_cal"

Toma un archivo con formato de nombre 'xxxkHz_yyA_zzz_Mss_TM_Ti.dat'
Siendo:
  xxx  = frecuencia en kHz del campo.
  yy   = valor de la corriente Idc (Internal Direct Current) en Ampere. 
  zzz  = valor de muestreo (i.e., cuantos puntos por segundo registradas), en 
  megamuestras por segundo.
  TM   = Tipo de Muestra (FF ferrofluido, FG ferrogel, TT tejido tumoral, TC tejido congelado, FC = ferrofluido congelado)
  Ti = Temperatura de la medicion i-ésima.
  
> El valor de frecuencia solo es utilizado como semilla para ajustar la frecuencia 
real, y para verificar que no haya ocurrido un error al nombrar los archivos.

> El valor de corriente es empleado para verificar que el campo obtenido a partir de 
la amplitud sea coherente.

> El valor de muestreo es usado para la base temporal.

> La concentracion de NP de la muestra (g/m^3) es desconocida a priori, se obtiene 
posteriormente en una medida destructiva.

> Cada archivo .dat especifica la unidad de medida (mV en este caso), y separa la 
informacion en 3 columnas: 
                            1º Frecuency: base temporal (00001:05000) 
                            2º CH1: muestra (notar orden de magnitud ~100 mV)
                            3º CH2 = campo (~10.000 mV)

> El programa busca los 2 arhivos: 
    'xxxkHz_yyydA_zzzMss_TM_fondo.dat', con los datos del fondo 
    'xxxkHx_yyydA_zzzMss_TM_cal.dat' , con los datos de la calibración.
    
    *(La extension puede modificarse a .txt segun sea necesario)
    
> Se ajustan funciones senoidales a las tres señales de referencia.
    Se elige el tiempo cero para que las tres señales comiencen en fase.
    Tambien se dilata el tiempo de señal y referencia de fodo multiplicandolo por el 
    cociente de frecuencias, para el caso de que haya alguna pequeña 
    discrepancia, tanto para la muestra como para la calibración  (Gd en este caso).

> Se resta señal de muestra y señal de fondo, asi como tambien se resta señal de 
    calibración  y señal de fondo correspondiente. 
    A partir de este punto se trabaja unicamente con estas restas y con las referencias.   

> Se filtra el ruido aislado de cada resta, discriminando puntos donde la derivada 
    (o menos la derivada) es alta en comparacion al resto de la señal y sus entornos. 
    En dichas regiones se ajusta un polinomio con los puntos sin ruido a ambos lados de 
    la zona ruidosa.
    Se hace lo propio en las mismas regiones temporales que su señal para las respectivas 
    referencias. 

> Se recortan las señales para tener un numero entero de periodos, y se omiten el 
    primer y ultimo semi-periodo.
    De esta manera se evita el ruido que no pudo ser tratado adecuadamente al pcipio. y 
    al final de cada medida. 

> Se promedian resta y referencia sobre todos los periodos y se integra. 

> Se lleva el campo a  unidades de A/m normalizando y multiplicando por el campo maximo
    medido para el valor correspondiente de Idc.

> Se grafica el ciclo en unidades de V*s para la calibración  y se ajusta una recta.
    Con la pendiente obtenida se lleva el eje de magnizacion a unidades de A/m para 
    la muestra.

> Se guarda un archivo con la imagen del grafico del ciclo de histéresis, y otro 
    con los puntos del ciclo en ASCII.

> Se calculan coercitividad y la remanencia. 

> Se calcula el valor de SAR integrando los ciclos.

> Se imprime en un archivo de salida la siguiente informacion:
            Nombre del archivo
            Frecuencia (kHz)
            Campo maximo (kA/m)
            Peor quita de ruido porcentual.



Parte I: inputs, seleccion de archivos, chequeo errores de nomenclatura.
"""
import numpy as np
import matplotlib.pyplot as plt 
import scipy as sc              
from scipy.signal import find_peaks 
#hay mas pkg importados mas abajo

'''Masa de NP sobre volumen de FF en g/l'''
concentracion_default = 10.0 #Por si el input esta vacio

'''¿Quiero quitarle ruido a la muestra? 2=Fourier 1=Actis 0=No'''
FILTRARMUESTRA = 0 

'''Permeabilidad magnetica del vacio en N/A^2'''
mu_0 = 4*np.pi*10**-7

'''Texto que identifica los archivos de fondo'''
textofondo = '_fondo.txt' 

'''Texto que identifica los archivos de calibración '''
textocalibracion = '_cal.txt'

'''Texto que identifica los archivos del gadolinio'''
textogadolinio = '_gadolinio.dat'

'''
Importante (razon de la nomenclatura): 
name = 'xxxkHz_yydA_zzz_Mss_TM_Ti.txt'
        
    > Caracteres [7:9] del nombre de archivo deben ser el valor numerico de Idc en Ampere.
    > Caracteres [16:19] los valores de muestreo.
    > Caracteres [20:22] tipo de muestra.
    > Caracteres [23:25] temperatura de la medicion i-ésima.

Calibración  de la bobina: constante que dimensionaliza al campo en A/m a partir de la 
calibración  realizada sobre la bobina del RF
'''
pendiente_HvsI = 43.18*79.77 
ordenada_HvsI = 2.73*79.77  

'''Suceptibilidad del patron de calibración '''

rho_bulk_Gd2O3 = 7.41e3   # kg/m^3
rho_patron_Gd2O3 = 2e3   # kg/m^3
xi_bulk_Gd2O3_masa = (1.35e-4)*4*np.pi*1e-3  #emu*m/g/A = m^3/kg
xi_patron_vol = xi_bulk_Gd2O3_masa*rho_patron_Gd2O3

#%% Apertura de archivos - Cuadro de dialogo para seleccionar archivos:

import tkinter as tk
from tkinter import filedialog
root = tk.Tk()
root.withdraw()
texto_encabezado = "Seleccionar archivos con las medidas de la muestra:"
filepaths = tk.filedialog.askopenfilenames(title=texto_encabezado,filetypes=(("Archivos .txt","*.txt"),("Archivos .dat","*.dat"),("Todos los archivos","*.*")))
#cada elemento de filepahts (tupla) es un str con el path de cada archivo 

'''Nombre para el archivo ASCII de salida''' 
nombre_salida = str(input('Elegir un nombre para el archivo de salida:') or 'Prueba')

#%%
nombres_archivos=[]      #lista de str para sacar info 

for item in filepaths:    
    nombres_archivos.append(item.split('/')[-1])

print('Archivos seleccionados: ')
for item in nombres_archivos:
    print(item)
#%%Para ingresar concentracion por teclado
concentracion = float(input('Concentración de la muestra (g/l): ') or concentracion_default)

#%% Para ingresar Temperaturas por teclado    
Temperaturas = []
for i in range(0,len(nombres_archivos)):
    Temperaturas.append(float(input('Ingresar temperatura de la muestra %i: ' %i) or 20.0))

#%% Para medir el tiempo de procesamiento
import time
start_time = time.time() 

#%% =============================================================================
# Listas que uso en el procesamiento
# =============================================================================
filenombres_muestra = [] #Nombres archivos muestra
filenombres_fondo = []   #Nombres archivos de fondo
filenombres_cal = []     #Nombres archivos calibración 

for item in nombres_archivos:
    filenombres_muestra.append(item)
    filenombres_fondo.append(item[0:20]+textofondo)
    filenombres_cal.append(item[0:20]+textocalibracion)

#%%
fa=len(filepaths)
rutas_de_carga = []      #c/ elem es una lista con el path=[m,f,c]  
for i in range(0,fa):    #tantos elementos como archivos seleccione 
    rutas_de_carga.append([filepaths[i],filepaths[i],filepaths[i]])
    rutas_de_carga[i][1]=rutas_de_carga[i][1].replace(filenombres_muestra[i],filenombres_fondo[i])
    rutas_de_carga[i][2]=rutas_de_carga[i][2].replace(filenombres_muestra[i],filenombres_cal[i])


#%%
''' 
Parte II: Procesamiento de datos
    Identifica a las dos señales de cada canal como señal y referencia, para muestra,
    fondo y calibración .
    Recibe datos en mV y acá pasan a V.
    La 1er columna se pasa a tiempo.

Glosario de variables:
                        t: tiempo       m: muestra      
                        v: voltaje      f: fondo
                        r: referencia   c: calibración

Definicion funcion de suavizado: fftsmooth()'''
def fft_smooth(data_v, freq_n):
    """fft low pass filter para suavizar la señal"""
    fft_data_v = sc.fft(data_v)
    s_fft_data_v = np.zeros(len(data_v),dtype=complex)
    s_fft_data_v[0:int(freq_n)] = fft_data_v[0:int(freq_n)]
    s_fft_data_v[-1 - int(freq_n): ] = fft_data_v[-1 - int(freq_n):] 
    s_data_v = np.real(sc.ifft(s_fft_data_v))
    
    return(s_data_v)
     
'''Defino funcion a ajustar offset(A), amplitud(B), frecuencia(C), fase(d)'''
def sinusoide(t,A,B,C,D):
    return(A + B*np.sin(2*np.pi*C*t - D))

#%% =============================================================================
# Defino las salidas que obtengo del script
# =============================================================================
#0:Corriente del resonador:
Corriente_A=[]
for item in nombres_archivos:
    Corriente_A.append(float(item[7:9]))
    
#Tipo de muestra:
tipo_muestra=[]    
for item in nombres_archivos:
    if item[18:20]=='FF':
        tipo_muestra.append('Ferrofluido')
    elif item[18:20]=='FG':
        tipo_muestra.append('Ferrogel')
    elif item[18:20]=='T':
        tipo_muestra.append('Tejido tumoral')
    elif item[18:20]=='TC':
        tipo_muestra.append('Tejido congelado')
    elif item[18:20]=='FC':
        tipo_muestra.append('Ferrofluido congelado')
    else:
        tipo_muestra.append('No especificado')
        
#1: Frecuencia de la referencia en la medida de la muestra
Frecuencia_muestra_kHz=[]
    
#2: Frecuencia de la referencia en la medida del fondo
Frecuencia_fondo_kHz=[]
    
#3: Specific Absorption Rate
SAR=[]
    
#4: Campo maximo en kA/m
Campo_maximo_kAm=[]
    
#5: Campo coercitivo en kA/m
Coercitividad_kAm=[]
    
#6: Magnetizacion remanente en kA/m
Magnetizacion_remanente_kAm=[]
#Peor quita de ruido porcentual
peor_dif=[]

Mag_a_Hmax = []


#%% =============================================================================
# Defino listas para almacenar multiples ciclos y comparar
# =============================================================================
Ciclos_eje_H = []
Ciclos_eje_M = []
Labels = []

#%% =============================================================================
# Ahora itero [i] y proceso completamente a cada uno de los archivos seleccionados '''    
# =============================================================================

delta_t = [] #Base temporal
Idc = []     #Internal direc current en el generador de RF

for i in range(0,fa): 
    
    delta_t.append(1e-6/float(filenombres_muestra[i][11:14]))
    Idc.append(float(filenombres_muestra[i][7:9]))
    
    #Defino los vectores a procesar
    
    #Importo archivo xxxkHz_yyA_zzzMss_TM.dat
    muestra = np.loadtxt(rutas_de_carga[i][0],skiprows=3,dtype=float)
    t_m = (muestra[:,0]-muestra[0,0])*delta_t[i]    #1er columna 
    v_m = muestra[:,1]*0.001                        #CH1 
    v_r_m = muestra[:,2]*0.001                      #CH2

    
    #Imp: archivo xxxkHz_yyA_zzzMss_TM.dat_fondo (fondo: bobinas sin muestra)
    fondo = np.loadtxt(rutas_de_carga[i][1],skiprows=3,dtype=float)
    t_f = (fondo[:,0]-fondo[0,0])*delta_t[i]
    v_f = fondo[:,1]*0.001      
    v_r_f = fondo[:,2]*0.001    

    #Imp: xxxkHz_yyA_zzzMss_TM_cal.dat(calibracion: material paramagnetico, sin histéresis)
    calibracion = np.loadtxt(rutas_de_carga[i][2],skiprows=3,dtype=float)
    t_c = (calibracion[:,0]-calibracion[0,0])*delta_t[i]      #
    v_c = calibracion[:,1]*0.001      #
    v_r_c = calibracion[:,2]*0.001    #    
 
    par_m=[np.mean(v_r_m),(np.max(v_r_m)-np.min(v_r_m))/2] #Muestra
    par_f=[np.mean(v_r_f),(np.max(v_r_f)-np.min(v_r_f))/2] #Fondo
    par_c=[np.mean(v_r_c),(np.max(v_r_c)-np.min(v_r_c))/2] #Calibracion
    '''
    Ajusta las 3 referencias con funciones seno: V(t)=V0+A*sin(2*np.pi*f*t - phi)

    Estimacion valores iniciales para los ajustes (seeds)
    par_x[0] = Valor medio de la señal/Offset (i.e. constante aditiva)
    par_X[1] = Amplitud 
    Despues anexo:
    par_X[2] = Frecuencia
    par_X[3] = Fase

    par_x=[offset,amplitud]    
    '''
    #Suavizo las señales con fft_smooth     
    suave_m = fft_smooth(v_r_m, np.around(int(len(v_r_m)*6/1000)))
    suave_f = fft_smooth(v_r_f, np.around(int(len(v_r_f)*6/1000)))
    suave_c = fft_smooth(v_r_c, np.around(int(len(v_r_c)*6/1000)))


    indices_m = find_peaks(suave_m,height=0) #tupla
    indices_f = find_peaks(suave_f,height=0)
    indices_c = find_peaks(suave_c,height=0)

    #Para evitar problemas por errores de nomenclatura en los archivos mide el tiempo 
    #entre picos
    t_entre_max_m = np.mean(np.diff(t_m[indices_m[0]]))
    t_entre_max_f = np.mean(np.diff(t_f[indices_f[0]]))
    t_entre_max_c = np.mean(np.diff(t_c[indices_c[0]]))
    #indices_m[1]['peak_heights']: acceso al diccionario
    # par_x[2] = Frecuencia Semilla
    par_m.append(1/t_entre_max_m) #Frec de la Muestra
    par_f.append(1/t_entre_max_f) #Frec del Fondo
    par_c.append(1/t_entre_max_c) #Frec de la Calibracion

    #par_x=[offset,amplitud,frecuencia] 
 
    '''
    Fase inicial, a partir del tiempo del primer maximo
    '''
    #par_x[3] = Fase Semilla
    par_m.append(2*np.pi*par_m[2]*t_m[indices_m[0][0]] - np.pi/2) #Fase inic Muestra
    par_f.append(2*np.pi*par_f[2]*t_f[indices_f[0][0]] - np.pi/2) #Fase inic Fondo
    par_c.append(2*np.pi*par_c[2]*t_c[indices_c[0][0]] - np.pi/2) #Fase inic Calibracion

    #par_x=[offset,amplitud,frecuencia,fase] 

    seeds_m = [par_m[0],par_m[1],par_m[2],par_m[3]]    
    seeds_f = [par_f[0],par_f[1],par_f[2],par_f[3]]
    seeds_c = [par_c[0],par_c[1],par_c[2],par_c[3]]

    from scipy.optimize import curve_fit
    t_1=t_m
    y_1=suave_m
    coef_m, cov_m = curve_fit(sinusoide,t_1,y_1,seeds_m) #A,B,C,D 

    t_2=t_f
    y_2=suave_f
    coef_f, cov_f = curve_fit(sinusoide,t_2,y_2,seeds_f)

    t_3=t_c
    y_3=suave_c
    coef_c, cov_c = curve_fit(sinusoide,t_3,y_3,seeds_c)
    
    # =============================================================================
    # Graficos señales, ajustes y restos
    # =============================================================================
    
    # Muestra
    n_m = len(t_m)
    y_m = np.empty(n_m)
    for k in range(n_m):
        y_m[k] = sinusoide(t_m[k],coef_m[0],coef_m[1],coef_m[2],coef_m[3])

    #Para calcular el R^2 del ajuste: Metodo a usar: sklearn,
    from sklearn.metrics import r2_score 
    R_2_m = r2_score(v_r_m,y_m)
    resto_m = v_r_m - y_m #Resto: diferencia entre funcion original y la ajustada

    # Identifico Frecuencia y la Fase de la referencia de muestra:
   
    offset_m = coef_m[0]
    amplitud_m = coef_m[1]
    frecuencia_m = coef_m[2]
    frecuencia_m_kHz = frecuencia_m/1000
    fase_m = coef_m[3]    

    #Fondo
    n_f = len(t_f)
    y_f = np.empty(n_f)
    for j in range(n_f):
        y_f[j] = sinusoide(t_f[j],coef_f[0],coef_f[1],coef_f[2],coef_f[3])

    R_2_f = r2_score(v_r_f,y_f)
    resto_f = v_r_f - y_f
        
    offset_f = coef_f[0]
    amplitud_f = coef_f[1]
    frecuencia_f = coef_f[2]
    frecuencia_f_kHz = frecuencia_f/1000
    fase_f = coef_f[3]

    #Calibracion
    n_c = len(t_c)
    y_c = np.empty(n_c)
    for l in range(n_c):
        y_c[l] = sinusoide(t_c[l],coef_c[0],coef_c[1],coef_c[2],coef_c[3])

    R_2_c = r2_score(v_r_c,y_c)
    resto_c = v_r_c - y_c
    
    offset_c = coef_c[0]
    amplitud_c = coef_c[1]
    frecuencia_c = coef_c[2]
    frecuencia_c_kHz = frecuencia_c/1000
    fase_c = coef_c[3]
 
    #%% =============================================================================
    #  Si la diferencia entre frecuencias es muy gde => error
    # =============================================================================
    if(abs(frecuencia_m-frecuencia_f)/frecuencia_f >0.02):
        print('Error: incompatibilidad de frecuencias')
        print(f'Muestra: {frecuencia_m:6.3f} Hz')
        print(f'Fondo: {frecuencia_f:6.3f} Hz')
        print(f'Calibración: {frecuencia_c:6.3f} Hz')
    elif(abs(frecuencia_c-frecuencia_f)/frecuencia_f >0.02):
        print('Error: incompatibilidad de frecuencias')
        print(f'Muestra: {frecuencia_m:6.3f} Hz')
        print(f'Fondo: {frecuencia_f:6.3f} Hz')
        print(f'Calibración: {frecuencia_c:6.3f} Hz')
    elif(abs(frecuencia_c-frecuencia_m)/frecuencia_f >0.02):
        print('Error: incompatibilidad de frecuencias')
        print(f'Muestra: {frecuencia_m:6.3f} Hz')
        print(f'Fondo: {frecuencia_f:6.3f} Hz')
        print(f'Calibración: {frecuencia_c:6.3f} Hz')    
    
    
#%% Grafico las señales y sus ajustes
    #Muestra
    
    plt.plot(t_m*1e6, v_r_m,'o',label='Referencia de muestra')
    plt.plot(t_m*1e6,y_m,'r-',label='Ajuste de ref. de muestra')
    plt.plot(t_m*1e6,resto_m,'.', label='Restos')
    plt.xlabel('t ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.legend(loc='upper left',framealpha=1.0)
    #plt.text(max(t_m),min(v_r_m),'$R^2$ = {}'.format(R_2_m),bbox=dict(alpha=1.0), ha='right',va='bottom')
    plt.text(max(t_m)*1e6,-10,'$R^2$ = {}'.format(R_2_m),
             bbox=dict(alpha=0.9), ha='right',
             va='top')
    plt.text(0,-10,'Frecuencia: %0.2f kHz '%frecuencia_m_kHz,
             bbox=dict(alpha=0.9), ha='left',
             va='top')
    
    plt.text(max(t_m)*1e6,10,'Amplitud: %0.2f V '%amplitud_m,
             bbox=dict(alpha=0.9), ha='right',
             va='bottom')
    plt.axhspan(amplitud_m, -amplitud_m, facecolor='g', alpha=0.2)
    
    plt.title(filenombres_muestra[i][:-4] + ' - Señal de Muestra, ajuste y restos')
    plt.ylim(-12,12)
    plt.grid()
    #plt.savefig(filenombres_muestra[i][:-4] + ' - ajuste_muestra.png',dpi=500,bbox_inches='tight')
    plt.show()

#%% Fondo
    plt.plot(t_f*1e6, v_r_f,'go',label='Referencia de fondo')
    plt.plot(t_f*1e6,y_f,'y-',label='Ajuste de ref. de fondo')
    plt.plot(t_f*1e6,resto_f,'.', label='Restos')
    plt.xlabel('t ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.text(max(t_f)*1e6,-10,'$R^2$ = {}'.format(R_2_f),
             bbox=dict(alpha=0.9), ha='right',va='top')
    plt.text(0,-10,'Frecuencia: %0.2f kHz '%frecuencia_f_kHz,
             bbox=dict(alpha=0.9), ha='left',va='top')
    plt.text(max(t_f)*1e6,10,'Amplitud: %0.2f V '%amplitud_f,
             bbox=dict(alpha=0.9), ha='right',va='bottom')
    plt.axhspan(amplitud_f, -amplitud_f, facecolor='g', alpha=0.2)
    plt.title(filenombres_fondo[i] + ' - Señal de Fondo, ajuste y restos')
    plt.ylim(-12,12)
    plt.grid()
    #plt.savefig(filenombres_muestra[i][:-4] + ' - ajuste_fondo.png',dpi=500,bbox_inches='tight')
    plt.show()

#%% Calibracion 
    plt.plot(t_c*1e6, v_r_c,'o',label='Referencia de calibración ')
    plt.plot(t_c*1e6,y_c,'-m', label='Ajuste de ref. de calibración ')
    plt.plot(t_c*1e6,resto_c,'.', label='Restos')
    plt.xlabel('t ($\mu$s)')
    plt.ylabel('Amplitud de señal(V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.text(max(t_c)*1e6,-10,'$R^2$ = {}'.format(R_2_c),
             bbox=dict(alpha=0.9), ha='right',va='top')
    plt.text(0,-10,'Frecuencia: %0.2f kHz '%frecuencia_c_kHz,
             bbox=dict(alpha=0.9), ha='left',va='top')   
    plt.text(max(t_c)*1e6,10,'Amplitud: %0.2f V '%amplitud_c,
             bbox=dict(alpha=0.9), ha='right',va='bottom')
    plt.axhspan(amplitud_c, -amplitud_c, facecolor='g', alpha=0.2)    
    plt.title(filenombres_cal[i] + ' -  Señal de calibración , ajuste y restos')
    plt.ylim(-12,12)
    plt.grid()
    #plt.savefig(filenombres_cal[i] + ' - ajuste_calibracion.png',dpi=500,bbox_inches='tight')
    plt.show()

#%%
    #Desplazamiento temporal para poner en fase las referencias
    #Obtengo el modulo 2 pi de las fases y calcula el tiempo de fase 0    
    t_fase_m = np.mod(fase_m,2*np.pi)/(2*np.pi*frecuencia_m)
    t_fase_f = np.mod(fase_f,2*np.pi)/(2*np.pi*frecuencia_f)
    t_fase_c = np.mod(fase_c,2*np.pi)/(2*np.pi*frecuencia_c)

    #Desplazo en tiempo para que haya coicidencia de fase e/referencias.
    t_m_1 = t_m - t_fase_m 
    t_f_1 = t_f - t_fase_f 
    t_c_1 = t_c - t_fase_c 
    
    #Correccion por posible diferencia de frecuencias dilatando el tiempo del fondo 
    t_fm = t_f_1*frecuencia_f/frecuencia_m
    t_fc = t_f_1*frecuencia_f/frecuencia_c
    
    #Resta offset a las referencias y a los ajustes
    v_r_f_1 = v_r_f - offset_f
    v_r_m_1 = v_r_m - offset_m
    v_r_c_1 = v_r_c - offset_c
    
    y_f_1 = y_f - offset_f   
    y_m_1 = y_m - offset_m
    y_c_1 = y_c - offset_c
#%%   
    #Comparacion de las ref. de muestra y fondo 
    #desplazadas en tiempo y offset
    plt.plot(t_fm*1e6,v_r_f_1, label='Referencia de fondo',lw=0.7)
    plt.plot(t_m_1*1e6, v_r_m_1, label='Referencia de muestra',lw=0.7)
    plt.plot(t_fm*1e6, y_f_1, label='Ajuste Fondo',lw=0.7)
    plt.plot(t_m_1*1e6, y_m_1,label='Ajuste Muestra',lw=0.7)

    plt.grid()
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal(V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Comparación de ref. desplazadas y restados sus offsets')
    #plt.savefig(filenombres_muestra[i][:-4] + ' - comparacion_muestra_fondo.png',dpi=500,bbox_inches='tight')
    plt.show()
#%% 
    #Comparacion de las ref. de calibración y fondo 
    #desplazadas en tiempo y offset 
    plt.plot(t_fc*1e6,v_r_f_1, label='Ref. de fondo')
    plt.plot(t_c_1*1e6, v_r_c_1, label='Ref. de calibración')
    plt.plot(t_fc*1e6, y_f_1, label='Ajuste Ref. Fondo')
    plt.plot(t_c_1*1e6, y_c_1,label='Ajuste Calibración')

    plt.grid()
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Comparación de ref. desplazadas y restados sus offsets')
    #plt.savefig(filenombres_muestra[i][:-4] + ' - comparacion_calibracion_fondo.png',dpi=500,bbox_inches='tight')
    plt.show()        
#%%    
    '''
    Recorte para que ambras medidas tengan igual tiempo inicial y num entero de periodos
    Resta medida y fondo interpolando para que corresponda mejor el tiempo. 
    No trabaja mas con las medidas individuales, solo con la resta.
    Se toma la precaucion para los casos de trigger distinto e/ fondo y medida.
    '''
    v_m_1 = v_m     # Muestra a medir SAR 
    v_f_1 = v_f     # Fondo
    v_c_1 = v_c     # Calibracion campo (H) y magnetizacion (M)
    
    #Muestra
    t_min_m = t_fm[0]
    t_max_m = t_fm[-1]
    t_aux_m = t_m_1[np.nonzero((t_m_1>=t_min_m) & (t_m_1<= t_max_m))]

    interp_aux_m = np.interp(t_aux_m,t_fm,v_f_1)
    interp_m = np.zeros_like(v_m_1)   

    for w in range(0,len(t_m_1),1):
        index_min_m = np.argmin(abs(t_aux_m - t_m_1[w])) #obtengo el indice donde esta el minimo
        interp_m[w] = interp_aux_m[index_min_m] #defino c/ elemento de interp_m

    Resta_m = v_m_1 - interp_m 

    #Ahora lo mismo, pero para la calibracion
    t_min_c = t_fc[0]
    t_max_c = t_fc[-1]
    t_aux_c = t_c_1[np.nonzero((t_c_1>=t_min_c) & (t_c_1<= t_max_c))]

    interp_aux_c = np.interp(t_aux_c,t_fc,v_f_1)
    interp_c = np.zeros_like(v_c_1)    

    for x in range(0,len(t_c_1),1):
        index_min_c = np.argmin(abs(t_aux_c - t_c_1[x])) #obtengo el indice donde esta el minimo
        interp_c[x] = interp_aux_c[index_min_c] #defino c/ elemento de interp_m

    Resta_c = v_c_1 - interp_c 
#%%
    # =============================================================================
    # Comparacion de las medidas con el tiempo que pone las referencias en fase
    # =============================================================================
    
    plt.plot(t_m_1*1e6,v_m_1,'-',label='Señal de muestra',linewidth=1)
    plt.plot(t_fm*1e6,v_f_1,'-',label='Señal de fondo',linewidth=1)
    plt.plot(t_m_1*1e6,interp_m,'-',label='Interpolación de la muestra en fase con el fondo',linewidth=1)
    plt.xlim(0.0, t_m_1[int(len(t_m_1)/8)]*1e6) #Para zoom 
    plt.legend(loc='upper left',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Comparación de señales: Muestra, Fondo e Interpolacion')
    plt.grid()
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Comparacion_Fondo_Muestra_ref_en_fase.png',dpi=500,bbox_inches='tight')
    plt.show()
 #%%   
    plt.plot(t_c_1/1e6,v_c_1,'-',label='Señal de calibración ',linewidth=1)
    plt.plot(t_fm/1e6,v_f_1,'-',label='Señal de fondo',linewidth=1)
    plt.plot(t_c_1/1e6,interp_c,'-',label='Interpolación de la calibración  en fase con el fondo',linewidth=1)
    plt.xlim(0.0, t_c_1[int(len(t_c_1)/8)]/1e6) #para zoom
    plt.legend(loc='best',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Comparación de señales: Calibración, Fondo e Interpolación')
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Comparacion_Fondo_Calibracion_ref_en_fase.png',dpi=500,bbox_inches='tight')
    plt.grid()
    plt.show()
#%% 
    #Resta de señales
    plt.plot(t_m_1*1e6,v_m_1,'-',label='Señal de muestra',linewidth=1)
    plt.plot(t_m_1*1e6,interp_m,'-',label='Interpolación de la muestra en fase con el fondo',linewidth=1)
    plt.plot(t_m_1*1e6,Resta_m,'-',label='$Muestra - Fondo$',linewidth=1)
    plt.legend(loc='best',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.title(filenombres_muestra[i][:-4] + ' - Resta de señales: Muestra')
    plt.xlim(0.0,t_m_1[int(len(t_m_1)/8)]*1e6)
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Resta_de_señales.png',dpi=500,bbox_inches='tight')
    plt.grid()
    plt.show() 
#%%    
    plt.plot(t_c_1*1e6,v_c_1,'-',label='Señal de calibración ',linewidth=1)
    plt.plot(t_c_1*1e6,interp_c,'-',label='Interpolación del calibración  en fase con el fondo',linewidth=1)
    plt.plot(t_c_1*1e6,Resta_c,'-',label='calibración  $-$ Fondo',linewidth=1)
    plt.legend(loc='best',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.title(filenombres_muestra[i][:-4] + ' - Resta de señales')
    plt.xlim(0.0,t_m_1[int(len(t_m_1)/8)]*1e6) #para zoom
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Resta_de_señales.png',dpi=500,bbox_inches='tight')
    plt.grid()
    plt.show()
#%%
    #Filtrado
    if FILTRARMUESTRA == 0:

        #Filtro por Fourier la muestra para suavizarla    
        freq_m = np.around(int(len(v_r_m)/5))
        Resta_m_3 = fft_smooth(Resta_m, freq_m)
        v_r_m_3 = fft_smooth(v_r_m,freq_m)
        t_m_3 = t_m_1
        ajuste_m_3 = y_m_1

        #Filtro por Fourier la calibracion
        freq_c = np.around(int(len(v_r_c)/5))
        Resta_c_3 = fft_smooth(Resta_c, freq_c)
        v_r_c_3 = fft_smooth(v_r_c,freq_c)
        t_c_3 = t_c_1
        ajuste_c_3 = y_c_1        
#%%     Controlo que el suavizado final sea satisfactorio
        
        ax1 = plt.subplot(211)
        plt.plot(t_m_1*1e6,v_r_m_1/max(v_r_m_3), '-',
                 label='Ref. de muestra',linewidth=1)
        plt.plot(t_m_3*1e6,v_r_m_3/max(v_r_m_3), '-',
                 label='Ref. de muestra filtrada',linewidth=1)
        plt.legend(loc='lower left',framealpha=1.0)
        plt.setp(ax1.get_xticklabels(),visible=False)
        #plt.xlabel('Tiempo (s)')
        plt.ylabel('u.a.')
        plt.grid()
        plt.title(filenombres_muestra[i][:-4] 
        + ' - Quita de ruido por FFT: Muestra')
        
        ax2 = plt.subplot(212,sharex=ax1)
        plt.plot(t_m_1*1e6,Resta_m/max(Resta_m),'-',label='Resta de señales: $(M-F)$',linewidth=1)
        plt.plot(t_m_3*1e6,Resta_m_3/max(Resta_m_3),'-',
                 label='$(M-F)$ filtrada',linewidth=1)        
        plt.legend(loc='lower left',framealpha=1.0)
        plt.xlabel('t ($\mu$s)')
        plt.ylabel('u.a.')
        plt.grid()
        plt.tight_layout()
        #plt.savefig(filenombres_muestra[i][:-4] + ' - Quita_de_ruido_por_FFT_Muestra.png',dpi=500,bbox_inches='tight')
#%%        
        ax3 = plt.subplot(211)
        plt.plot(t_c_1*1e6,v_r_c_1/max(v_r_c_3), 
                 label='Ref. de calibración ',linewidth=1)
        plt.plot(t_c_3*1e6,v_r_c_3/max(v_r_c_3),
                 label='Ref. de calibración  filtrada',linewidth=1)
        plt.legend(loc='lower left',framealpha=1.0)
        plt.setp(ax3.get_xticklabels(),visible=False)
        plt.ylabel('u.a.')
        plt.title(filenombres_cal[i] 
        + ' - Quita de ruido por FFT: Ref. de calibración ')
        plt.grid()
        
        ax4 = plt.subplot(212,sharex=ax3)
        plt.plot(t_c_1*1e6,Resta_c/max(Resta_c_3),
                 label='Resta de señales: $(C-F)$',linewidth=1)
        plt.plot(t_c_3*1e6,Resta_c_3/max(Resta_c_3),
                 label='$(C-F)$ filtrada',linewidth=1)
        plt.xlabel('Tiempo ($\mu$s)')
        plt.ylabel('u.a.')
        plt.legend(loc='lower left',framealpha=1.0)
        plt.grid()
        plt.tight_layout()
        #plt.savefig(filenombres_muestra[i][:-4] + ' - Quita_de_ruido_por_FFT_Calibracion.png',dpi=500,bbox_inches='tight')
        plt.show()
#%%          
        #Diferencia entre señal sin ruido y señal. Guarda el peor valor
        dif_resta_m = Resta_m_3 - Resta_m
        dif_resta_c = Resta_c_3 - Resta_c
        peor_dif.append(max([np.mean(abs(dif_resta_m))/max(Resta_m),np.mean(abs(dif_resta_c))/max(Resta_c)]))
#%%        
    # =============================================================================
    # Recorta un número entero de períodos o ciclos
    # =============================================================================
    #Numero de ciclos
    N_ciclos_m =  np.floor((t_m_3[-1] - t_m_3[0])*frecuencia_m)   
    N_ciclos_c =  np.floor((t_c_3[-1] - t_c_3[0])*frecuencia_c)   
    
    #Indices ciclo
        #Muestra
    indices_ciclo_m = np.argwhere(t_m_3 < t_m_3[0] + N_ciclos_m/frecuencia_m)  
    largo_m = indices_ciclo_m[-1][0]

    if np.mod(largo_m,N_ciclos_m) == 0:
        largo_m = largo_m - np.mod(largo_m,N_ciclos_m)
    elif np.mod(largo_m,N_ciclos_m) <= 0.5:
        largo_m = largo_m - np.mod(largo_m,N_ciclos_m)
    
    else:
        largo_m = largo_m + N_ciclos_m - np.mod(largo_m,N_ciclos_m)
 
        #Calibracion
    indices_ciclo_c = np.argwhere(t_c_3 < t_c_3[0] + N_ciclos_c/frecuencia_c)    
    largo_c = indices_ciclo_c[-1][0]

    if np.mod(largo_c,N_ciclos_c) == 0:
        largo_c = largo_c - np.mod(largo_c,N_ciclos_c)
    elif np.mod(largo_c,N_ciclos_c) <= 0.5:
        largo_c = largo_c - np.mod(largo_c,N_ciclos_c)
    
    else:
        largo_c = largo_c + N_ciclos_c - np.mod(largo_c,N_ciclos_c)
         
    #Recorto 1/2 periodo al ppio y al final
    indices_recortado_m = np.arange(np.ceil(largo_m/N_ciclos_m/2),np.ceil(largo_m-largo_m/N_ciclos_m/2),1,dtype=int)    
    N_ciclos_m = N_ciclos_m -1    
    indices_recortado_c = np.arange(np.ceil(largo_c/N_ciclos_c/2),np.ceil(largo_c-largo_c/N_ciclos_c/2),1,dtype=int)    
    N_ciclos_c = N_ciclos_c -1    

    #Se recortan los vectores    
    t_m_4 = t_m_3[indices_recortado_m]
    v_r_m_4 = v_r_m_3[indices_recortado_m]
    Resta_m_4 = Resta_m_3[indices_recortado_m]
    ajuste_m_4 = ajuste_m_3[indices_recortado_m]

    t_c_4 = t_c_3[indices_recortado_c]
    v_r_c_4 = v_r_c_3[indices_recortado_c]
    Resta_c_4   = Resta_c_3[indices_recortado_c]
    ajuste_c_4 = ajuste_c_3[indices_recortado_c]    
 #%%     
    '''Referencia de muestra sin valor medio y su ajuste y resta de señal de muestra y fondo,
    recortadas para que tengan un numero entero de ciclos'''
   
    plt.plot(t_m_4*1e6,v_r_m_4,label='Referencia de muestra',linewidth=1)
    plt.plot(t_m_4*1e6,Resta_m_4,label='($M-F$)',linewidth=1)
    plt.legend(loc='lower left',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.grid()
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Ref. de muestra desplazada s/ valor medio, y resta de señales.\n Nº entero de periodos') 
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Ref_de_Muestra_desplazada_sin_vm_resta_de_señales.png',dpi=500,bbox_inches='tight')
    plt.show()    
    #%%
    plt.plot(t_c_4*1e6,v_r_c_4,label='Referencia de calibración ',linewidth=1)
    plt.plot(t_c_4*1e6,Resta_c_4,label='Resta de calibración ',linewidth=1)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.5)
    plt.legend(loc='lower left',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Señal (V)')
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Ref. de calibración  desplazada s/ valor medio y resta de señales\n N entero de periodos')
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Ref_de_Callbracion_desplazada_sin_vm_resta_de_señales.png',dpi=500,bbox_inches='tight')
    plt.show()    

    #%% =============================================================================
    # Campo y Magnetizacion: se integran las funciones de referencia y resta
    # =============================================================================

    #Ultimos ajustes de las referencias: par_x_f= [Offset, Amplitud, Frecuencia, Fase inicial]    
    par_m_f =[np.mean(v_r_m_4),(np.max(v_r_m_4)-np.min(v_r_m_4))/2,frecuencia_m,2*np.pi*frecuencia_m*t_m_4[0]-(np.pi/2)] #Muestra
    par_c_f =[np.mean(v_r_c_4),(np.max(v_r_c_4)-np.min(v_r_c_4))/2,frecuencia_c,2*np.pi*frecuencia_c*t_c_4[0]-(np.pi/2)] #Calibración 

    #Ajuste y obtencion de coeficientes
    ajuste_final_m , cov_final_m  = curve_fit(sinusoide,t_m_4,v_r_m_4,par_m_f) #A,B,C,D 
    frecuencia_final_m = ajuste_final_m[2]

    ajuste_final_c , cov_final_c  = curve_fit(sinusoide,t_c_4,v_r_c_4,par_c_f)  
    frecuencia_final_c = ajuste_final_c[2]

    '''Promedio los ciclos antes de integrar'''

    t_final_m = t_m_4[np.nonzero(t_m_4<=t_m_4[0]+(1/frecuencia_final_m))]
    t_final_c = t_c_4[np.nonzero(t_c_4<=t_c_4[0]+(1/frecuencia_final_c))]

    fondo_1_m = np.zeros_like(t_final_m)
    fondo_1_c = np.zeros_like(t_final_c)  

    muestra_0 = np.zeros_like(t_final_m)
    cal_0 = np.zeros_like(t_final_c)

    ajuste_m = np.zeros_like(t_final_m) 
    ajuste_c = np.zeros_like(t_final_c) 
    # Cambie la interpolacion de numpy a la de scipy para poder hacer interpolacion CubicSpline
    
    '''Muestra'''
    for m in range(1,int(N_ciclos_m)+1,1):
        if t_final_m[-1] + (m-1)/frecuencia_final_m < t_m_4[-1]:
            interp_m_2 = sc.interpolate.interp1d(t_m_4,v_r_m_4)
            fondo_1_m = fondo_1_m + interp_m_2(t_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m
            
            interp_m_3 = sc.interpolate.interp1d(t_m_4,Resta_m_4)
            muestra_0 = muestra_0 + interp_m_3(t_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m
            
            interp_m_4 = sc.interpolate.interp1d(t_m_4,ajuste_m_4)
            ajuste_m = ajuste_m + interp_m_4(t_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m
        else:
            interp_m_5 = sc.interpolate.CubicSpline(t_m_4,v_r_m_4) #CubicSpline
            fondo_1_m = fondo_1_m + interp_m_5(t_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m
        
            interp_m_6 = sc.interpolate.CubicSpline(t_m_4,Resta_m_4)
            muestra_0 = muestra_0 + interp_m_6(t_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m
            
            interp_m_7 = sc.interpolate.CubicSpline(t_m_4,ajuste_m_4)
            ajuste_m = ajuste_m + interp_m_7(t_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m
            
    '''Calibración '''
    for m in range(1,int(N_ciclos_c)+1,1):
        if t_final_c[-1] + (m-1)/frecuencia_final_c < t_c_4[-1]:
            interp_c_2 = sc.interpolate.interp1d(t_c_4,v_r_c_4)
            fondo_1_c = fondo_1_c + interp_c_2(t_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c
            
            interp_c_3 = sc.interpolate.interp1d(t_c_4,Resta_c_4)
            cal_0 = cal_0 + interp_c_3(t_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c
            
            interp_c_4 = sc.interpolate.interp1d(t_c_4,ajuste_c_4)
            ajuste_c = ajuste_c + interp_c_4(t_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c
                
        else:
            interp_c_5 = sc.interpolate.CubicSpline(t_c_4,v_r_c_4) #CubicSpline
            fondo_1_c = fondo_1_c + interp_c_5(t_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c
        
            interp_c_6 = sc.interpolate.CubicSpline(t_c_4,Resta_c_4)
            cal_0 = cal_0 + interp_c_6(t_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c
        
            interp_c_7 = sc.interpolate.CubicSpline(t_c_4,ajuste_c_4)
            ajuste_c = ajuste_c + interp_c_7(t_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c
    
    #Quito valor medio a la resta y al fondo
    Rm = muestra_0 - np.mean(muestra_0)
    fondo_1_m = fondo_1_m - np.mean(fondo_1_m)

    Rc = cal_0 - np.mean(cal_0)
    fondo_1_c = fondo_1_c - np.mean(fondo_1_c)        
        
    #Paso temporal 
    delta_t_m = (t_final_m[-1] -t_final_m[0])/(len(t_final_m))        
    delta_t_c = (t_final_c[-1] -t_final_c[0])/(len(t_final_c))   
    #%% 
    plt.plot(t_final_m*1e6,fondo_1_m,label='Muestra')
    plt.plot(t_final_c*1e6,fondo_1_c,label='Calibración ')
    plt.legend(loc='lower left',framealpha=1.0)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('Amplitud de la señal (V)')
    plt.grid()
    plt.title(filenombres_muestra[i][:-4] 
    + ' - Señal de Fondo: Muestra y calibración')
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Señal_Fondo_Muestra_Calibracion.png',dpi=500,bbox_inches='tight')
    plt.show()
    #%%      
    '''Calibración :
    Calcula sumas acumuladas y convierte a campo(H) y magnetizacion(M) 
    Cte que dimensionaliza al campo en A/m a partir de la calibración realizada 
    sobre la bobina del RF
    '''
    C_norm_campo = Idc[i]*pendiente_HvsI + ordenada_HvsI
    
    #Campo en volt*segundo, falta llevar a la amplitud conocida.
    campo_ua0_c = delta_t_c*sc.integrate.cumtrapz(fondo_1_c,initial=0)
    # Centrado en cero
    campo_ua_c = campo_ua0_c-np.mean(campo_ua0_c)
    
    # Lo mismo con el campo obtenido del ajuste de la referencia
    campo_fit_ua_c = delta_t_c*sc.integrate.cumtrapz(ajuste_c - np.mean(ajuste_c),initial=0)
    campo_fit_ua_c = campo_fit_ua_c - np.mean(campo_fit_ua_c)
    
    #Magnetización en volt*segundo, falta el factor geométrico
    magnetizacion_ua0_c = delta_t_c*sc.integrate.cumtrapz(Rc,initial=0)
    
    #Centrado en cero
    magnetizacion_ua_c = magnetizacion_ua0_c - np.mean(magnetizacion_ua0_c)  
    #Doy unidades al campo con cte: C_norm_campo
    campo_c = campo_ua_c*C_norm_campo/max(campo_ua_c)  
    campo_fit_ua_c = campo_fit_ua_c*C_norm_campo/max(campo_fit_ua_c)
    
    #%% Campo y magnetizacion normalizados - Calibración
    plt.plot(t_final_c*1e6,campo_ua_c/max(campo_ua_c),label='Campo')
    plt.plot(t_final_c*1e6,magnetizacion_ua_c/max(magnetizacion_ua_c),
             label='Magnetización')
    plt.legend(loc='best', borderaxespad=0.5)
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('u.a.')
    plt.title(filenombres_muestra[i][:-4] + ' - Campo y Magnetización normalizados - Calibración ')
    #plt.savefig(filenombres_muestra[i][:-4] + ' - Campo y magnetizacion normalizados_Calibracion.png',dpi=500,bbox_inches='tight')
    plt.show()
    #%%   
    #Ajusto recta a la magnetizacion del paramagneto
    recta = np.polyfit(campo_c,magnetizacion_ua_c,1)
    polaridad = np.sign(recta[0])
    pendiente = recta[0]
    ordenada = recta[1]
    t_recta = np.linspace(min(campo_c),max(campo_c),len(campo_c))
    
    #Calibracion u.a. a Ampere
    calibracion = xi_patron_vol/pendiente
    #%%
    # =============================================================================
    # Ciclo de Histeresis del paramagneto (idealmente sin area)   
    # =============================================================================
    plt.plot(campo_c,magnetizacion_ua_c,label='Gd$_2$O$_3$')
    plt.plot(t_recta,pendiente*t_recta+ordenada,label='Ajuste lineal')
    plt.legend(loc='best', borderaxespad=0.5)
    plt.xlabel('H (Oe)')
    plt.ylabel('u.a.')
    plt.grid()
    plt.title(filenombres_cal[i] + ' - Ciclo de histéresis del la calibración ')
    #plt.savefig(filenombres_cal[i] + 'Ciclo_Histeresis_calibracion.png',dpi=500,bbox_inches='tight')
    plt.show()
    #%%  
    '''Muestra (repito rutina):
    Calcula sumas acumuladas y convierte a campo(H) y magnetizacion(M) 
    Cte que dimensionaliza al campo en A/m a partir de la calibración  realizada 
    sobre la bobina del RF
    '''
    #Campo en volt*segundo, falta llevar a la amplitud conocida.
    campo_ua0_m = delta_t_m*sc.integrate.cumtrapz(fondo_1_m,initial=0)
    # Centrado en cero
    campo_ua_m = campo_ua0_c - np.mean(campo_ua0_c)
   
    # Lo mismo con el campo obtenido del ajuste de la referencia (por eso 'fit')
    campo_fit_ua_m = delta_t_m*sc.integrate.cumtrapz(ajuste_m - np.mean(ajuste_m),initial=0)
    campo_fit_ua_m = campo_fit_ua_m - np.mean(campo_fit_ua_m)
   
    '''Magnetización en volt*segundo, falta el factor geométrico'''
    magnetizacion_ua0_m = delta_t_m*sc.integrate.cumtrapz(Rm,initial=0)
   
    '''Centrado en cero'''
    magnetizacion_ua_m = magnetizacion_ua0_m - np.mean(magnetizacion_ua0_m)  
   
    '''Doy unidades al campo, al campo ajustado de la ref, y a la magnetizacion'''
    campo_m = C_norm_campo*(campo_ua_m/max(campo_ua_m) ) 
    
    campo_fit_m = C_norm_campo*(campo_fit_ua_m/max(campo_fit_ua_m))
    
    magnetizacion_m = calibracion*magnetizacion_ua_m

    #%% =============================================================================
    #Para salvar posibles incongruencia en largo de arrays (agregado 08-01-2021)
    # =============================================================================
    if len(t_final_m) > len(campo_ua_m):
        t_final_m = np.resize(t_final_m,len(campo_ua_m))
    elif len(t_final_m) < len(campo_ua_m):
        campo_ua_m = np.resize(campo_ua_m,len(t_final_m))
    if len(magnetizacion_m) > len(t_final_m) :
        magnetizacion_m = np.resize(magnetizacion_m,len(t_final_m))
    # elif len(t_final_m) < len(magnetizacion_m):
    #     campo_ua_m = np.resize(magnetizacion_m,len(t_final_m))        

    #%% =============================================================================
    #Grafico campo y magnetizacion vs. t (desplazo tiempo hasta valores positivos)
    # =============================================================================
    if min(t_final_m)<0:
        t_desp_m =  (t_final_m + abs(min(t_final_m)))*1e6 #(paso a microseg)
    elif min(t_final_m)>0:
        t_desp_m =  (t_final_m - abs(min(t_final_m)))*1e6 
        
    plt.plot(t_desp_m,campo_c/max(campo_c),label='Campo Magnetico H')
    plt.plot(t_desp_m,magnetizacion_m/max(magnetizacion_m),label='Magnetización M')
    plt.legend( loc='best', borderaxespad=0.5)
    plt.title(filenombres_muestra[i][:-4] + ' - Campo y magnetización normalizados')
    plt.xlabel('Tiempo ($\mu$s)')
    plt.ylabel('u.a.')
    plt.grid()
    plt.show()
# =============================================================================
#     #%%
#     '''Campo y Magnetizacion finales:
#     Se mantiene un campo adicional generado a partir del ajuste senoidal de la referencia
#     para hacer comparaciones (campo_fit).  
#     '''  
# 
#     plt.plot(campo_m,magnetizacion_m, label=filenombres_muestra[i][:-4])
#     
#     plt.title('Ciclo de histéresis - ' + str( filenombres_muestra[i][:-4]))
#     plt.xlim(-1.1*max(abs(campo_m)),1.1*max(abs(campo_m)))
#     plt.ylim(-1.1*max(abs(magnetizacion_m)),1.1*max(abs(magnetizacion_m)))
#     plt.xlabel('Campo (A/m)')
#     plt.ylabel('Magnetizacion (A/m)')
#     plt.text(0,0,'T= {}'.format(Temperaturas[i]),
#               bbox=dict(alpha=0.9), ha='center',
#               va='top')
#     plt.grid()
#     plt.legend()
# =============================================================================
    #%%    
    '''Exporto ciclos de histeresis en ascii: tiempo | campo_m | magnetizacion_m |'''
    from astropy.io import ascii
    from astropy.table import Table, Column, MaskedColumn
    col0 = t_desp_m
    col1 = campo_m
    col2 = magnetizacion_m
    ciclo_out = Table([col0, col1, col2])
    
    encabezado = ['Tiempo (us)','Campo (A/m)', 'Magnetización (A/m)']
    formato = {'Tiempo (us)':'%f' ,'Campo (A/m)':'%f','Magnetización (A/m)':'%f'} 
    ascii.write(ciclo_out,filenombres_muestra[i][:-4] + '  - Ciclo de histéresis.dat',
                names=encabezado,overwrite=True,delimiter='\t',formats=formato)
    
    #%%
    # =============================================================================
    # Calculo Campo Coercitivo (Hc) y Magnetizacion Remanente (Mr) 
    # =============================================================================
    m = magnetizacion_m
    h = campo_m

    Hc = [] #np.zeros(15) #Coercitivo  
    Mr = []#np.zeros(15) #Remanente   
    
    for z in range(0,len(m)-1):
        if ((m[z]>0 and m[z+1]<0) or (m[z]<0 and m[z+1]>0)): #M remanente
            Hc.append(abs(h[z] - m[z]*(h[z+1] - h[z])/(m[z+1]-m[z])))
            
        if((h[z]>0 and h[z+1]<0) or (h[z]<0 and h[z+1]>0)):  #H coercitivo
            Mr.append(abs(m[z] - h[z]*(m[z+1] - m[z])/(h[z+1]-h[z])))
            
    Hc_mean = np.mean(Hc)
    Hc_error = np.std(Hc)
    
    Mr_mean = np.mean(Mr)
    Mr_error = np.std(Mr)
    
    print(f'Hc = {Hc_mean:.2f} (+/-) {Hc_error:.2f} (A/m)')
    print(f'Mr = {Mr_mean:.2f} (+/-) {Mr_error:.2f} (A/m)')
    
    #%%=============================================================================
    # Grafico el 1er ciclo de histéresis de la muestra
    # =============================================================================
    plt.plot(campo_m,magnetizacion_m, label='Archivo: %s' %filenombres_muestra[i])
    plt.axhline(y=0, color='k', linestyle='-')
    plt.axvline(x=0, color='k', linestyle='-')
    plt.plot([0,0],[Mr_mean,-Mr_mean],'D', label='Magnetizacion remanente: %.2f A/m' %Mr_mean)
    plt.plot([Hc_mean,-Hc_mean],[0,0],'D', label='Campo coercitivo: %.2f A/m' %Hc_mean)
    #plt.plot(campo_m[0],magnetizacion_m[0],'o', label=' %.2f A/m' %magnetizacion_m[0])
    #plt.fill_between(campo_m,magnetizacion_m, color='#539ecd')
    plt.title('Ciclo de histéresis - '+ str( filenombres_muestra[i] ))
    plt.xlim(-1.1*max(abs(campo_m)),1.1*max(abs(campo_m)))
    plt.ylim(-1.1*max(abs(magnetizacion_m)),1.1*max(abs(magnetizacion_m)))#es en terminos de 110% del valor maximo en cada eje, redondeado
    plt.xlabel('Campo (A/m)')
    plt.ylabel('Magnetizacion (A/m)')
    plt.grid()
    plt.legend()
    plt.show()
#%%
    '''Determinacion de areas: 
    Armo vector auxiliar desplazado a valores positivos'''
    magnetizacion_m_des = magnetizacion_m + 2*abs(min(magnetizacion_m))
    
    '''Area del lazo'''
    Area_ciclo = abs(sc.integrate.trapz(magnetizacion_m_des,campo_m)) 
    #print('Area del ciclo: %f' %Area_ciclo)
    
    '''Calculo de potencia disipada SAR'''
    sar = mu_0*abs(Area_ciclo)*frecuencia_final_m/(1000*concentracion) 
    #multiplico por 1000 para pasar de concentracion en g/l a g/m^3
    print(f'El SAR de la muestra es: {sar:f} (W/g)' )
      
    '''Salidas importantes: las listas tienen tantos elmentos como archivos seleccionados'''    
    
    ''' Corriente del resonador: definida en linea 198'''
        
    ''' Frecuencia de la referencia en la medida de la muestra'''
    Frecuencia_muestra_kHz.append(frecuencia_final_m/1000)
    
    ''' Frecuencia de la referencia en la medida del fondo'''
    Frecuencia_fondo_kHz.append(frecuencia_f/1000)
    
    ''' Specific Absorption Rate'''
    SAR.append(sar)
    
    ''' Campo maximo en kA/m'''
    Campo_maximo_kAm.append(C_norm_campo/1000)
    
    ''' Campo coercitivo en kA/m'''
    Coercitividad_kAm.append(Hc_mean/1000)
    
    ''' Magnetizacion remanente en kA/m'''
    Magnetizacion_remanente_kAm.append(Mr_mean/1000)

    '''Grafico con SAR y demas info'''
    #cierro ciclo
    campo_m[-1] = campo_m[1]
    magnetizacion_m[-1] = magnetizacion_m[1]
    
    
    Mag_a_Hmax.append(magnetizacion_m[np.argwhere(campo_m==max(campo_m))])
#%% Ciclo de Histereris final individual
    
    #fig=plt.figure(figsize=(8,8))
    plt.plot(campo_m,magnetizacion_m,'-')
    plt.axhline(y=0, color='k', linestyle='-',linewidth=0.5)
    plt.axvline(x=0, color='k', linestyle='-',linewidth=0.5)
    plt.plot([0,0],[Mr_mean,-Mr_mean],'D', label='Magnetizacion remanente: %.2f A/m' %Mr_mean)
    plt.plot([Hc_mean,-Hc_mean],[0,0],'D', label='Campo coercitivo: %.2f $A/m$' %Hc_mean)
    #plt.fill_between(campo_m[-2:0:-1],magnetizacion_m[-2:0:-1], color='#539ecd',alpha=1)
    plt.title('Ciclo de histéresis - Muestra: '+ filenombres_muestra[i][:-4] ) 
    plt.grid()    
    #plt.xlim(-1.1*max(abs(campo_m)),1.1*max(abs(campo_m)))
    #plt.ylim(-1.1*max(abs(magnetizacion_m)),1.1*max(abs(magnetizacion_m)))
    
    plt.xlabel('Campo ($A/m$)')
    plt.ylabel('Magnetización ($A/m$)')
    plt.legend(loc='best',framealpha=1.0)
    plt.ylabel('Magnetizacion (A/m)')
    plt.text(min(campo_m),0,'T: {} $ºC$ \n SAR: {:.2f} $W/g$'.format(Temperaturas[i], sar),        bbox=dict(boxstyle="round", alpha=0.8, lw=0.7),
              ha='left', va='center')

    plt.savefig(filenombres_muestra[i] + ' - Ciclo de histéresis.png',dpi=500,bbox_inches='tight')
    plt.show()   
    
    #%% =============================================================================
    # Agrego los ciclos a las listas definidas para graficarlos juntos a continuacion
    # =============================================================================
    Ciclos_eje_H.append(campo_m)
    Ciclos_eje_M.append(magnetizacion_m)
    
#%% Colormap para el grafico con todos los ciclos
import matplotlib.colors as colors
import matplotlib.cm as cmx    
cm = plt.get_cmap('rainbow') 
cNorm  = colors.Normalize(vmin=0, vmax=fa)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    
#%%    
for  i in range(0,fa):
    colorVal = scalarMap.to_rgba(i)   
    
    plt.plot(Ciclos_eje_H[i],Ciclos_eje_M[i],'-', lw=1.0, color=colorVal,
             label='%.1f ºC' %(Temperaturas[i]))
    plt.axhline(y=0, color='k', ls='-',lw=0.5)
    plt.axvline(x=0, color='k', ls='-',lw=0.5)
    plt.title(nombre_salida +' - Ciclos de histéresis a distintas temperaturas'  )  
    #plt.xlim(-1.1*max(abs(campo_m)),1.1*max(abs(campo_m)))
    #plt.ylim(-1.1*max(abs(magnetizacion_m)),1.1*max(abs(magnetizacion_m)))
    plt.xlabel('Campo ($A/m$)')
    plt.ylabel('Magnetización ($A/m$)')
    plt.legend(loc='best',bbox_to_anchor=(1,1),framealpha=1.0,ncol=2)
    plt.ylabel('Magnetizacion (A/m)')
plt.grid()
#plt.savefig(nombre_salida +' - Ciclos de histéresis a distintas T.png',dpi=500,bbox_inches='tight')
plt.show()
#%%
# =============================================================================
# Coercitividad y remanencia en funcion de la Temperatura
# =============================================================================
ax8=plt.subplot(311)
plt.plot(Temperaturas,Coercitividad_kAm,'-D',c='#2ca02c')
plt.setp(ax8.get_xticklabels(),visible=False) 
plt.ylabel('$H_c$ ($kA/m$)')
plt.grid() 
plt.title(nombre_salida +' - Campo coercitivo y Magnetizacion vs. Temperatura')  
    
ax9=plt.subplot(312,sharex=ax8)
plt.setp(ax9.get_xticklabels(),visible=False)
plt.plot(Temperaturas,Magnetizacion_remanente_kAm,'-D', c='#ff7f0e' )
plt.ylabel('$M_r$ ($kA/m$)')
plt.xlabel('Temperatura ($ºC$)')
plt.grid()    

ax10=plt.subplot(313,sharex= ax9)
plt.plot(Temperaturas,SAR,'-D', c='#d62728' )
plt.setp(ax8.get_xticklabels(),visible=False)
plt.ylabel('$SAR$ ($W/g$)')
plt.xlabel('Temperatura ($ºC$)')
plt.grid() 

plt.tight_layout()
#plt.savefig(nombre_salida +' - Hc, Mr y SAR vs. Temperatura.png',dpi=500,bbox_inches='tight')
plt.show() 
plt.show()

#%% Para poner la fecha del procesamiento
from datetime import datetime
fecha = datetime.today().strftime('%d - %m - %Y - %H:%M:%S')



#%%

vector = []
for i in range(0,fa):
    vector.append(Mag_a_Hmax[i][0][0])
plt.plot(Temperaturas,vector,'.-')
plt.ylabel('$Magnetizacion a campo maximo$ ($W/g$)')
plt.xlabel('Temperatura ($ºC$)')
plt.xlim(-20,10)
plt.grid() 



plt.savefig('Nueva.png',dpi=500,bbox_inches='tight')
#%%

'''
Archivo de salida: aca se itera sobre los k archivos que seleccione 
Salidas de interes: armo las listas y les voy agregando los datos al final de cada 
iteracion
'''    
#Encabezado del archivo de salida 
encabezado_salida = ['Nombre del archivo procesado','Tipo de muestra','Frecuencia (kHz)',
                      'Temperatura (ºC)' ,'SAR (W/g)','Campo Maximo (kA/m)',
                      'Coercitividad (kA/m)','Magnetizacion Remanente (kA/m)', 
                      'Peor quita de ruido porcentual']
     
#Las columnas deben ser listas con largo=num de archivos seleccionados 
col_0 = nombres_archivos
col_1 = tipo_muestra
col_2 = Frecuencia_muestra_kHz
col_3 = SAR          
col_4 = Campo_maximo_kAm   
col_5 = Coercitividad_kAm            
col_6 = Magnetizacion_remanente_kAm  
col_7 = peor_dif
col_8 = Temperaturas
#Armo la tabla    
salida = Table([col_0, col_1, col_2, col_8, col_3, col_4, col_5, col_6, col_7]) 
formato_salida = {'Nombre del archivo procesado':'%s','Tipo de muestra': '%s',
                  'Frecuencia (kHz)':'%f','Temperatura (ºC)': '%.1f','SAR (W/g)':'%f',
                  'Campo Maximo (kA/m)':'%f','Coercitividad (kA/m)': '%f',
                  'Magnetizacion Remanente (kA/m)': '%f', 
                  'Peor quita de ruido porcentual': '%f'} 
#Agrego la fecha y la concentracion
salida.meta['comments'] = [fecha, 'Concentración de la muestra: '
           +str(concentracion)+ ' g/l']

ascii.write(salida,'Resultados_Esar_'+ nombre_salida +'.dat',names=encabezado_salida,
            overwrite=True,delimiter='\t',formats=formato_salida)


#%% Tiempo de procesamiento
end_time = time.time()
print(f'Concentración: {(concentracion):.2f} g/l')
print(f'Tiempo de ejecución del script: {(end_time-start_time):6.3f} s.')
print(fecha)



