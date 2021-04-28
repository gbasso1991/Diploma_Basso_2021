# -*- coding: utf-8 -*-
"""
Psyco Holiday - Giuliano Andrés Basso - 19 Nov 2020 - 

Este script permite obtener ver la relacion de proporcionalidad entre la corriente Idc
y la amplitud de las señales de referencia analizados con Planet Caravan.
La idea es obtener el valor de Idc directamente de los archivos procesados y asi independizar 
los resultados de la nomenclatura de los mencionados archivos. Esto aportaria mayor robustez 
al procesamiento de datos pues minimiza la posibilidad de error humano.

Toma un archivo con formato de nombre 'xxxkHz_yyA_zzz_Mss_TM.dat'
Siendo:
        xxx  = frecuencia en kHz del campo.
        yy   = valor de la corriente Idc (Internal Direct Current) en Ampere. 
        zzz  = valor de muestreo (i.e., cuantos puntos por segundo registradas), en 10**6 muestras por segundo.
        TM   = Tipo de Muestra (FF ferrofluido, FG ferrogel, TT tejido tumoral, TC tejido congelado) 

Almacena el valor de Idc provisto por el nombre del archivo y lo compara con la amplitud de la señal de referencia.
Esperamos ver relacion lineal.
  

"""
#%% Cambiar extension segun sea necesario
'''Texto que identifica los archivos de fondo'''
textofondo = '_fondo.txt' 

'''Texto que identifica los archivos de calibración '''
textocalibracion = '_cal.txt'

'''calibración  de la bobina: constante que dimensionaliza al campo en A/m a partir de la 
calibración  realizada sobre la bobina del RF
'''
#%%
import numpy as np
import matplotlib.pyplot as plt 
import scipy as sc              
from scipy.signal import find_peaks 

''' Apertura de archivos - Cuadro de dialogo para seleccionar archivos:''' 
import tkinter as tk
from tkinter import filedialog
root = tk.Tk()
root.withdraw()
filepaths=tk.filedialog.askopenfilenames(title="Seleccionar archivos con las medidas de la muestra IMPORTANTE: ordenar por fecha de modificación",filetypes=(("Archivos .txt","*.txt"),("Archivos .dat","*.dat"),("Todos los archivos","*.*")))

'''Nombre para el archivo ASCII de salida''' 
nombre_salida = str(input('Elegir un nombre para el archivo de salida:') or 'Prueba')

#%%
nombres_archivos=[]      #lista de str para sacar info 

for item in filepaths:    
    nombres_archivos.append(item.split('/')[-1])

print('Archivos seleccionados: ')
for item in nombres_archivos:
    print(item)

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
     
'''Defino funcion a ajustar con parametros: offset(A), amplitud(B), frecuencia(C), fase(D)'''
def sinusoide(t,A,B,C,D):
    return(A + B*np.sin(2*np.pi*C*t - D))
    
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
#%%
'''Ahora itero [i] y proceso compleamente a cada uno de los archivos seleccionados '''    

delta_t = [] #Base temporal
Idc = []     #Internal direc current en el generador de RF
Amplitudes_m = [] #Lista para almacenar amplitud de señal de referencia para m, f y c.
Amplitudes_f = []
Amplitudes_c = []
FactCal_m = []    #Lista para almacenar el factor de calibracion universal   
FactCal_f = []
FactCal_c = []
Amp_Frec = []
for i in range(0,fa): 
    
    delta_t.append(1e-6/float(filenombres_muestra[i][11:14]))
    Idc.append(float(filenombres_muestra[i][7:9]))
    
    #Importo archivo nomenclado: xxxkHz_yyA_zzzMss_TM.dat
    muestra = np.loadtxt(rutas_de_carga[i][0],skiprows=3,dtype=float)
    t_m = (muestra[:,0]-muestra[0,0])*delta_t[i]    #1er columna 
    v_m = muestra[:,1]*0.001                        #CH1 
    v_r_m = muestra[:,2]*0.001                      #CH2
   
# =============================================================================
#     Agregado para correr cualquier archivo (ya no carga los _fondo y _cal correspondiente)
# =============================================================================
    
    #Imp archivo: xxxkHz_yyA_zzzMss_TM.dat_fondo (fondo: bobinas sin muestra)
    fondo = np.loadtxt(rutas_de_carga[i][1],skiprows=3,dtype=float)
    t_f = (fondo[:,0]-fondo[0,0])*delta_t[i]
    v_f = fondo[:,1]*0.001      
    v_r_f = fondo[:,2]*0.001    

    #Imp: xxxkHz_yyA_zzzMss_TM_cal.dat(calibracion: material paramagnetico, sin histéresis)
    calibracion = np.loadtxt(rutas_de_carga[i][2],skiprows=3,dtype=float)
    t_c = (calibracion[:,0]-calibracion[0,0])*delta_t[i]
    v_c = calibracion[:,1]*0.001      
    v_r_c = calibracion[:,2]*0.001        
    
    #par_x=[promedio de v_r_x, amplitud de v_r_x]
    par_m=[np.mean(v_r_m),(np.max(v_r_m)-np.min(v_r_m))/2] #Muestra
    par_f=[np.mean(v_r_f),(np.max(v_r_f)-np.min(v_r_f))/2] #Fondo
    par_c=[np.mean(v_r_c),(np.max(v_r_c)-np.min(v_r_c))/2] #Calibracion
    
    '''
    Ajusta las 3 referencias con funciones seno: V(t)=V0+A*sin(2*np.pi*f*t - phi)

    Los valores iniciales para los ajustes (seeds):
    par_x[0] = Valor medio de la señal/Offset (i.e. constante aditiva)
    par_X[1] = Amplitud 
    Despues anexo:
    par_X[2] = Frecuencia
    par_X[3] = Fase

    par_x=[offset,amplitud]    
    '''
    #Suavizo las señales de referencia:
    #   fft_smooth(señal a suavizar, nº de armonicos) 
    
    suave_m = fft_smooth(v_r_m, np.around(int(len(v_r_m)*6/1000)))
    suave_f = fft_smooth(v_r_f, np.around(int(len(v_r_f)*6/1000)))
    suave_c = fft_smooth(v_r_c, np.around(int(len(v_r_c)*6/1000)))
       
    '''
    Para evitar problemas por errores de nomenclatura de la frecuencia de los archivos,
    mido el tiempo entre picos, i.e. el periodo. 
    Invierto para obtener la frecuencia semilla y la agrego a la lista par_m  
    ''' 
    indices_m = find_peaks(suave_m,height=0) #tupla
    indices_f = find_peaks(suave_f,height=0)
    indices_c = find_peaks(suave_c,height=0)
    #Acceso al diccionario: indices_m[1]['peak_heights'] 
    
    t_entre_max_m = np.mean(np.diff(t_m[indices_m[0]]))
    t_entre_max_f = np.mean(np.diff(t_f[indices_f[0]]))
    t_entre_max_c = np.mean(np.diff(t_c[indices_c[0]]))
    
    # par_x[2] = Frecuencia Semilla
    par_m.append(1/t_entre_max_m) #Frec de la Muestra
    par_f.append(1/t_entre_max_f) #Frec del Fondo
    par_c.append(1/t_entre_max_c) #Frec de la Calibracion

    #Fase inicial, a partir del tiempo del primer maximo:
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
    #Graficos señales, ajustes y restos
    # Muestra
    n_m = len(t_m)
    y_m = np.empty(n_m)
    for k in range(n_m):
        y_m[k] = sinusoide(t_m[k],coef_m[0],coef_m[1],coef_m[2],coef_m[3])

    #Para calcular el R^2 del ajuste: Metodo a usar: sklearn
    from sklearn.metrics import r2_score 
    R_2_m = r2_score(v_r_m,y_m)
    resto_m = v_r_m - y_m #Resto: diferencia entre funcion original y la ajustada
#%%
#Identifico Offset, Amplitud, Frecuencia y Fase de cada referencia:
    offset_m = coef_m[0]
    amplitud_m = coef_m[1]
    frecuencia_m = coef_m[2]
    frecuencia_m_kHz=coef_m[2]/1000
    fase_m = coef_m[3]
    
    offset_f = coef_f[0]
    amplitud_f = coef_f[1]
    frecuencia_f = coef_f[2]
    fase_f = coef_f[3]

    offset_c = coef_c[0]
    amplitud_c = coef_c[1]
    frecuencia_c = coef_c[2]
    fase_c = coef_c[3]
          
#%%    
    
    plt.plot(t_m, v_r_m,'o',label='Referencia de muestra')
    plt.plot(t_m,y_m,'r-',label='Ajuste de ref. de muestra')
    plt.plot(t_m,resto_m,'.', label='Restos')    
    plt.xlabel('t (s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.text(max(t_m),0,'$R^2$ = {}'.format(R_2_m),
             bbox=dict(alpha=0.9), ha='right',
             va='top')
    plt.text(0,0,'Frecuencia: %0.2f kHz '%frecuencia_m_kHz,
             bbox=dict(alpha=0.9), ha='left',va='top')
    
    
    plt.text(max(t_m),0,'Amplitud: %0.2f V '%amplitud_m,
             bbox=dict(alpha=0.9), ha='right',va='bottom')
    plt.axhspan(amplitud_m, -amplitud_m, facecolor='g', alpha=0.2)
    plt.title(filenombres_muestra[i] + ' - Señal de Muestra, ajuste y restos')
    plt.grid()
    #plt.ylim(-12,12)
    #plt.savefig(filenombres_muestra[i] + ' - ajuste_muestra.png',dpi=500,bbox_inches='tight')
    plt.show()
    #%%
    
    #Fondo
    n_f = len(t_f)
    y_f = np.empty(n_f)
    for j in range(n_f):
        y_f[j] = sinusoide(t_f[j],coef_f[0],coef_f[1],coef_f[2],coef_f[3])

    R_2_f = r2_score(v_r_f,y_f)
    resto_f = v_r_f - y_f

    plt.plot(t_f, v_r_f,'go',label='Referencia de fondo')
    plt.plot(t_f,y_f,'y-',label='Ajuste de ref. de fondo')
    plt.plot(t_f,resto_f,'.', label='Restos')
    plt.xlabel('t (s)')
    plt.ylabel('Amplitud de señal (V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.text(max(t_f),min(v_r_f),'$R^2$ = {}'.format(R_2_f),bbox=dict(alpha=1.0), ha='right',va='bottom')
    plt.title(filenombres_fondo[i] + ' - Señal de Fondo, ajuste y restos')
    #plt.savefig(filenombres_muestra[i] + ' - ajuste_fondo.png',dpi=500,bbox_inches='tight')
    plt.show()

    #Calibracion

    n_c = len(t_c)
    y_c = np.empty(n_c)
    for l in range(n_c):
        y_c[l] = sinusoide(t_c[l],coef_c[0],coef_c[1],coef_c[2],coef_c[3])

    R_2_c = r2_score(v_r_c,y_c)
    resto_c = v_r_c - y_c

    plt.plot(t_c, v_r_c,'o',label='Referencia de calibración ')
    plt.plot(t_c,y_c,'-m', label='Ajuste de ref. de calibración ')
    plt.plot(t_c,resto_c,'.', label='Restos')

    plt.xlabel('t (s)')
    plt.ylabel('Amplitud de señal(V)')
    plt.legend(loc='upper left',framealpha=1.0)
    plt.text(max(t_c),min(v_r_c),'$R^2$ = {}'.format(R_2_c),bbox=dict(alpha=1.0), ha='right',va='bottom')
    plt.title(filenombres_cal[i] + ' -  Señal de calibración , ajuste y restos')
    #plt.savefig(filenombres_cal[i] + ' - ajuste_calibracion.png',dpi=500,bbox_inches='tight')
    plt.show()
    
#%%    
    #Estos arrays contienen los datos que me interesan    
    Frecuencia_muestra_kHz.append(frecuencia_m/1000)
    Amplitudes_m.append(amplitud_m) 
    Amp_Frec.append(amplitud_m/(frecuencia_m/1000))
    Amplitudes_f.append(amplitud_f)
    Amplitudes_c.append(amplitud_c)
      
    FactCal_m.append(amplitud_m/(Idc[i]*frecuencia_m))
    FactCal_f.append(amplitud_f/(Idc[i]*frecuencia_f))
    FactCal_c.append(amplitud_c/(Idc[i]*frecuencia_c))
    FactCal_total= FactCal_m+FactCal_f+FactCal_c
    
    FactCal_total= FactCal_m
    
ValorMedio_m=np.mean(FactCal_m)
VM_m_norm = ValorMedio_m/max(FactCal_m)
Varianza_m=np.std(FactCal_m)    
ValorMedio_f=np.mean(FactCal_f)
VM_f_norm = ValorMedio_f/max(FactCal_f)
Varianza_f=np.std(FactCal_f)      
ValorMedio_c=np.mean(FactCal_c)
VM_c_norm = ValorMedio_c/max(FactCal_c)
Varianza_c=np.std(FactCal_c)  

ValorMedio_total=np.mean(FactCal_total)  
VM_t_norm = ValorMedio_total/max(FactCal_total) 

#%%
labels_frec = []
labels_amp = []
for item in Frecuencia_muestra_kHz:
    labels_frec.append('%.2f kHz' %item)
for item in Amplitudes_m:
    labels_amp.append('%.2f V' %item)

#%% SORTED solo para cuando agrego por nombre, si es por fecha no va

#.fig=plt.figure(figsize=(10,10))
plt.scatter(Idc,Amp_Frec,label='5 Feb 2021')  

plt.xlim(0,16.5)
#plt.ylim(0,16)
#plt.text(15,2,nombre_salida,bbox=dict(alpha=0.8), ha='center',va='top')
plt.xlabel('Idc (A)')
plt.ylabel('Amplitud de señal/Frecuencia ')
for i in range(0,len(Idc)):
        plt.annotate(labels_frec[i],(Idc[i],Amp_Frec[i]-0.005),size=8,
                     bbox=dict(boxstyle="round", alpha=0.5, lw=0.2),ha='center', va='bottom')
plt.legend(loc='upper left',framealpha=1.0)
plt.title(nombre_salida + ' - Amplitud de Señal de referencia vs. Idc')
plt.grid()
#plt.savefig(nombre_salida + ' - Idc vs Amplitud de señal Frecuencia.png',dpi=500,bbox_inches='tight')
#plt.show()

#%%
'''
XX = sorted(Amplitudes_m + Amplitudes_c + Amplitudes_f) #Junta y ordena los 3 arrays
YY = sorted(Idc + Idc + Idc)

#Ajusto recta 
recta = np.polyfit(YY,XX,1)
pendiente = recta[0]
ordenada = recta[1]
t_recta = np.linspace(2,15,1000)
y_recta = pendiente*t_recta+ordenada

fig=plt.figure(figsize=(8,8))
plt.plot(YY,XX,'-o',color='#ff7f0e',label=nombre_salida)
plt.plot(t_recta,y_recta,label='Ajuste lineal')
plt.xlim(0,16)
plt.ylim(0,16)
plt.text(15,2,'%s \nPendiente: %f\nOrdenada: %f' %(nombre_salida,pendiente, ordenada) ,bbox=dict(alpha=0.8), ha='right',va='bottom')
plt.xlabel('Idc (A)')
plt.ylabel('Amplitud de señal (V)')
plt.legend(loc='upper left',framealpha=1.0)
plt.title('Amplitud de Señal de referencia vs. Idc (con ajuste)')
plt.grid()
#plt.savefig(nombre_salida + ' - Idc vs Amplitud de señal (con ajuste).png',dpi=500,bbox_inches='tight')
plt.show()
'''
#%%
import os.path, datetime
fecha_m = []
FMT = '%y-%m-%d %H:%M'
for item in filepaths:
    fecha_m.append(datetime.datetime.fromtimestamp(os.path.getmtime(item)).strftime(FMT))

fecha_min = []    #para horario en minutos para graficar por dia
for date in fecha_m:
    fecha_min.append(date[9:14])
#%% #para fecha absoluta
fecha_m_abs = []    
for item in filepaths: 
    fecha_m_abs.append((os.path.getmtime(item))) #.strftime(FMT))    
fecha_m_abs = np.diff(fecha_m_abs) #.strftime(FMT) 

#dias = [fecha_m[0],fecha_m[44],fecha_m[123],fecha_m[129],fecha_m[130],fecha_m[142],fecha_m[155],fecha_m[-1]]    
    #%%
    
FactCal_m_norm = FactCal_m/max(FactCal_m)     
    
fig=plt.figure(figsize=(11,5))
plt.plot(fecha_min,FactCal_m_norm,'-o',label='Muestra')
plt.plot(fecha_min,FactCal_f/max(FactCal_f),'-s',label='Fondo')
plt.plot(fecha_min,FactCal_c/max(FactCal_c),'-^',label='Calibración')
plt.tick_params(axis='x', which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
plt.axvspan(fecha_min[0],fecha_min[37],  color='g',alpha=0.4)
plt.axvspan(fecha_min[37],fecha_min[42] ,color='r',alpha=0.4)
plt.axvspan(fecha_min[42],fecha_min[42], color='#8c564b',alpha=0.4)
plt.axvspan(fecha_min[42],fecha_min[55], color='c', alpha=0.4)
plt.axvspan(fecha_min[55],fecha_min[60], color='#7f7f7f',alpha=0.4)
plt.axvspan(fecha_min[60],fecha_min[77], color='orange',alpha=0.4)

plt.title('Factor de calibración normalizado al máximo valor - Marzo 2020',fontsize = 18)
plt.ylabel('u.a.')
#plt.xlabel('Horario')
#plt.text(10,min(FactCal_m),'Valor Medio muestra: %0.3e - Std Muestra: %0.3e\nValor Medio fondo: %0.3e - Std fondo: %0.3e\nValor Medio cal: %0.3e - Std cal: %0.3e' %(ValorMedio_m,Varianza_m,ValorMedio_f,Varianza_f,ValorMedio_c,Varianza_c),bbox=dict(alpha=0.8), ha='left',va='bottom')

plt.gcf().autofmt_xdate()
plt.xlim(fecha_min[0],fecha_min[-1])
plt.grid()
plt.hlines(VM_t_norm,fecha_min[0],fecha_min[-1],color='r',linestyles='dashdot', label='Valor medio: %0.2f' %VM_t_norm)
plt.legend(loc='best',framealpha=0.9)
#plt.ylim(0.8,1.01)


plt.savefig(nombre_salida + ' - Factor de calibracion.png',dpi=200,bbox_inches='tight')

plt.show()



