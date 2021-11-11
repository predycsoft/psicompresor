from rutinas.proter import *
from sympy import *
from prettytable import PrettyTable
import numpy as np
"""
Entradas:

  Constantes:
  GC: No se que es

  Variables:
  NSEC: numero de secciones
  FLUJO: Flujo de la operacion
  TAU: Relacion de transmision
  NRPM: Velocidad [RPM]
  CAIPRES: Caida de presion
  gases: arreglo de cromatografía en cada seccion
  diametros: arreglo de diametro de cada seccion
  surge: arreglo de limite de surge de cada seccion
  stonew = arreglo de limite de stone wall de cada seccion
  expocp = arreglo de exponente de coeficiente politropico de cada seccion (polinomio)
  cc1 = arreglo de c1 de coeficiente politropico de cada seccion (polinomio)
  cc2 = arreglo de c2 de coeficiente politropico de cada seccion (polinomio)
  cc3 = arreglo de c3 de coeficiente politropico de cada seccion (polinomio)
  cc4 = arreglo de c4 de coeficiente politropico de cada seccion (polinomio)
  expoep = arreglo de exponente de eficiencia politropica de cada seccion (polinomio)
  ce1 = arreglo de c1 de eficiencia politropica de cada seccion (polinomio)
  ce2 = arreglo de c2 de eficiencia politropica de cada seccion (polinomio)
  ce3 = arreglo de c3 de eficiencia politropica de cada seccion (polinomio)
  ce4 = arreglo de c4 de eficiencia politropica de cada seccion (polinomio)
  tsuc = arreglo de temperatura de succion de cada seccion
  psuc = arreglo de presion de succion de cada seccion
  divflj = arreglo de extraccion de gas de cada seccion
  relvel = arreglo de relacion de velocidad c/r al tren de cada seccion
  nimpuls = arreglo de numero de impulsores de cada seccion
  # Datos del fabricante
  DIAM = Diametro impulsor de seccion a evaluar
  SURGE = Limite de surge de seccion a evaluar
  STONEW = Limite de stone wall de seccion a evaluar
  EXPOCP = Exponente de coeficiente politropico de seccion a evaluar
  CC1 = C1 de coeficiente politropico de seccion a evaluar
  CC2 = C2 de coeficiente politropico de seccion a evaluar
  CC3 = C3 de coeficiente politropico de seccion a evaluar
  CC4 = C4 de coeficiente politropico de seccion a evaluar
  EXPOEP = Exponente de eficiencia politropica de seccion a evaluar
  CE1 = C1 de eficiencia politropica de seccion a evaluar
  CE2 = C2 de eficiencia politropica de seccion a evaluar
  CE3 = C3 de eficiencia politropica de seccion a evaluar
  CE4 = C4 de eficiencia politropica de seccion a evaluar
  z = mezcla de gas
  TSUC = temperatura de succion
  PSUC = presion de succion
  DIVFLJ = extraccion de gas
  RELVEL = relacion de velocidad
  NIMPULS = numero de impulsores

Calculos:
  CFHEAD: coeficiente politropico
  EFIC: eficiencia politropica
  QN: Q/N
  n: exponente politropico
  dc: Densidad Calculada en [lbmol/pie3]
  hm: Entalía Calculada
  sm: Entropía Calculada
  DG: Densidad Calculada en [lbmol/pie3]
  ymw: Peso Molecular de la mezcla
  deltaHComp: cambio de entalpia de gas
  deltaHGas: cambio de entalpia de gas calculada con exponente politropico estimado
  delta: diferencia entre cambios de entalpia para determinar si se consiguio la presion de succion correcta
  XTAN, YTAN: variables de tanteo
  P2: presion de descarga
  T2: Temperatura de descarga

  Auxiliares de funciones:
    XTAN1 XTAN2 YTAN1 YTAN2: Variables para el tanteo isentrópico

Funciones:
  punto
  temperatura
  presion
  polinomio
  simulacion
  cambioentalpiacomp
  cambioentalpiagas
  sumariocompresion

Esta rutina retorna los valores de:
    'Presion[psia]': P1
    'Temperatura[°F]:' T1
    'Densidad[lbm/ft3]: DG1
    'Entalpia[btu/lbmol]: HG
    'Q/N Surge:': Surge
    'Q/N:':  QN])
    'Q/N stw':  STONEW])
    'Coef. Head':  CFHEAD])
    'Head:':  HEAD])
    'Eficiencia Politrópica:':  EFIC])
    'HP gas:':  HP])
    'Exp. polit:':  POLLY])
    'FLUJO[MMSCFD]':  CAUDAL/2634.5974])
    'Velocidad[RPM]':  NRPM])
    'Total HP Gas':  tothpsec])
    'Temperatura de Succión [°F]':  t0 - 459.69])
    'Temperatura de Descarga [°F]':  tdesc])
    'Presion de Succion [psig]':  p0 - 14.7])
    'Presion de Descarga [psig]': pdesc])
"""

def punto(z, TSUC, PSUC, FLUJO, DIAM, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, SURGE, STONEW, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM):
    # Cambios de dimensiones
    #DIAMETRO
    if (DEQ_DIM == '[pulg]'):
        DIAM = DIAM
    elif (DEQ_DIM == '[mm]'):
        DIAM = DIAM / 25.4
    elif (DEQ_DIM == '[pies]'):
        DIAM = DIAM * 12
    #Temperatura en la succión.
    if (T_DIM== '[°F]'):
        TSUC = TSUC
    elif (T_DIM== '[°K]'):
        TSUC = TSUC * 1.8 - 459.67
    elif (T_DIM== '[°C]'):
        TSUC = TSUC * 1.8 + 32
    #Presión en la succión.
    if (P_DIM == '[psia]'):
        PSUC = PSUC - 14.7
    elif (P_DIM == '[psig]'):
        PSUC = PSUC
    elif (P_DIM == '[KPag]'):
        PSUC = PSUC * 0.14503773773
    elif (P_DIM == '[barg]'):
        PSUC = PSUC * 14.503773773
    ## Calculo de propiedades
    P1 = PSUC + 14.7
    T1 = TSUC + 459.69
    dc1, HG, sg1, DG, ymw = PROTER(z, T1, P1)
    vg = 1/DG
    # Cambiar unidades de flujo
    if FLUJO_DIM == "[lbm/min]":
        FLUJO = FLUJO *1440/ymw
    elif  FLUJO_DIM == "[MMSCFD]":
        FLUJO = FLUJO * 2634.5974
    elif  FLUJO_DIM == "[lbmol/dia]":
        FLUJO = FLUJO
    elif  FLUJO_DIM == "[ACFM]":
        FLUJO = FLUJO*dc1*1440
    elif  FLUJO_DIM == "[m3/min]":
        FLUJO = FLUJO*dc1*1440*35.31466672
    CFHEAD, QN = polinomio(FLUJO, ymw, vg, NRPM, CC1, CC2, CC3, CC4, EXPOCP)
    EFIC, QN = polinomio(FLUJO, ymw, vg, NRPM, CE1, CE2, CE3, CE4, EXPOEP)
    if QN < 0 and QN > 100000:
        print("En el impulsor Nº el flujo ha sido modificado ya que sale del rango recomendado por el fabricante")
        return
    # valor inicial del exponente politrópico "n" se toma 1.34868
    n = 1.34868
    # Se calcula la presión de descarga
    PDES = presion(P1, CFHEAD, n, DIAM, NRPM, vg, GC)
    TDES = temperatura(T1, P1, PDES, n)
    # Se calucla el cambio de entalpia del gas
    deltaHComp = cambioentalpiacomp(P1, PDES, vg, EFIC, n)
    deltaHGas, hg2 = cambioentalpiagas(z, TDES, PDES, HG)
    # Se verifica si era correcto
    delta = deltaHComp - deltaHGas
    ## Inicia Tanteo
    XTAN1 = n
    YTAN1 = delta
    XTAN2 = n - 0.01
    PDES = presion(P1, CFHEAD, XTAN2, DIAM, NRPM, vg, GC)
    TDES = temperatura(T1, P1, PDES, XTAN2)
    deltaHComp = cambioentalpiacomp(P1, PDES, vg, EFIC, n)
    deltaHGas, hg2 = cambioentalpiagas(z, TDES, PDES, HG)
    delta = deltaHComp - deltaHGas
    YTAN2 = delta
    if abs(delta) > 0.06:
        while 1:
            n = (XTAN1 * YTAN2 - XTAN2 * YTAN1)/( YTAN2 - YTAN1)
            PDES = presion(P1, CFHEAD, n, DIAM, NRPM, vg, GC)
            TDES = temperatura(T1, P1, PDES, n)
            deltaHComp = cambioentalpiacomp(P1, PDES, vg, EFIC, n)
            deltaHGas, hg2 = cambioentalpiagas(z, TDES, PDES, HG)
            delta = deltaHComp - deltaHGas
            if (abs(n - XTAN1) >= abs(n - XTAN2)):
                XTAN1 = n
                YTAN1 = delta
            else:
                XTAN2 = n
                YTAN2 = delta
            print("...convergiendo entalpia teorica", delta)
            if abs(delta) <= 0.006:
                break
    PDES = presion(P1, CFHEAD, n, DIAM, NRPM, vg, GC)
    TDES = temperatura(T1, P1, PDES, n)
    PDES = PDES - 14.7
    TDES = TDES - 459.69
    HP = deltaHComp*ymw*FLUJO/61080
    tothpsec = HP
    POLLY = n
    PROD = (3.14159*DIAM*NRPM)**2
    HEAD = CFHEAD * PROD/518400/GC
    sumario_compresion(P1, T1, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, FLUJO, NRPM, tothpsec, TDES, PDES)
    return PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY

def sumario_compresion(p0, t0, DG, hg, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, CAUDAL, NRPM, tothpsec, tdesc, pdesc):
    table = PrettyTable(['Propiedad', 'Valor'])
    table.add_row(['Presion[psia]', p0])
    table.add_row(['Temperatura[°F]:', t0 - 459.69])
    table.add_row(['Densidad[lbm/ft3]:', DG])
    table.add_row(['Entalpia[btu/lbmol]:', hg])
    table.add_row(['Q/N Surge:', SURGE])
    table.add_row(['Q/N:', QN])
    table.add_row(['Q/N stw', STONEW])
    table.add_row(['Coef. Head', CFHEAD])
    table.add_row(['Head:', HEAD])
    table.add_row(['Eficiencia Politrópica:', EFIC])
    table.add_row(['HP gas:', HP])
    table.add_row(['Exp. polit:', POLLY])
    table.add_row(['FLUJO[MMSCFD]', CAUDAL/2634.5974])
    table.add_row(['Velocidad[RPM]', NRPM])
    table.add_row(['Total HP Gas', tothpsec])
    table.add_row(['Temperatura de Succión [°F]', t0 - 459.69])
    table.add_row(['Temperatura de Descarga [°F]', tdesc])
    table.add_row(['Presion de Succion [psig]', p0 - 14.7])
    table.add_row(['Presion de Descarga [psig]',pdesc])
    #print(table)

def cambioentalpiacomp(P1, P2, V1, EFIC, n):
    Parte1 = (0.185028*P1*V1)/(EFIC*((n-1)/n))
    Parte2 = ((P2/P1)**((n-1)/n)) - 1
    deltaH = Parte1*Parte2
    return deltaH

def cambioentalpiagas(z, T2, P2, HG):
    dc2, hg2, sg2, DG2, ymw = PROTER(z, T2, P2)
    deltaHGas = (hg2 - HG)/ymw
    return deltaHGas, hg2

def presion(P1, CFHEAD, n, DIAM, NRPM, V1, GC):
    polym = n/(n-1)
    prod = (3.14159 * DIAM * NRPM )**2
    ctep2 = CFHEAD*prod/(P1*V1*GC*74649600)
    P2 = P1*(ctep2/polym + 1)**polym
    return P2

def temperatura(T1, P1, P2, n):
    T2 = T1*(P2/P1)**((n-1)/n)
    return T2

def polinomio(CAUDAL, ymw, vg, NRPM, C1, C2, C3, C4, EXP):
    QACTUA = CAUDAL * ymw*vg/1440
    QN = QACTUA/NRPM
    value = C1 + C2*QN + C3*QN**2 + C4*QN**EXP
    return value, QN

# def simulacion(gases, tsuc, psuc, FLUJO, diametros, NSEC, surge, stonew, expocp, cc1,cc2,cc3,cc4,expoep,ce1,ce2,ce3,ce4,divflj,relvel,nimpuls, NRPM):
#     for i in range(NSEC):
#         # Datos del fabricante
#         DIAM = diametros[i]
#         SURGE = surge[i]
#         STONEW = stonew[i]
#         EXPOCP = expocp[i]
#         CC1 = cc1[i]
#         CC2 = cc2[i]
#         CC3 = cc3[i]
#         CC4 = cc4[i]
#         EXPOEP = expoep[i]
#         CE1 = ce1[i]
#         CE2 = ce2[i]
#         CE3 = ce3[i]
#         CE4 = ce4[i]
#         # Datos del proceso
#         print(gases.shape[1])
#         z = transpose(gases.row(i))
#         #z = Matrix((0.480000,0.066000,0.040100,0.005000,0.018000,0.005000,0.006500,0.007400,0.000000,0.000000,0.000000,0.000000,0.330000,0.021000,0.021000))
#         TSUC = tsuc[i]
#         PSUC = psuc[i]
#         ### Revisar extrac de gas combustible
#         DIVFLJ = divflj[i]
#         RELVEL = relvel[i]
#         NIMPULS = nimpuls[i]
#         punto(z, TSUC, PSUC, FLUJO, DIAM, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, SURGE, STONEW)
    return



