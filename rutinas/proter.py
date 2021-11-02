from sympy import exp, ones, ln, log
import math
from decimal import *

from sympy.matrices.expressions.transpose import transpose
from rutinas.constantes import CMW, TC, cd, B, CKIJ, acf, PC, CI
import time

flag = 0
R = 10.7335

"""
Entradas:
  z = Mezcla del gas
  CMW: Pesos Moleculares
  TC: Temperaturas Críticas
  cd: Densidades Críticas
  acf: Factores Acéntricos
  CKIJ: Parametros de interacción entro componentes "K(i,j)" de la ecuación de estado
  CI: Coeficientes polinómicos para el cálculo de entalpía y entropía
  B: Valores de los coeficientes "A(j)" y "B(j)" de la ecuacion de estado
  T: Temperatura
  P: Presión
  R: Constante
  CONVR: Factor de conversión en entalpia y entropia de BTU/lb a psia.ft3/lbmol

Calculos:
  dc: Densidad Calculada en [lbmol/pie3]
  hm: Entalía Calculada
  sm: Entropía Calculada
  dg: Densidad Calculada en [lbmol/pie3]
  ymw: Peso Molecular de la mezcla
  B0, A0, C0, D0, E0, gamma, a, b,c, d, alpha: Parametros de la ecuación BWRS

  Auxiliares de funciones:
    esw1: es el calculo del tanteo que debe ser cero para conseguir la densidad en la ecuación BWRS
    punto1: es el primer punto evaluado del tanteo de debsidades para usar el metodo Régula-Falsi
    punto2: es el segundo punto evaluado del tanteo de debsidades para usar el metodo Régula-Falsi
    DeltaH: es el factor de desviación de la entalpía
    Hideal: es la entalpía calculada para gas ideal
    Deltas: es el factor de desviación de la entropía
    Hideal: es la entropía calculada para gas ideal

Funciones:
  cambiosgas: elimina de las matrices de coeficientes las contantes de los componentes que no estan presentes en la mezcla
  bwrc: calcula los parametros de la ecuacion BWRS (Benedic-Webb-Rubin-Starling)
  densidad: calcula la densidad a través de la ecuación BWRS (Benedic-Webb-Rubin-Starling) usando el método Régula-Falsi
  esw: es la evaluación del tanteo para hallar la densidad que resuelve la ecuacion BWRS
  entalpia: calcula entalpía mediante la ecuación del Anexo 4
  entropia: calcula entropía mediante la ecuación del Anexo 4
  pesomolecular: Calcula el peso molecular de la mezcla y realiza la transformacion de unidades de la densidad

Esta rutina retorna los valores de:
- dc: Densidad en [lbm/lbmol]
- hm: Entalpia
- sm: Entropia
- dg: Densidad en [lbm/pie3]
- ymw: Peso Molecular de la mezcla
"""
AUXCMW = CMW
AUXTC = TC
AUXcd = cd
AUXPC = PC
AUXacf = acf
AUXCKIJ = CKIJ
AUXCI = CI
firstrun = 0
B0 = 0
A0 = 0
C0 = 0
D0 = 0
E0 = 0
gamma = 0
a = 0 
b = 0
c = 0
d = 0
alpha = 0

def PROTER(z,T,P):
    global firstrun
    global AUXCMW, AUXTC, AUXcd, AUXPC, AUXacf, AUXCKIJ, AUXCI
    CMW = AUXCMW; TC = AUXTC; cd = AUXcd; acf = AUXacf; PC = AUXPC; CKIJ = AUXCKIJ; CI = AUXCI
    z, NC, cd, acf, TC, CKIJ, CMW, CI = cambiogas(z, CMW, TC, cd, acf, PC, CKIJ, CI)
    if firstrun == 0:
        global B0, A0, C0, D0, E0, gamma, a, b,c, d, alpha
        B0, A0, C0, D0, E0, gamma, a, b,c, d, alpha = bwrc(z,B,acf,cd,R,TC,CKIJ, NC)
        firstrun = 1
    dc = densidad(z,acf,B0,A0,C0,D0,E0,gamma,a,b,alpha,c,d,R,T,P)
    hm = entalpia(B0, A0, C0, D0, E0,a,b,c,d,alpha,gamma,R,dc,T, CI, CMW, z)
    sm = entropia(B0, A0, C0, D0, E0,a,b,c,d,alpha,gamma,R,dc,T, CI, CMW, z)
    dg, ymw = pesomolecular(CMW,z,dc)
    return dc, hm, sm, dg, ymw

def cambiogas(z, CMW, TC, cd, acf, PC, CKIJ, CI):
    # Reviso que gases estan presentes en la mezcla
    z = transpose(z)
    m = z.shape[0]
    rows = [i for i in range(m) if any(z[i, j] != 0 for j in range(1))]
    # Simplifico Z
    z = z[rows, 0]
    # Numero de componentes
    NC = z.shape[0]
    ## Ahora simplifico los pesos moleculares
    CMW = CMW[rows, 0]
    # Simplifico las temperaturas criticas y los llevo a rankine?
    TC = TC[rows,0] + 459.69*ones(NC,1)
    # Simplificar densidades criticas
    cd = cd[rows, 0]
    # Simplificar factores acéntricos
    acf = acf[rows, 0]
    # Simplificar presiones criticas
    PC = PC[rows,0]
    # Coeficientes entalpia y entropia
    CI = CI[rows,0:7]
    # Simplificacion y completación de Coeficientes de interaccion
    CKIJ = CKIJ[rows, rows]
    for i in range(NC-1):
        for j in range(i+1, NC):
            CKIJ[j,i] = CKIJ[i,j]
    # Calculo de Peso Molecular
    PMOLEC = z.dot(CMW)
    return z, NC, cd, acf, TC, CKIJ, CMW, CI

def bwrc(z,B,acf,cd,R,TC,CKIJ, NC):
    star = time.time()
    B0 = 0; A0 = 0; C0 = 0; D0 = 0; E0 = 0; gamma = 0; a = 0; b = 0;c = 0; d = 0; alpha = 0
    for i in range(0,NC):
        gamma = gamma + (z[i,0]*((B[0,3] + B[1,3]*acf[i,0])/cd[i,0]**2)**(1/2))
        a = a+ (z[i,0]*((B[0,5]+B[1,5]*acf[i,0])*R*TC[i,0]/cd[i,0]**2)**(1/3))
        b = b + (z[i,0]*((B[0,4]+B[1,4]*acf[i,0])/cd[i,0]**2)**(1/3))
        alpha = alpha + (z[i,0]*((B[0,6]+B[1,6]*acf[i,0])/cd[i,0]**3)**(1/3))
        c = c + (z[i,0]*((B[0,7]+B[1,7]*acf[i,0])*R*TC[i,0]**3/cd[i,0]**2)**(1/3))
        d = d + (z[i,0]*((B[0,9]+B[1,9]*acf[i,0])*R*TC[i,0]**2/cd[i,0]**2)**(1/3))
        B0 = B0 + (z[i,0]*(B[0,0]+B[1,0]*acf[i,0])/cd[i,0])
        for j in range(0,NC):
            A0 = A0 + (z[i,0]*z[j,0]*((B[0,1]+B[1,1]*acf[i,0])*R*TC[i,0]/cd[i,0])**(1/2)*((B[0,1]+B[1,1]*acf[j,0])*R*TC[j,0]/cd[j,0])**(1/2)*(1-CKIJ[i,j]))
            C0 = C0 + (z[i,0]*z[j,0]*((B[0,2]+B[1,2]*acf[i,0])*R*TC[i,0]**3/cd[i,0])**(1/2)*((B[0,2]+B[1,2]*acf[j,0])*R*TC[j,0]**3/cd[j,0])**(1/2)*(1-CKIJ[i,j])**3)
            D0 = D0 +(z[i,0]*z[j,0]*((B[0,8]+B[1,8]*acf[i,0])*R*TC[i,0]**4/cd[i,0])**(1/2)*((B[0,8]+B[1,8]*acf[j,0])*R*TC[j,0]**4/cd[j,0])**(1/2)*(1-CKIJ[i,j])**4)
            E0 = E0 + (z[i,0]*z[j,0]*((B[0,10]+B[1,10]*acf[i,0]*exp(-3.8*acf[i,0]))*R*TC[i,0]**5/cd[i,0])**(1/2)*((B[0,10]+B[1,10]*acf[j,0]*exp(-3.8*acf[j,0]))*R*TC[j,0]**5/cd[j,0])**(1/2)*(1-CKIJ[i,j])**5)
    gamma = gamma**2
    a = a**3
    b = b**3
    alpha = alpha**3
    c = c**3
    d = d**3
    end = time.time()
    print("tiempo proter", end - star)
    return B0, A0, C0, D0, E0, gamma, a, b,c, d, alpha

def densidad(z,acf,B0,A0,C0,D0,E0,gamma,a,b,alpha,c,d,R,T,P):
    acfm = z.dot(acf)
    # calculo de coeficientes para calculo de densidad critica
    C1 = (R*T)
    C2 = (B0 * R * T - A0 - C0 / T **2 + D0 / T**3 - E0 / T**4)
    C3 = (b * R * T - a - d / T)
    C6 = (alpha * (a + d / T))
    # Asignaciones iniciales del tanteo
    global flag
    if flag == 0:
        D1 = 0
        D2 = D1 + 0.00005
    else:
        D1 = 0.01
        D2 = 0.001
    DC = D1
    flag = flag + 1
    esw1 = esw(C1, C2, C3, C6, c, gamma, P, DC, T)
    punto1 = esw1
    DC = D2
    esw1 = esw(C1, C2, C3, C6, c, gamma, P, DC, T)
    punto2 = esw1
    # Inicia el tanteo de busqueda de densidad critica
    while 1:
        # Calculo de intercepcion
        DC = ((D1*punto2 - D2*punto1)/(punto2 - punto1))
        esw1 = esw(C1, C2, C3, C6, c, gamma, P, DC, T)
        # chqueo de convergencia SE DISMINUYO LA TOLERACIA DE CONVERGENCIA PONER EN OBSERVACION EN FUNCION AL TIEMPO 0.0001
        if abs(esw1) <= 0.0001:
            return DC
        # Nuevas asignaciones por la no convergencia
        if abs(D1-DC) >= abs(D2 -DC):
            D1 = DC
            punto1 = esw1
        else:
            D2 = DC
            punto2 = esw1

def esw(C1, C2, C3, C6, c, gamma, P, DC, T):
    #print("C1:", C1, "C2", C2,"C3",C3,"C6",C6, "c", c,"gamma",gamma,"P",P,"DC",DC, "T",T)
    esw1 = C1 * DC + C2 * (DC**2) + C3 * (DC**3) + C6 * (DC**6 )+  c * (DC**3) / (T**2) * (1 + gamma * (DC**2)) * math.exp(-gamma * DC**2) - P
    #print("esw1:", esw1, "DC", DC)
    return esw1

def entalpia(B0, A0, C0, D0, E0,a,b,c,d,alpha,gamma,R,DC,T, CI, CMW, z):
    CONVR = 0.185057
    DeltaH1 = (B0*R*T - 2*A0 - 4*C0/(T**2) + 5*D0/(T**3) - 6*E0/(T**4)) * DC
    DeltaH2 = (1/2)*(2*b*R*T -3*a -4*d/T) * (DC**2)
    DeltaH3 = (1/5)*alpha*(6*a + 7*d/T) * (DC**5)
    DeltaH4 = (c/(gamma*(T**2))) * (3 - (3+ (1/2)*gamma*(DC**2) - (gamma**2)*(DC**4)) * math.exp(-gamma*(DC**2)))
    DeltaH = (DeltaH1 + DeltaH2 + DeltaH3 + DeltaH4)*CONVR
    CI = CI[0:CI.shape[0],0:6]
    Hideal = 0
    for i in range(0, CI.shape[0]):
        Hideal = Hideal + CMW[i,0]*z[i,0]*(CI[i,0]+CI[i,1]*T+CI[i,2]*T**2+CI[i,3]*T**3+CI[i,4]*T**4+CI[i,5]*T**5)
    hc = DeltaH + Hideal
    return hc

def entropia(B0, A0, C0, D0, E0,a,b,c,d,alpha,gamma,R,DC,T, CI, CMW, z):
    CONVR = 0.185057
    # Pregunta por qué es 14.7?
    DeltaS1 = -R * ln(DC * R * T /14.7)
    DeltaS2 = -(B0*R + 2*C0/(T**3) - 3*D0/ (T**4) + 4*E0/(T**5)) * DC
    DeltaS3 = -(1/2) * (b*R + d/(T**2)) * (DC**2)
    DeltaS4 = alpha*d*(DC**5)/(T**5)
    DeltaS5 = ( (2*c) / (gamma*(T**3)) )  *( 1- ( 1 + (1/2)*gamma*(DC**2))* math.exp(-gamma*(DC**2)))
    DeltaS  = (DeltaS1 + DeltaS2 + DeltaS3 + DeltaS4 + DeltaS5)*CONVR
    XLNX = 0
    Sideal = 0
    for i in range(0, CI.shape[0]):
        Sideal = Sideal + CMW[i,0]*z[i,0]*(CI[i,1]*log(T)+2*CI[i,2]*(T-1)+1.5*CI[i,3]*(T**2-1)+(1.3333333333)*CI[i,4]*(T**3-1)+ 1.25* CI[i,5]*(T**4-1)+CI[i,6])
        XLNX = XLNX + z[i,0] * log(z[i,0])
    Sideal = Sideal - R * XLNX * CONVR
    sc = DeltaS + Sideal
    return sc

def pesomolecular(CMW, z, DG):
    ymw = CMW.dot(z).evalf(20)
    DG = ymw*DG
    return DG, ymw

