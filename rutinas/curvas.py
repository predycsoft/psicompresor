import time
from sympy import *
from prettytable import PrettyTable
import numpy as np
from simulacionteorica import punto
from adimensional import ADIM

def curvaoperacion(mapa, NRPMS, SURGE, STONEW, z, TSUC, PSUC, DEQ, CC1, CC2, CC3, CC4, EXPOCP, CE1,CE2,CE3, CE4, EXPOEP, GC, step, NRPM, T_DIM, P_DIM, DEQ_DIM, FLUJO_DIM):
    if mapa == True:
        for i in range(len(NRPMS)):
            start = time.time()
            NRPM = NRPMS[i]
            curvas(mapa, SURGE, STONEW,z, TSUC, PSUC, DEQ, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, step, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM)
            end = time.time()
            print("tiempo de grafica de curva", end - start)
    else:
        curvas(mapa, SURGE, STONEW,z, TSUC, PSUC, DEQ, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, step, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM)

def calculosurge():
    return

def calculostw():
    return

def curvas(mapa, SURGE, STONEW,z, TSUC, PSUC, DEQ, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, step, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM):
    rangomin = SURGE*NRPM
    rangomax = STONEW*NRPM
    npunto = 100
    step = (rangomax - rangomin)/(npunto - 1)
    if mapa == True:
        table = PrettyTable(['Caudal [ACFM]','Presion Descarga', 'Temperatura Descarga', 'Q/N', 'Potencia'])
    else:
        table = PrettyTable(['Caudal', 'Q/N', 'Eficiencia', 'C.Head', 'C. Work In', 'Head'])
    for adan in range(0,npunto):
        start = time.time()
        #punto
        caudalreal = rangomin + (adan)*step
        QSBRN, TDES, PDES = punto(z, TSUC, PSUC, caudalreal, DEQ, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, SURGE, STONEW, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM)
        NP, SI, MI, WP, gashp = ADIM(z, DEQ, TSUC, PSUC, TDES, PDES, caudalreal, NRPM, 1, DEQ_DIM, T_DIM, T_DIM, P_DIM, P_DIM, FLUJO_DIM)
        if mapa == True:
            table.add_row([caudalreal, PDES, TDES, caudalreal/NRPM,gashp])
        else:
            table.add_row([caudalreal, QSBRN, NP, MI, SI, WP])
        end = time.time()
        print("tiempo de analisis de punto de curva", end - start)
    print(table)
    return

if __name__ == "__main__":
    curvaoperacion()
