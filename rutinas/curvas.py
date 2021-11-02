import time
import numpy as np
from rutinas.adimensional import *
from rutinas.simulacionteorica import *


def curvas(mapa, SURGE, STONEW,z, TSUC, PSUC, DEQ, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, step, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM):
    #dataset = ["Q/N","Temperatura Descarga", "Presion Descarga","Potencia", "Caudal [ACFM]", "Eficiencia", "C. head", "C. Work In", "Head"]
    dataset = []
    rangomin = SURGE*NRPM
    rangomax = STONEW*NRPM
    npunto = 10
    step = (rangomax - rangomin)/(npunto - 1)
    # if mapa == True:
    #     table = PrettyTable(['Caudal [ACFM]','Presion Descarga', 'Temperatura Descarga', 'Q/N', 'Potencia'])
    # else:
    #     table = PrettyTable(['Caudal', 'Q/N', 'Eficiencia', 'C.Head', 'C. Work In', 'Head'])
    for adan in range(0,npunto):
        start = time.time()
        #punto
        caudalreal = rangomin + (adan)*step
        PDES, TDES, DG, HG, QN = punto(z, TSUC, PSUC, caudalreal, DEQ, NRPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4,EXPOEP, GC, SURGE, STONEW, DEQ_DIM, T_DIM, P_DIM, FLUJO_DIM)[0:5]
        NP, SI, MI, WP, gashp = ADIM(z, DEQ, TSUC, PSUC, TDES, PDES, caudalreal, NRPM, 1, DEQ_DIM, T_DIM, T_DIM, P_DIM, P_DIM, FLUJO_DIM)[0:5]
        if (adan == 0):
            dataset = np.array([NRPM, QN, TDES, PDES, gashp, caudalreal, NP, MI, SI, WP])
        else:
            dataset = np.c_[dataset, [NRPM, QN, TDES, PDES, gashp, caudalreal, NP, MI, SI, WP]]
        end = time.time()
        print("tiempo de analisis de punto de curva", end - start)
        # if mapa == True:
        #     table.add_row([caudalreal, PDES, TDES, caudalreal/NRPM,gashp])
        # else:
        #     table.add_row([caudalreal, QSBRN, NP, MI, SI, WP])
    # print(table)
    return dataset

