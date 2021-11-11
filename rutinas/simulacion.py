from sympy import Matrix, pprint, transpose
import numpy as np
from constantes import GC

from simulacionteorica import punto
from adimensional import ADIM

def simulacion():
    tiposim = 3
    ## Simulacion teorica
    if tiposim == 1:
        salida = simulacionTeorica()
        for i in range(salida.shape[0]):
            row = []
            for j in range(salida.shape[1]):
                row.append(salida[i][j])
            print(row)
    if tiposim == 2:
        salida = simulacionCampo()
        for i in range(salida.shape[0]):
            row = []
            for j in range(salida.shape[1]):
                row.append(salida[i][j])
            print(row)
    if tiposim == 3:
        pruebaEficiencia()

def pruebaEficiencia():
    salidaCampo = simulacionCampo()
    pdescampo = salidaCampo[15,1:]
    NSEC = pdescampo.shape[0]
    # print(NSEC)
    salida =np.array(["tipo", "ID", "PSUC", "PDES", "TSUC", "TDES", "Densidad", "Entalpia", "SURGE", "QN", "STONEW", "CFHEAD", "HEAD", "EFIC", "HP", "POLLY", "Q", "RPM"])
    salida = np.c_[salida]
    from entradateorica import entradas
    comp, secc, impul, impuleq, gas, diam, rpm, tsuc, psuc, q, dimT, dimQ, dimP, dimD, caipres, relvel, stonew, surge, expocp, cc1, cc2, cc3, cc4, expoep, ce1, ce2, ce3, ce4, qExt, NIMPULS = buildVariablesTeorica(entradas)
    PSUC = psuc[0]
    TSUC = tsuc[0]
    for i in range(NSEC):
        COMP = comp[i]
        SECC = secc[i]
        IMPUL = impul[i]
        z = transpose(gas[0:15,i])
        DIAM = diam[i]
        RPM = rpm[i]
        Q = q[i]
        DIMT = dimT[i]
        DIMQ = dimQ[i]
        DIMP = dimP[i]
        DIMD = dimD[i]
        CAIPRES = caipres[i]
        RELVEL = relvel[i]
        STONEW = stonew[i]
        SURGE = surge[i]
        EXPOCP = expocp[i]
        CC1 = cc1[i]
        CC2 = cc2[i]
        CC3 = cc3[i]
        CC4 = cc4[i]
        EXPOEP = expoep[i]
        CE1 = ce1[i]
        CE2 = ce2[i]
        CE3 = ce3[i]
        CE4 = ce4[i]
        QEXT = qExt[i]
        PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, Q, DIAM, RPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DIMD, DIMT, DIMP, DIMQ)
        XTAN1 = RPM
        YTAN1 = PDES
        XTAN2 = RPM*0.95
        PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, Q, DIAM, XTAN2, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DIMD, DIMT, DIMP, DIMQ)
        YTAN2 = PDES
        while 1:
            RPM = (XTAN1 * (YTAN2 - PDES) - XTAN2 * (YTAN1 - PDES)) / (YTAN2 - YTAN1)
            PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, Q, DIAM, XTAN2, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DIMD, DIMT, DIMP, DIMQ)
             # Asignacion de nuevos valores del tante
            if abs(XTAN1 - RPM) >= abs(XTAN2 - RPM):
                XTAN1 = RPM
                YTAN1 = PDES
            else:
                XTAN2 = RPM
                YTAN2 = PDES
            ## Chequeo convergencia
            if abs(pdescampo[i]-PDES) <= 0.001:
                break
        obj = np.array(["s", i, PSUC, PDES, TSUC, TDES, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, Q, RPM])
        print("entre a sumario")
        salida = np.c_[salida,obj]
    print(salida)

def simulacionTeorica():
    salida =np.array(["tipo", "ID", "PSUC", "PDES", "TSUC", "TDES", "Densidad", "Entalpia", "SURGE", "QN", "STONEW", "CFHEAD", "HEAD", "EFIC", "HP", "POLLY", "Q", "RPM"])
    salida = np.c_[salida]
    from entradateorica import entradas
    comp, secc, impul, impuleq, gas, diam, rpm, tsuc, psuc, q, dimT, dimQ, dimP, dimD, caipres, relvel, stonew, surge, expocp, cc1, cc2, cc3, cc4, expoep, ce1, ce2, ce3, ce4, qExt, NIMPULS = buildVariablesTeorica(entradas)
    PSUC = psuc[0]
    TSUC = tsuc[0]
    HPSEC = 0
    TSUCSEC = TSUC
    PSUCSEC = PSUC
    HPCOMP = 0
    TSUCCOMP = TSUC
    PSUCCOMP = PSUC
    for i in range(NIMPULS):
        COMP = comp[i]
        SECC = secc[i]
        IMPUL = impul[i]
        z = transpose(gas[0:15,i])
        DIAM = diam[i]
        RPM = rpm[i]
        Q = q[i]
        DIMT = dimT[i]
        DIMQ = dimQ[i]
        DIMP = dimP[i]
        DIMD = dimD[i]
        CAIPRES = caipres[i]
        RELVEL = relvel[i]
        STONEW = stonew[i]
        SURGE = surge[i]
        EXPOCP = expocp[i]
        CC1 = cc1[i]
        CC2 = cc2[i]
        CC3 = cc3[i]
        CC4 = cc4[i]
        EXPOEP = expoep[i]
        CE1 = ce1[i]
        CE2 = ce2[i]
        CE3 = ce3[i]
        CE4 = ce4[i]
        QEXT = qExt[i]
        PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, Q, DIAM, RPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DIMD, DIMT, DIMP, DIMQ)
        salida = sumarioimpulsor(salida, IMPUL, SECC, COMP, PSUC, TSUC, PDES, TDES, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, Q, RPM)
        PSUC = PDES
        TSUC = TDES
        HPSEC = HPSEC + HP
        HPCOMP = HPCOMP + HP
        print("if", i, NIMPULS)
        if i < NIMPULS -1:
            print("if", secc[i], secc[i+1])
            if secc[i] != secc[i+1]:
                TDESSEC = TDES
                PDESSEC = PDES
                TDESCOMP = TDESSEC
                PDESCOMP = PDESSEC
                EFICSEC, COEFWORKINP, CFHEADADIM, WORKPOLY, HPADIM = ADIM(z, DIAM, TSUCSEC, PSUCSEC, TDESSEC, PDESSEC, Q, RPM, NIMPULS, DIMD, DIMT, DIMT, DIMP, DIMP, DIMQ)[0:5]
                salida =  sumarioseccion(salida, IMPUL, SECC, COMP, HPSEC, TSUCSEC, TDESSEC, PDESSEC, PSUCSEC, EFICSEC, HPADIM)

                TSUC = tsuc[i+1]
                HPSEC = 0
                TSUCSEC = TSUC
                PSUCSEC = PSUC
                if comp[i] != comp[i+1]:
                    salida =  sumariocompresor(salida, IMPUL, SECC, COMP, TSUCCOMP, TDESCOMP, PSUCCOMP, PDESCOMP, HPCOMP)
                    TSUCCOMP = TSUCSEC
                    PSUCCOMP = PSUCSEC
                    HPCOMP = 0
        else:
            salida =  sumarioimpulsor(salida, IMPUL, SECC, COMP, PSUC, TSUC, PDES, TDES, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, Q, RPM)
            TDESSEC = TDES
            PDESSEC = PDES
            EFICSEC, COEFWORKINP, CFHEADADIM, WORKPOLY, HPADIM = ADIM(z, DIAM, TSUCSEC, PSUCSEC, TDESSEC, PDESSEC, Q, RPM, NIMPULS, DIMD, DIMT, DIMT, DIMP, DIMP, DIMQ)[0:5]
            salida = sumarioseccion(salida, IMPUL, SECC, COMP, HPSEC, TSUCSEC, TDESSEC, PDESSEC, PSUCSEC, EFICSEC, HPADIM)
            TDESCOMP = TDESSEC
            PDESCOMP = TDESSEC
            salida = sumariocompresor(salida, IMPUL, SECC, COMP, TSUCCOMP, TDESCOMP, PSUCCOMP, PDESCOMP, HPCOMP)
            return salida

def simulacionCampo():
    salida = np.array(["Seccion","Potencia Gas","Flujo Másico","Eficiencia Politrópica","Coef. de Flujo Q/N","Coef.Head Politrópico","Coeficiente Work Input","Trabajo Politrópico","Relación de Compresión","Relación de Volumen","Volumen ACFM"," Temp. Succión","Temp. Descarga","Temp. Isentrópica","Presión Succión ","Presión Descarga ","Presión Isentrópica","Densidad Succión","Densidad Descarga","Densidad Isentrópica","Volumen Succión ","Volumen Descarga ","Volumen isentrópico","Entalpía Succión","Entalpía Descarga","Entalpía Isentrópica","Entropía Succión","Entropía Descarga","Entropía Isentrópica","Calidad Succión","Calidad Descarga","Calidad isentrópica","Fac. Comp. Succión","Fac. Comp. Desc.","Fac. Comp. Isent.","Peso molecular"])
    salida = np.c_[salida]
    from entradateorica import entradacampo
    comp, secc, gas, diam, rpm, tsuc, tdes, psuc, pdes, q, dimT, dimQ, dimP, dimD, NSECC = buildVariablesCampo(entradacampo)
    for i in range(NSECC):
        SECC = secc[i]
        z = transpose(gas[0:15,i])
        DIAM = diam[i]
        RPM = rpm[i]
        TSUC = tsuc[0][i]
        TDES = tdes[i]
        PSUC = psuc[i]
        PDES = pdes[i]
        Q = q[i]
        DIMT = dimT[i]
        DIMQ = dimQ[i]
        DIMP = dimP[i]
        DIMD = dimD[i]
        EFIC, COEFWORKIN, CFHEAD, WP, HPGAS, QW, RELCOMP, RELVOL, TISENT, PISENT, DSUC, DDES, DISENT, VOLSUC, VOLDES, VOLISEN, HSUC, HDES, HISENT,SSUC,SDES,SISENT, FACCOMPSUC, FACCOMPDES, FACCOMPISENT, PESOMOL, QN  = ADIM(z, DIAM, TSUC, PSUC, TDES, PDES, Q, RPM, 1, DIMD, DIMT, DIMT, DIMP, DIMP, DIMQ)
        salida = sumariocampo(salida, SECC, HPGAS, QW, EFIC, QN, CFHEAD, COEFWORKIN, WP, RELCOMP, RELVOL, "Volumen", TSUC, TDES, TISENT, PSUC, PDES, PISENT, DSUC, DDES, DISENT, VOLSUC, VOLDES, VOLISEN, HSUC, HDES, HISENT,SSUC,SDES,SISENT, 1, 1, 1, FACCOMPSUC, FACCOMPDES, FACCOMPISENT, PESOMOL)
    return salida


def sumariocampo(salida, SECC, HPGAS, QW, EFIC, QN, CFHEAD, COEFWORKIN, WP, RELCOMP, RELVOL, VOL, TSUC, TDES, TISENT, PSUC, PDES, PISENT, DSUC, DDES, DISENT, VOLSUC, VOLDES, VOLISEN, HSUC, HDES, HISENT,  SSUC, SDES, SISENT, CALSUC, CALDES, CALISENT, FACCOMPSUC, FACCOMPDES, FACCOMPISENT, PESOMOL):
    obj = np.array([SECC, HPGAS , QW , EFIC , QN , CFHEAD, COEFWORKIN, WP, RELCOMP, RELVOL, VOL, TSUC, TDES, TISENT, PSUC, PDES, PISENT, DSUC, DDES, DISENT, VOLSUC, VOLDES, VOLISEN, HSUC, HDES, HISENT, SSUC, SDES, SISENT, CALSUC, CALDES, CALISENT, FACCOMPSUC, FACCOMPDES, FACCOMPISENT, PESOMOL])
    salida = np.c_[salida,obj]
    return salida

def sumarioimpulsor(salida, IMPUL, SECC, COMP, PSUC, TSUC, PDES, TDES, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, Q, RPM):
    tipo = "I"
    ID = str(IMPUL) + str(SECC) + str(COMP)
    PSUC = PSUC
    PDES = PDES
    TSUC = TSUC
    TDES = TDES
    DG = DG
    HG = HG
    SURGE = SURGE
    QN = QN
    STONEW = STONEW
    CFHEAD = CFHEAD
    HEAD = HEAD
    EFIC = EFIC
    HP = HP
    POLLY = POLLY
    Q = Q
    RPM = RPM
    obj = np.array([tipo, ID, PSUC, PDES, TSUC, TDES, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, Q, RPM])
    print("entre a sumario")
    salida = np.c_[salida,obj]
    return salida

def sumarioseccion(salida, IMPUL, SECC, COMP, HPSEC, TSUCSEC, TDESSEC, PDESSEC, PSUCSEC, EFICSEC, HPADIM):
    tipo = "S"
    ID = str(IMPUL) + str(SECC) + str(COMP)
    PSUC = ""
    PDES = ""
    TSUC = ""
    TDES = ""
    DG = ""
    HG = ""
    SURGE = ""
    QN = ""
    STONEW = ""
    CFHEAD = ""
    HEAD = ""
    EFIC = ""
    HP = ""
    POLLY = ""
    Q = ""
    RPM = HPADIM
    HPSEC = HPSEC
    TSUCSEC = TSUCSEC
    TDESSEC = TDESSEC
    PSUCSEC = PSUCSEC
    PDESSEC = PDESSEC
    EFICSEC = EFICSEC
    obj = np.array([tipo, ID, PSUCSEC, PDESSEC, TSUCSEC, TDESSEC, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFICSEC, HPSEC, POLLY, Q, RPM])
    salida = np.c_[salida,obj]
    return salida

def sumariocompresor(salida, IMPUL, SECC, COMP, TSUCCOMP, TDESCOMP, PSUCCOMP, PDESCOMP, HPCOMP):
    tipo = "C"
    ID = str(IMPUL) + str(SECC) + str(COMP)
    PSUC = ""
    PDES = ""
    TSUC = ""
    TDES = ""
    DG = ""
    HG = ""
    SURGE = ""
    QN = ""
    STONEW = ""
    CFHEAD = ""
    HEAD = ""
    EFIC = ""
    HP = ""
    POLLY = ""
    Q = ""
    RPM = ""
    TSUCCOMP = TSUCCOMP
    TDESCOMP = TDESCOMP
    PSUCCOMP = PSUCCOMP
    PDESCOMP = PDESCOMP
    obj = np.array([tipo, ID, PSUCCOMP, PDESCOMP, TSUCCOMP, TDESCOMP, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, Q, RPM])
    salida = np.c_[salida,obj]
    return salida

def buildVariablesTeorica(entradas: Matrix):
    datacol = entradas.shape[1]
    comp = entradas[0,1:datacol].astype(int)
    secc = entradas[1,1:datacol].astype(int)
    impul = entradas[2,1:datacol].astype(int)
    impuleq = entradas[3,1:datacol].astype(bool)
    gas = Matrix(entradas[4:19,1:datacol].astype(float))
    diam = entradas[19,1:datacol].astype(float)
    rpm = entradas[20,1:datacol].astype(float)
    tsuc = entradas[21,1:datacol].astype(float)
    psuc = entradas[22,1:datacol].astype(float)
    q = entradas[23,1:datacol].astype(float)
    dimT = entradas[24,1:datacol]
    dimQ = entradas[25,1:datacol]
    dimP = entradas[26,1:datacol]
    dimD = entradas[27,1:datacol]
    caipres = entradas[28,1:datacol].astype(float)
    relvel = entradas[29,1:datacol].astype(float)
    stonew = entradas[30,1:datacol].astype(float)
    surge = entradas[31,1:datacol].astype(float)
    expocp = entradas[32,1:datacol].astype(int)
    cc1 = entradas[33,1:datacol].astype(float)
    cc2 = entradas[34,1:datacol].astype(float)
    cc3 = entradas[35,1:datacol].astype(float)
    cc4 = entradas[36,1:datacol].astype(float)
    expoep = entradas[37,1:datacol].astype(int)
    ce1 = entradas[38,1:datacol].astype(float)
    ce2 = entradas[39,1:datacol].astype(float)
    ce3 = entradas[40,1:datacol].astype(float)
    ce4 = entradas[41,1:datacol].astype(float)
    qExt = entradas[42,1:datacol].astype(float)
    NIMPULS = impul.shape[0]
    suma = 1
    return comp, secc, impul, impuleq, gas, diam, rpm, tsuc, psuc, q, dimT, dimQ, dimP, dimD, caipres, relvel, stonew, surge, expocp, cc1, cc2, cc3, cc4, expoep, ce1, ce2, ce3, ce4, qExt, NIMPULS

def buildVariablesCampo(entradas: Matrix):
    datacol = entradas.shape[1]
    comp = entradas[0,1:datacol].astype(int)
    secc = entradas[1,1:datacol].astype(int)
    gas = Matrix(entradas[2:17,1:datacol].astype(float))
    diam = entradas[17,1:datacol].astype(float)
    rpm = entradas[18,1:datacol].astype(float)
    tsuc = entradas[19,1:datacol].astype(float),
    tdes = entradas[20,1:datacol].astype(float)
    psuc = entradas[21,1:datacol].astype(float)
    pdes = entradas[22,1:datacol].astype(float)
    q = entradas[23,1:datacol].astype(float)
    dimT = entradas[24,1:datacol]
    dimQ = entradas[25,1:datacol]
    dimP = entradas[26,1:datacol]
    dimD = entradas[27,1:datacol]
    NSECC = secc.shape[0]
    return comp, secc, gas, diam, rpm, tsuc, tdes, psuc, pdes, q, dimT, dimQ, dimP, dimD, NSECC

if  __name__ == "__main__":
    simulacion()