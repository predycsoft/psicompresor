from rutinas.proter import *
from prettytable import PrettyTable
import math

"""
Entradas:

  DEQ: Diametro Equivalente
  TSUC: Temperatura de succion
  PSUC: Presion de succion
  TDES: Temperatura de descarga
  PDES: Presion de descarga
  qw: Flujo
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

Calculos:
  dc: Densidad Calculada en [lbmol/pie3]
  hm: Entalía Calculada
  sm: Entropía Calculada
  dg: Densidad Calculada en [lbmol/pie3]
  ymw: Peso Molecular de la mezcla

  Auxiliares de funciones:
    XTAN1 XTAN2 YTAN1 YTAN2: Variables para el tanteo isentrópico

Funciones:
  Proter (Ver rutina Proter): Calcula propiedades termodinámicas a través de la ecuación BWRS

Esta rutina retorna los valores de:
    'Flujo Masico QW': qw
    'Potencia Gas Hp QW:': gashp
    'Eficiencia Politropica Np:': Np
    'Coeficiente de Flujo QN:': qn
    'Coeficiente de Cabezal Politropico:': mi
    'Coeficiente Work Input:': si
    'Trabajo Politrópico:': wp
    'Relación de compresión:': relacion_de_compresion
    'Relación de Volumen:': relacion_de_volumen
    'Temperatura de Succión [°F]:': t1s
    'Temperatura de Descarga [°F]:': t2d
    'Temperatura de Isentrópica [°F]:': t3i
    'Presión de Succión [psia]': p1s
    'Presión de Descarga [psia]': p2d
    'Presión de Isentrópica [psia]': p3i
    'Densidad de Succión [lbm/pie3]': dens_suc
    'Densidad de Descarga [lbm/pie3]': dens_des
    'Densidad de Isentrópica [lbm/pie3]': dens_isen
    'Volumen de Succión [pie3/lbm]': vm1s
    'Volumen de Descarga  [pie3/lbm]': vm2d
    'Volumen de Isentrópica  [pie3/lbm]': vg3i
    'Entalpía de Succión [BTU/lbm]': hm1s
    'Entalpía de Descarga [BTU/lbm]': hm2d
    'Entalpía de Isentrópica [BTU/lbm]': hm3i
    'Entropía de Succión BTU/lbmol°R': sg1s
    'Entropía de Descarga BTU/lbmol°R': sg2d
    'Entropía de Isentrópica BTU/lbmol°R': sg3i
    'Calidad de Succión': v1s
    'Calidad de Descarga': v2d
    'Calidad de Isentrópica': v3i
    'Factor de Compresión de Succión': zcomp
    'Factor de Compresión de Descarga':  zsal
    'Factor de Compresión de Isentrópica': fact_comp_isen
    'Peso Molecular': ymw
"""

def ADIM(z, DEQ, TSUC, PSUC, TDES, PDES, qw, NREV, NIMPULS, DEQ_DIM, TSUC_DIM, TDES_DIM, PSUC_DIM, PDES_DIM, FLUJO_DIM):
    V = 1
    ############################################
    ######## AJUSTE DIMENSIONAL ################
    #DIAMETRO
    if (DEQ_DIM == '[pulg]'):
        DEQ = DEQ
    elif (DEQ_DIM == '[mm]'):
        DEQ = DEQ / 25.4
    elif (DEQ_DIM == '[pies]'):
        DEQ = DEQ * 12
    #Temperatura en la succión.
    if (TSUC_DIM== '[ºF]'):
        TSUC = TSUC
    elif (TSUC_DIM== '[ºK]'):
        TSUC = TSUC * 1.8 - 459.67
    elif (TSUC_DIM== '[ºC]'):
        TSUC = TSUC * 1.8 + 32
    #Temperatura en la descarga.
    if (TDES_DIM == '[ºF]'):
        TDES = TDES
    elif (TDES_DIM == '[ºK]'):
        TDES = TDES * 1.8 - 459.67
    elif (TDES_DIM == '[ºC]'):
        TDES = TDES * 1.8 + 32
    #Presión en la succión.
    if (PSUC_DIM == '[psia]'):
        PSUC = PSUC - 14.7
    elif (PSUC_DIM == '[psig]'):
        PSUC = PSUC
    elif (PSUC_DIM == '[KPag]'):
        PSUC = PSUC * 0.14503773773
    elif (PSUC_DIM == '[barg]'):
        PSUC = PSUC * 14.503773773
    #Presión en la descarga
    if (PDES_DIM == '[psia]'):
        PDES = PDES - 14.7
    elif (PDES_DIM == '[psig]'):
        PDES = PDES
    elif (PDES_DIM == '[KPag]'):
        PDES = PDES * 0.14503773773
    elif (PDES_DIM == '[barg]'):
        PDES = PDES * 14.503773773
    ############################################
    ###### Etapa de succion ####################
    P = PSUC + 14.7
    T = TSUC + 459.69
    dc, hg, sg, dg, ymw = PROTER(z, T, P)
    p1s = P
    t1s = T - 459.69
    hm1s = hg / ymw
    vm1s = 1 / dg
    sg1s = sg
    v1s = V
    ############################################
    ###### Etapa de descarga ###################
    P = PDES + 14.7
    T = TDES + 459.69
    DC, hg, sg, dg, ymw = PROTER(z, T, P)
    p2d = P
    t2d = T - 459.69
    hm2d = hg/ymw
    vm2d = 1/dg
    sg2d = sg
    v2d = V
    ############################################
    ##### Calculo isentropico ##################
    XTANT1 = T
    YTANT1 = sg
    T = T*0.95
    DC, hg, sg, dg, ymw = PROTER(z, T, P)
    XTANT2 = T
    YTANT2 = sg
    if abs(sg - sg1s) >= 0.00000001:
        while 1:
            T = (XTANT1 * (YTANT2 - sg1s) - XTANT2 * (YTANT1 - sg1s)) / (YTANT2 - YTANT1)
            DC, hg, sg, dg, ymw = PROTER(z, T, P)
            # Asignacion de nuevos valores del tante
            if abs(XTANT1 - T) >= abs(XTANT2 - T):
                XTANT1 = T
                YTANT1 = sg
            else:
                XTANT2 = T
                YTANT2 = sg
            ## Chequeo convergencia
            if abs(sg-sg1s) <= 0.00000001:
                break
    sg3i = sg
    t3i = T - 459.69
    p3i = P
    hm3i = hg / ymw
    vg3i = 1 / dg
    v3i = V
    ############################################
    ###### Preparo la salida de ADIM ###########
    if FLUJO_DIM == "[lbm/min]":
        qw = qw
    elif  FLUJO_DIM == "[MMSCFD]":
        qw = qw * 2634.5974 * ymw / 1440
    elif  FLUJO_DIM == "[lbmol/dia]":
        qw = qw * ymw / 1440
    elif  FLUJO_DIM == "[ACFM]":
        qw = qw / vm1s
    elif  FLUJO_DIM == "[m3/min]":
        qw = qw / vm1s * 35.31466672
    elif  FLUJO_DIM == "[m3/seg]":
        qw = qw  / vm1s * 35.31466672 * 60

    Nyo = (math.log(p2d/p1s))/(math.log(vm1s/vm2d))
    Yo = Nyo / (Nyo - 1)
    M1 = (math.log(p2d/p1s))/(math.log(vm1s/vg3i))
    X1 = M1 / (M1 -1)
    FG = (hm3i - hm1s) * 778.16 / (X1 * (p2d * vg3i - p1s * vm1s) * 144)
    wp = FG * Yo * (p2d * vm2d - p1s * vm1s) * 144
    Np = wp / ((hm2d - hm1s) * 778.16)
    mi = 32.174 * wp / NIMPULS / (3.1415927 * DEQ * NREV / 720) ** 2
    qn = qw * vm1s/NREV
    si = mi/Np
    gashp = qw * 60 * (hm2d - hm1s) * 0.2931 * 0.001341
    zcomp = p1s * vm1s * ymw / ((t1s + 459.69) * 10.7335)
    zsal = p2d * vm2d * ymw / ((t2d + 459.69) * 10.7335)
    relacion_de_compresion = p2d/p1s
    relacion_de_volumen = vm1s/vm2d
    dens_suc = 1/vm1s
    dens_des = 1/vm2d
    dens_isen = 1/vg3i
    fact_comp_isen = p3i * vg3i * ymw / ((t3i + 459.69) * 10.7335)

    table = PrettyTable(['Propiedad', 'Valor'])
    table.add_row(['Flujo Masico QW', qw])
    table.add_row(['Potencia Gas Hp QW:', gashp])
    table.add_row(['Eficiencia Politropica Np:', Np])
    table.add_row(['Coeficiente de Flujo QN:', qn])
    table.add_row(['Coeficiente de Cabezal Politropico:', mi])
    table.add_row(['Coeficiente Work Input:', si])
    table.add_row(['Trabajo Politrópico:', wp])
    table.add_row(['Relación de compresión:', relacion_de_compresion])
    table.add_row(['Relación de Volumen:', relacion_de_volumen])
    table.add_row(['Temperatura de Succión [°F]:', t1s])
    table.add_row(['Temperatura de Descarga [°F]:', t2d])
    table.add_row(['Temperatura de Isentrópica [°F]:', t3i])
    table.add_row(['Presión de Succión [psia]', p1s])
    table.add_row(['Presión de Descarga [psia]', p2d])
    table.add_row(['Presión de Isentrópica [psia]', p3i])
    table.add_row(['Densidad de Succión [lbm/pie3]', dens_suc])
    table.add_row(['Densidad de Descarga [lbm/pie3]', dens_des])
    table.add_row(['Densidad de Isentrópica [lbm/pie3]', dens_isen])
    table.add_row(['Volumen de Succión [pie3/lbm]', vm1s])
    table.add_row(['Volumen de Descarga  [pie3/lbm]', vm2d])
    table.add_row(['Volumen de Isentrópica  [pie3/lbm]', vg3i])
    table.add_row(['Entalpía de Succión [BTU/lbm]', hm1s])
    table.add_row(['Entalpía de Descarga [BTU/lbm]', hm2d])
    table.add_row(['Entalpía de Isentrópica [BTU/lbm]', hm3i])
    table.add_row(['Entropía de Succión BTU/lbmol°R', sg1s])
    table.add_row(['Entropía de Descarga BTU/lbmol°R', sg2d])
    table.add_row(['Entropía de Isentrópica BTU/lbmol°R', sg3i])
    table.add_row(['Calidad de Succión', v1s])
    table.add_row(['Calidad de Descarga', v2d])
    table.add_row(['Calidad de Isentrópica', v3i])
    table.add_row(['Factor de Compresión de Succión', zcomp])
    table.add_row(['Factor de Compresión de Descarga',  zsal])
    table.add_row(['Factor de Compresión de Isentrópica', fact_comp_isen])
    table.add_row(['Peso Molecular', ymw])
    #print(table)
    return Np, si, mi, wp, gashp, qw, relacion_de_compresion, relacion_de_volumen, t3i, p3i,dens_suc, dens_des, dens_isen, vm1s, vm2d, vg3i, hm1s, hm2d, hm3i, sg1s, sg2d, sg3i, zcomp, zsal, fact_comp_isen, ymw, qn
