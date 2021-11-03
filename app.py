from flask import Flask, json, request, jsonify
from flask.json import JSONEncoder
from time import time
from sympy.matrices import Matrix
import numpy as np
from rutinas.ajustecurvas import *
from rutinas.curvas import *
import decimal

from rutinas.entradas_punto import NIMPULS

app = Flask(__name__)



@app.route('/ajustecurva/', methods=['POST'])
def hello():
    headers = {'Access-Control-Allow-Origin': '*'}
    req_data = json.loads(request.data.decode('utf-8'))
    print(req_data)
    datos = np.asanyarray(req_data)
    SELECCION = int(datos[1,1])
    ORDEN = int(datos[1,3])
    X1 = np.array(datos[2:datos.shape[0],0], dtype=np.float64)
    Y1 = np.array(datos[2:datos.shape[0],1] , dtype=np.float64)
    X2 = np.array(datos[2:datos.shape[0],2] , dtype=np.float64)
    Y2 = np.array(datos[2:datos.shape[0],3] , dtype=np.float64)
    X1 = X1[X1 != 0]
    Y1 = Y1[Y1 != 0]
    X2 = X2[X2 != 0]
    Y2 = Y2[Y2 != 0]
    ARR1 = AJUSTE(2,SELECCION,X1,Y1,datos.shape[0]-2,ORDEN)
    ARR2 = AJUSTE(2,SELECCION,X2,Y2,datos.shape[0]-2,ORDEN)
    polinomios = np.c_[ARR1, ARR2]
    response = jsonify(polinomios.tolist())
    response.headers = headers
    response.status_code = 200
    return response


@app.route('/adimensional/', methods=['POST'])
def adimensional():
    headers = {'Access-Control-Allow-Origin': '*'}
    req_data = json.loads(request.data.decode('utf-8'))
    print(req_data)
    datos = np.asanyarray(req_data)
    varAdim = np.array(["Eficiencia Politropica Np", "Coeficiente de Work Input", "Coeficiente de Cabezal Politropico", "Trabajo Politropico", "Potencia de gas", "Flujo Másico" ,"Relacion de Compresion", "Relacion de Volumen", "Temperatura Isentropica", "Presion Isentropica", "Densidad de Succion","Densidad de descarga","Densidad Isentropica", " Volumen en succion", "Volumen en descarga", "Volumen Isentropico", "Entalpia de succion", "Entalpia de descarga", "Entalpia isentropica", "Entropia de succion", "Entropia de descarga", "Entropia isentropica","Factor de compresion succion","factor de compresion descarga", "factor de compresion isentropica", "peso molecular", "Coeficiente de flujo QN"])
    z = Matrix((datos[1,0:15])).T
    DEQ = (datos[1,15]).astype(np.float64)
    TSUC = (datos[1,16]).astype(np.float64)
    PSUC = (datos[1,17]).astype(np.float64)
    TDES = (datos[1,18]).astype(np.float64)
    PDES = (datos[1,19]).astype(np.float64)
    FLUJO = (datos[1,20]).astype(np.float64)
    RPM = (datos[1,21]).astype(np.float64)
    NIMPULS = 1
    DDIM = datos[1,22]
    TDIM = datos[1,23]
    PDIM = datos[1,24]
    QDIM = datos[1,25]
    arrAdim = ADIM(z, DEQ, TSUC, PSUC, TDES, PDES, FLUJO, RPM, NIMPULS, DDIM, TDIM, TDIM, PDIM, PDIM, QDIM)
    varAdim = np.c_[varAdim, np.array(arrAdim)]
    response = jsonify(varAdim.astype(str).tolist())
    response.headers = headers
    response.status_code = 200
    return response

@app.route('/pruebaEficiencia/', methods=['POST'])
def pruebaEficiencia():
    headers = {'Access-Control-Allow-Origin': '*'}
    req_data = json.loads(request.data.decode('utf-8'))
    print(req_data)
    datos = np.asanyarray(req_data)
    salida =np.array(["tipo", "ID", "PSUC", "PDES", "TSUC", "TDES", "Densidad", "Entalpia", "SURGE", "QN", "STONEW", "CFHEAD", "HEAD", "EFIC", "HP", "POLLY", "Q", "RPM"])
    salida = np.c_[salida]
    z = Matrix((datos[1,0:15])).T
    TSUC = (datos[1,15]).astype(np.float64)
    PSUC = (datos[1,16]).astype(np.float64)
    FLUJO = (datos[1,17]).astype(np.float64)
    DIAM = (datos[1,18]).astype(np.float64)
    RPM = (datos[1,19]).astype(np.float64)
    CC4 = (datos[1,20]).astype(np.float64)
    CC3 = (datos[1,21]).astype(np.float64)
    CC2 = (datos[1,22]).astype(np.float64)
    CC1 = (datos[1,23]).astype(np.float64)
    EXPOCP = (datos[1,24]).astype(np.float64)
    CE4 = (datos[1,25]).astype(np.float64)
    CE3 = (datos[1,26]).astype(np.float64)
    CE2 = (datos[1,27]).astype(np.float64)
    CE1 = (datos[1,28]).astype(np.float64)
    EXPOEP = (datos[1,29]).astype(np.float64)
    SURGE = (datos[1,30]).astype(np.float64)
    STONEW = (datos[1,31]).astype(np.float64)
    GC = 32.174
    DDIM = datos[1,32]
    TDIM = datos[1,33]
    PDIM = datos[1,34]
    QDIM = datos[1,35]
    PDESCAMPO = (datos[1,36]).astype(np.float64)
    PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, FLUJO, DIAM, RPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DDIM, TDIM, PDIM, QDIM)
    XTAN1 = RPM
    YTAN1 = PDES
    XTAN2 = RPM*0.95
    PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, FLUJO, DIAM, XTAN2, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DDIM, TDIM, PDIM, QDIM)
    YTAN2 = PDES
    while 1:
        RPM = (XTAN1 * (YTAN2 - PDESCAMPO) - XTAN2 * (YTAN1 - PDESCAMPO)) / (YTAN2 - YTAN1)
        PDES, TDES, DG, HG, QN, CFHEAD, HEAD, EFIC, HP, POLLY = punto(z, TSUC, PSUC, FLUJO, DIAM, RPM, CC1, CC2, CC3, CC4, EXPOCP, CE1, CE2, CE3, CE4, EXPOEP, GC, SURGE, STONEW, DDIM, TDIM, PDIM, QDIM)
        # Asignacion de nuevos valores del tante
        if abs(XTAN1 - RPM) >= abs(XTAN2 - RPM):
            XTAN1 = RPM
            YTAN1 = PDES
        else:
            XTAN2 = RPM
            YTAN2 = PDES
        ## Chequeo convergencia
        if abs(PDESCAMPO-PDES) <= 0.001:
            break
    obj = np.array(["s", "id", PSUC, PDES, TSUC, TDES, DG, HG, SURGE, QN, STONEW, CFHEAD, HEAD, EFIC, HP, POLLY, FLUJO, RPM])
    salida = np.c_[salida,obj]
    response = jsonify(salida.astype(str).tolist())
    response.headers = headers
    response.status_code = 200
    return response

@app.route('/generarMapa/', methods=['POST'])
def generarMapa():
    headers = {'Access-Control-Allow-Origin': '*'}
    req_data = json.loads(request.data.decode('utf-8'))
    datos = np.asanyarray(req_data)
    # for i in range(datos.shape[0]):
    #         row = []
    #         for j in range(datos.shape[1]):
    #             row.append(datos[i][j])
    #         print(row)
    mapa = np.array(["RPM","Q/N","Temperatura Descarga", "Presion Descarga","Potencia", "Caudal [ACFM]", "Eficiencia", "C. head", "C. Work In", "Head"])
    startCurva = time.time()
    for i in range(1,datos.shape[0]):
        start = time.time()
        z = Matrix((datos[i,0:15])).T
        TSUC = (datos[i,15]).astype(np.float64)
        PSUC = (datos[i,16]).astype(np.float64)
        RPM = (datos[i,17]).astype(np.float64)
        CC1 = (datos[i,18]).astype(np.float64)
        CC2 = (datos[i,19]).astype(np.float64)
        CC3 = (datos[i,20]).astype(np.float64)
        CC4 = (datos[i,21]).astype(np.float64)
        EXPOCP = (datos[i,22]).astype(np.float64)
        CE1 = (datos[i,23]).astype(np.float64)
        CE2 = (datos[i,24]).astype(np.float64)
        CE3 = (datos[i,25]).astype(np.float64)
        CE4 = (datos[i,26]).astype(np.float64)
        EXPOEP = (datos[i,27]).astype(np.float64)
        SURGE = (datos[i,28]).astype(np.float64)
        STONEW = (datos[i,29]).astype(np.float64)
        DEQ = (datos[i,30]).astype(np.float64)
        GC = 32.174
        DDIM = datos[i,31]
        TDIM = datos[i,32]
        PDIM = datos[i,33]
        QDIM = datos[i,34]
        step = 10
        print("TSUC", TSUC, "PSUC", PSUC, "RPM", RPM, "CC1", CC1, "CC2", CC2, "CC3", CC3, "CC4", CC4, "EXPOCP", EXPOCP)
        print("CE1", CE1, "CE2", CE2, "CE3", CE3, "CE4", CE4, "EXPOEP", EXPOEP, "SURGE", SURGE, "STONEW", STONEW, "DEQ", DEQ)
        print("DDIM", DDIM, "TDIM", TDIM, "PDIM", PDIM, "QDIM", QDIM)
        print(z)
        ARR = curvas(True, SURGE, STONEW, z, TSUC, PSUC, DEQ, RPM, CC4, CC3, CC2, CC1, EXPOCP, CE4, CE3, CE2, CE1,EXPOEP, GC, step, DDIM, TDIM, PDIM, "[ACFM]")
        mapa = np.c_[mapa, ARR]
        end = time.time()
        print("tiempo de grafica de curva", end - start)
    endCurva = time.time()
    print("tiempo de grafica del mapa", endCurva - startCurva)

    s = jsonify(mapa.astype(str).tolist())
    print(s)
    response = s
    response.headers = headers
    response.status_code = 200
    return response

app.run(debug=True)