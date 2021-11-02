from flask import Flask, json, request, jsonify
from flask.json import JSONEncoder
from time import time
from sympy.matrices import Matrix
import numpy as np
from rutinas.ajustecurvas import *
from rutinas.curvas import *
import decimal

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