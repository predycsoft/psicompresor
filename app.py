from flask import Flask, json, request, jsonify

from rutinas.ajustecurvas import AJUSTE
import numpy as np

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

