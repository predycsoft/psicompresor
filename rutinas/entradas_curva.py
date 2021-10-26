from sympy import Matrix, Array
import numpy as np

z = Matrix((0.480000,0.066000,0.040100,0.005000,0.018000,0.005000,0.006500,0.007400,0.000000,0.000000,0.000000,0.000000,0.330000,0.021000,0.021000))
TSUC = 92.48
PSUC = 106.12
NRPM = 10570
NRPMS = np.array([7400,8450,9500,10570,11100])
CC1= -1.4004
CC2= 21.503
CC3 = -35.491
CC4 = 15.915
EXPOCP = 3
CE1 = 0.7846
CE2 = -1.1637
CE3 = 4.2487
CE4 = -3.9903
EXPOEP = 3
GC = 32.174
SURGE = 0.39463
STONEW = 0.61821
DEQ = 17.7
contadorcurva = 0
caudaltanteo = Array(np.zeros([5,5,15]))
#curvasurge = calculosurge()
#curvastw = calculostw()
step = 10

DEQ_DIM = '[pulg]'
T_DIM = '[°F]'
P_DIM = '[psig]'
FLUJO_DIM = '[ACFM]'