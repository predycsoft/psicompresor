import numpy as np
from prettytable import PrettyTable

def AJUSTE(OPCION, SELECCION, X, Y ,NUMERO, EXP):
    print("SELECCION", SELECCION)
    print("X",X)
    print("Y",Y)
    print("NUMERO",NUMERO)
    print("EXP",EXP)

    OPCION = 2
    ORDEN = 2
    M = 0
    P = 100
    if SELECCION == 1:
        P = EXP
    cp1 = 0
    cp2 = 0
    cp3 = 0
    cp4 = 0
    ee = 0
    ERROR1 = 0
    AMN2 = np.zeros(20)
    CMN2 = np.zeros((20,20))
    BMN2 = np.zeros(20)
    DMN2 = np.zeros((20,20))
    YCALCCR = np.zeros(200)
    ERR1 = np.zeros(200)
    AMN2 = np.zeros(20)
    CT = 0
    nm = 0
    NUMERO = NUMERO
    YTOT = np.zeros(NUMERO)
    XTOT = np.zeros(NUMERO)
    XTOT = X
    SELCHP = np.zeros(NUMERO)
    SELEP = np.zeros(NUMERO)
    YMN2 = np.zeros(NUMERO)
    XMN2 = np.zeros(NUMERO)
    YX = np.zeros(50)
    YTOT = Y
    NPROV = 0
    for i in range(NUMERO):
        XMN2[i] = XTOT[i]
        YMN2[i] = YTOT[i]
    NPROV = NUMERO
    P1 = NPROV
    PASO = 1
    XC = np.zeros(M*2)
    if ORDEN == 1:
        MP1 = M + 1
        M2 = M * 2
        for K in range(M2):
            XC[K] = 0
            for i in range(NUMERO):
                XC[K] = XC[K] + XMN2[i] ** K
        YC = 0
        for i in range(NUMERO):
            YC = YC + YMN2[i]
        for K in range(M):
            YX[K] = 0
            for i in range(NUMERO):
                YX[K] = YX[K] + YMN2[i] * XMN2[i] ** K
        for i in range(MP1):
            for j in range(MP1):
                IPJM2 = i + j
                if IPJM2 == 0:
                    CMN2[0, 0] = NUMERO
                else:
                    CMN2[i, j] = XC[IPJM2]
        BMN2[0] = YC
        for i in range(1, MP1):
            BMN2[i] = YX[i - 1]
        INVMATRIZ()
        AMN2 = np.zeros(20)
        for i in range(MP1):
            AMN2[i] = 0
            for j in range(MP1):
                AMN2[i] = AMN2[i] + DMN2[i, j] * BMN2[j]
    #'En este caso el ajuste se hace con un polinomio de grado
    #'mayor a tres.
    if ORDEN == 2:
        NCR = NUMERO
        XCR = np.zeros(NUMERO)
        YCR = np.zeros(NUMERO)
        SUM_1 = np.zeros(NCR)
        SUM_2 = np.zeros(NCR)
        SUM_3 = np.zeros(NCR)
        SUM_4 = np.zeros(NCR)
        SUM_5 = np.zeros(NCR)
        SUM_6 = np.zeros(NCR)
        SUM_7 =np.zeros(NCR)
        SUM_8 = np.zeros(NCR)
        SUM_9 = np.zeros(NCR)
        SUM_10 = np.zeros(NCR)
        SUM_11 =np.zeros(NCR)
        SUM_12 = np.zeros(NCR)
        SUM_13 = np.zeros(NCR)
        for i in range(NCR):
            XCR[i] = XTOT[i]
            YCR[i] = YTOT[i]
        for F in range(NCR):
            if F > 0:
                SUM_1[F] = (XCR[F] ** 2) * YCR[F] + SUM_1[F-1]
                SUM_2[F] = (XCR[F] ** 4) + SUM_2[F-1]
                SUM_3[F] = (XCR[F] ** 3) + SUM_3[F-1]
                SUM_4[F] = (XCR[F] ** 2) + SUM_4[F-1]
                SUM_5[F] = XCR[F] + SUM_5[F-1]
                SUM_6[F] = XCR[F] * YCR[F] + SUM_6[F-1]
                SUM_7[F] = XCR[F] + SUM_7[F-1]
                SUM_8[F] = YCR[F] + SUM_8[F-1]
            else:
                SUM_1[F] = (XCR[F] ** 2) * YCR[F]
                SUM_2[F] = (XCR[F] ** 4)
                SUM_3[F] = (XCR[F] ** 3)
                SUM_4[F] = (XCR[F] ** 2)
                SUM_5[F] = XCR[F]
                SUM_6[F] = XCR[F] * YCR[F]
                SUM_7[F] = XCR[F]
                SUM_8[F] = YCR[F]
        ERROR1 = 1
        for H in range(0,P):
            if SELECCION == 1:
                H = P
            S2 = 0
            B2 = 0
            C2 = 0
            B3 = 0
            D2 = 0
            C3 = 0
            B4 = 0
            S3 = 0
            D3 = 0
            C4 = 0
            S4 = 0
            D4 = 0
            S2 = SUM_1[NCR-1]
            B2 = SUM_2[NCR-1]
            C2 = SUM_3[NCR-1]
            B3 = C2
            D2 = SUM_4[NCR-1]
            C3 = D2
            B4 = D2
            S3 = SUM_6[NCR-1]
            D3 = SUM_7[NCR-1]
            C4 = D3
            S4 = SUM_8[NCR-1]
            D4 = NCR
            if H > 1:
                for G in range(NCR):
                    if G > 0:
                        SUM_9[G] = XCR[G] ** H * YCR[G] + SUM_9[G - 1]
                        SUM_10[G] = XCR[G] ** (2 * H) + SUM_10[G - 1]
                        SUM_11[G] = XCR[G] ** 2 * \
                            XCR[G] ** H + SUM_11[G - 1]
                        SUM_12[G] = XCR[G] * XCR[G] ** H + SUM_12[G - 1]
                        SUM_13[G] = XCR[G] ** H + SUM_13[G - 1]
                    else:
                        SUM_9[G] = XCR[G] ** H * YCR[G]
                        SUM_10[G] = XCR[G] ** (2 * H)
                        SUM_11[G] = XCR[G] ** 2 * XCR[G] ** H
                        SUM_12[G] = XCR[G] * XCR[G] ** H
                        SUM_13[G] = XCR[G] ** H
            S1 = 0
            A1 = 0
            B1 = 0
            C1 = 0
            D1 = 0
            A2 = 0
            A3 = 0
            A4 = 0
            S1 = SUM_9[NCR-1]
            A1 = SUM_10[NCR-1]
            B1 = SUM_11[NCR-1]
            A2 = B1
            C1 = SUM_12[NCR-1]
            A3 = C1
            D1 = SUM_13[NCR-1]
            A4 = D1
            DET1 = A1 * (B2 * C3 * D4 + C2 * D3 * B4 + B3 * C4 *
                            D2 - D2 * C3 * B4 - D3 * C4 * B2 - C2 * B3 * D4)
            DET2 = -B1 * (A2 * C3 * D4 + C2 * D3 * A4 + A3 * C4 *
                            D2 - D2 * C3 * A4 - C2 * A3 * D4 - D3 * C4 * A2)
            DET3 = C1 * (A2 * B3 * D4 + A3 * B4 * D2 + B2 * D3 *
                            A4 - D2 * B3 * A4 - B2 * D4 * A3 - A2 * B4 * D3)
            DET4 = -D1 * (A2 * B3 * C4 + A3 * B4 * C2 + B2 * C3 *
                            A4 - C2 * B3 * A4 - B2 * A3 * C4 - B4 * C3 * A2)
            DETER = DET1 + DET2 + DET3 + DET4
            if H > 1:
                if A1 == 0:
                    if A2 != 0:
                        S5 = S1
                        A5 = A1
                        B5 = B1
                        C5 = C1
                        D5 = D1
                        S1 = S2
                        A1 = A2
                        B1 = B2
                        C1 = C2
                        D1 = D2
                        S2 = S5
                        A2 = A5
                        B2 = B5
                        C2 = C5
                        D2 = D5
                    else:
                        if A3 != 0:
                            S5 = S1
                            A5 = A1
                            B5 = B1
                            C5 = C1
                            D5 = D1
                            S1 = S3
                            A1 = A3
                            B1 = B3
                            C1 = C3
                            D1 = D3
                            S3 = S5
                            A3 = A5
                            B3 = B5
                            C3 = C5
                            D3 = D5
                        else:
                            S5 = S1
                            A5 = A1
                            B5 = B1
                            C5 = C1
                            D5 = D1
                            S1 = S4
                            A1 = A4
                            B1 = B4
                            C1 = C4
                            D1 = D4
                            S4 = S5
                            A4 = A5
                            B4 = B5
                            C4 = C5
                            D4 = D5
                if A2 != 0:
                    S1 = -A2 * S1 / A1
                    B1 = -A2 * B1 / A1
                    C1 = -A2 * C1 / A1
                    D1 = -A2 * D1 / A1
                    A1 = -A2
                    S2 = S1 + S2
                    A2 = A1 + A2
                    B2 = B1 + B2
                    C2 = C1 + C2
                    D2 = D1 + D2
                if A3 != 0:
                    S1 = -A3 / A1 * S1
                    B1 = -A3 / A1 * B1
                    C1 = -A3 / A1 * C1
                    D1 = -A3 / A1 * D1
                    A1 = -A3
                    S3 = S1 + S3
                    A3 = A1 + A3
                    B3 = B1 + B3
                    C3 = C1 + C3
                    D3 = D1 + D3
                if A4 != 0:
                    S1 = -A4 / A1 * S1
                    B1 = -A4 / A1 * B1
                    C1 = -A4 / A1 * C1
                    D1 = -A4 / A1 * D1
                    A1 = -A4
                    S4 = S1 + S4
                    A4 = A1 + A4
                    B4 = B1 + B4
                    C4 = C1 + C4
                    D4 = D1 + D4
            if H > 0:
                if B2 == 0:
                    if B3 != 0:
                        S5 = S2
                        B5 = B2
                        C5 = C2
                        D5 = D2
                        S2 = S3
                        B2 = B3
                        C2 = C3
                        D2 = D3
                        S3 = S5
                        B3 = B5
                        C3 = C5
                        D3 = D5
                    else:
                        S5 = S2
                        B5 = B2
                        C5 = C2
                        D5 = D2
                        S2 = S4
                        B2 = B4
                        C2 = C4
                        D2 = D4
                        S4 = S5
                        B4 = B5
                        C4 = C5
                        D4 = D5
                if B3 != 0:
                    S2 = -B3 / B2 * S2
                    C2 = -B3 / B2 * C2
                    D2 = -B3 / B2 * D2
                    B2 = -B3
                    S3 = S2 + S3
                    B3 = B2 + B3
                    C3 = C2 + C3
                    D3 = D2 + D3
                if B4 != 0:
                    S2 = -B4 / B2 * S2
                    C2 = -B4 / B2 * C2
                    D2 = -B4 / B2 * D2
                    B2 = -B4
                    S4 = S2 + S4
                    B4 = B2 + B4
                    C4 = C2 + C4
                    D4 = D2 + D4
            if C3 == 0:
                S5 = S3
                C5 = C3
                D5 = D3
                S3 = S4
                C3 = C4
                D3 = D4
                S4 = S5
                C4 = C5
                D4 = D5
            if C4 != 0:
                S3 = -C4 / C3 * S3
                D3 = -C4 / C3 * D3
                C3 = -C4
                S4 = S3 + S4
                C4 = C3 + C4
                D4 = D3 + D4
            S4 = S4 / D4
            S3 = (S3 - S4 * D3) / C3
            if H > 0:
                S2 = (S2 - D2 * S4 - C2 * S3) / B2
            if H > 1:
                S1 = (S1 - D1 * S4 - C1 * S3 - B1 * S2) / A1
            if H == 0:
                S1 = 0
            if H == 0:
                S2 = 0
            if H == 1:
                S1 = 0
            for L in range(NCR):
                YCALCCR[L] = S1 * XCR[L] ** H + S2 * \
                    XCR[L] ** 2 + S3 * XCR[L] + S4
            ERR1[0] = (YCALCCR[0] - YCR[0]) ** 2
            for T in range(1, NCR):
                ERR1[T] = (YCALCCR[T] - YCR[T]) ** 2 + ERR1[T - 1]
            #'OJO AQUI HICE UNOS CAMBIOS EN LOS NOMBRES DE UNAS VARIABLES  VER SUB ORIGINAL!!!!
            #'aa;BB;CC
            if ERR1[NCR-1] < ERROR1:
                ERROR1 = ERR1[NCR-1]
                cp1 = S1
                cp2 = S2
                cp3 = S3
                cp4 = S4
                ee = H
    AMN2[3] = cp1
    AMN2[2] = cp2
    AMN2[1] = cp3
    AMN2[0] = cp4
    EXPONENTE = ee
    DIF2TOT = ERROR1
    AJUSTE = PrettyTable(['Coeficiente','Valor'])
    AJUSTE.add_row(['cp1', cp1])
    AJUSTE.add_row(['cp2', cp2])
    AJUSTE.add_row(['cp3', cp3])
    AJUSTE.add_row(['cp4', cp4])
    AJUSTE.add_row(['Exponente', EXPONENTE])
    AJUSTE.add_row(['DIF2TOT', ERROR1])
    print(AJUSTE)
    obj = np.array([cp1,cp2,cp3,cp4,ee,ERROR1])
    return obj


def INVMATRIZ():
    return

if __name__ == "__main__":
    AJUSTE()