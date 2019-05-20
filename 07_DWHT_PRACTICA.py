# -*- coding: utf-8 -*-


########################################################

import numpy as np
import matplotlib.pyplot as plt


"""
Implementar una funcion H_WH(N) que devuelva la matriz NxN asociada a la transformación de Walsh-Hadamard

H_WH(4)=
      [[ 0.5,  0.5,  0.5,  0.5],
       [ 0.5,  0.5, -0.5, -0.5],
       [ 0.5, -0.5, -0.5,  0.5],
       [ 0.5, -0.5,  0.5, -0.5]]
"""
def calculocambios(vec):
    count=0
    for i in range(1,len(vec)):
        if vec[i-1]!=vec[i]:
            count=count+1
    return count

def multmat(scalar,matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            matrix[i][j]=matrix[i][j]*scalar
    return matrix

def fusiondevec(m1,m2):
    mres=[]
    for i in range(len(m1)):
        aux=m1[i]+m2[i]
        mres=mres+[aux]
    return mres

def compose(matrix):
    mat1=matrix
    mat2=matrix
    mat3=matrix
    mat4=multmat(-1,matrix)
    m12=fusiondevec(mat1,mat2)
    print("m12:",str(m12))
    m34=fusiondevec(mat3,mat4)
    print("m34:",str(m34))
    return m12+m34

print(compose([[1,1],[1,-1]]))

def defmataux(N):
    if N==1 :
        return [[1,1],[1,-1]]
    else:
        mataux=defmataux(N-1)
        return compose(mataux)

def H_WH(N):
    maux=defmataux(N)# calculo de mat auxiliar
    scores = []
    for itm in range(N):
        scores.append( (calculocambios(maux[itm]),itm) )
    scores=sorted(scores)#calculo de la diferencias de cada vector para ordenar las filas
    result=[]
    for s in scores:
        result.append(mataux[s[1]])#reordenacion de la matrix auxiliar
    return multmat(1/sqrt(N),result)

        



"""
Implementar la DWHT (Discrete Walsh-Hadamard Transform) y su inversa
para bloques NxN

dwht_bloque(p) 
idwht_bloque(p) 

p bloque NxN

dwht_bloque(
            [[217,   8, 248, 199],
             [215, 189, 242,  10],
             [200,  65, 191,  92],
             [174, 239, 237, 118]]
            )=
            [[ 661,   -7.5, -48.5, 201],
             [   3,  -27.5,  25.5,  57],
             [  59,  -74.5,  36.5, -45],
             [ -51, -112.5, 146.5,  45]]

"""
#def dwht_bloque(p):




#def idwht_bloque(p):


"""
Reproducir los bloques base de la transformación para los casos N=4,8,16
Ver imágenes adjuntas
"""


