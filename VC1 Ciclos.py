# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 23:13:20 2021

@author: Juliana
"""

#VC Diffusion
#Effect of number of cycles 

import numpy as np
import matplotlib.pyplot as plt


num=int(input('Enter the number of cycles:'))


#Constants
R=8.314 
T=298
F=96485
E0, Ei, Ef = 0.0, 0.3, -0.3    #Standard, initial and final potential 
Dm=0.45                        #Dimensionless constant of diffusion Dm=D*dt/dx**2
dE=(Ei-Ef)/50                  #Potential step 
n=1                            #Number of e-
k0=0.5                         #Charge transfer constant
alfa=0.5                       #Symmetry factor
points=100*num


#Vectors
fn=np.ones(100)                #f_new time j+1 (concentration)
fo=np.ones(100)                #f_old time j ; f[i] position
cor=np.zeros(points)              #Current
pot=np.zeros(points)              #Potential
pot[0]=Ei


nc=0                           #Cycle counter

for j in range(1, points):                              #j time counter (dt)
    for i in range(1, 98):                      #i space counter (dx)
        fn[i]=fo[i]+Dm*(fo[i+1]-2*fo[i]+fo[i-1])    
    if j%100==0:
        nc+=1
    if (j-(nc*100))<50:                              #Cathodic sweep
        pot[j]=Ei-(j-(nc*100))*dE  
    else:                                            #Anodic sweep 
        pot[j]=Ef+(j-(nc*100)-50)*dE
    sp=pot[j]-E0                                     #Overpotential
    kf=k0*np.exp(-alfa*n*F*sp/(R*T))                 #Forward constant
    kb=k0*np.exp((1-alfa)*n*F*sp/(R*T))              #Backward constant
    fn[0]=(kb+fn[1])/(kb+kf+1)                       #Acople difusión-transferencia carga en la superficie
    cor[j]=fn[1]-fn[0]                               #Corriente adimensional
    fn=fo                                            #Actualización de la función


# Graph  
plt.plot(pot, cor, color='limegreen')
plt.xlabel('$E\ (V)$')
plt.ylabel('$i\ s(A)$')
plt.title("Cyclic Voltammetry")
plt.text(-0.3, -0.15, 'Number of cycles: %d' %num, bbox={'facecolor': 'mintcream', 'alpha': 0.4, 'pad': 10})
plt.show()
