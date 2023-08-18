# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:48:36 2020

@author: Juliana
"""

#VC electroadsorción

#Librerias
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

#Constantes  Unidades 
alpha=0.5                    #Factor de simetría
n=1                          #Número de e-
F=96470        #C/mol
R=8.314        #J/(mol*K)
T=300          #K
v=0.05         #V/s          #Velocidad de barrido
Cb=0.001                     #Relación adimensional concentraciones C*/C0
E0=0.15        #V            #Potencial estandar
Ei=0.3         #V            #Potencial inicial
Ef=-0.3        #V            #Potencial final
p=n*F/(R*T)
tf=(Ei-Ef)/v   #s            #Tiempo un barrido del experimento

#Variables interactivas
g=0.0          #V            #Parámetro de interacción Frumkin
k0=100         #s**-1        #Cte transferencia de carga
C=1e-5         #C/(V*cm**2)  #Capacitancia doble capa eléctrica
Rs=1           #ohm          #Resistencia de la sln
gam=2e-9       #mol/cm**2    #Máximo moles adsorbidas superficie


def sol(g, k0, C, Rs, gam):
    def Adsorp(var,t):
        iT = var[0]
        theta = var[1]
        if t < tf:
            E = Ei - v*t       #Barrido Catódico
            s=1
        else:
            E = Ef + v*(t-tf)  #Barrido Anódico
            s=-1 
        eta=E-E0
        kred=k0*np.exp(-alpha*p*(eta+g*theta+Rs*iT))
        kox=k0*np.exp((1-alpha)*p*(eta+g*theta+Rs*iT))
        dtheta = Cb*(1-theta)*kred-theta*kox
        di = s*v/Rs-iT/(Rs*C)+(n*F*gam*dtheta)/(Rs*C)
        return [di,dtheta]    
    #Vectores
    t = np.linspace(0,2*tf,1000)               #Tiempo
    var0 = [0,0]                               #Variables iniciales
    var = odeint(Adsorp, var0, t)              #Solución
    E1 = np.linspace(0.3, -0.3, 500)           
    E2 = np.linspace(-0.3, 0.3, 500)           
    Ev = np.concatenate((E1,E2), axis=0)       #Potencial
    return Ev, var[:,0]*1e6

y = sol(g, k0, C, Rs, gam)

#Gráfica
fig, ax = plt.subplots()    
plt.subplots_adjust(left=0.15, bottom=0.4)
L1,=ax.plot(y[0], y[1], lw=2, color='royalblue')
ax.set(xlabel='E (V)', ylabel='i ($\mu$A)')


#Sliders
axg = plt.axes([0.20, 0.25, 0.65, 0.03], facecolor='lemonchiffon')
sg = Slider(axg, 'g (V)', -0.1, 0.5, valinit=0, valstep=0.05, color='gold')

axk0 = plt.axes([0.20, 0.20, 0.65, 0.03], facecolor='beige')
sk0 = Slider(axk0, '$log_{10} k_0$', -2, 2, valinit=2, valstep=0.5, color='peru')

axC = plt.axes([0.20, 0.15, 0.65, 0.03], facecolor='mistyrose')
sC = Slider(axC, '$C (\mu F)$', 10, 1010, valinit=50, valstep=50, color='lightcoral')

axRs = plt.axes([0.20, 0.1, 0.65, 0.03], facecolor='lavender')
sRs = Slider(axRs, '$R (\Omega)$', 50, 1000, valinit=1, valstep=50, color='cornflowerblue')

axGa = plt.axes([0.20, 0.05, 0.65, 0.03], facecolor='lightcyan')
sGa = Slider(axGa, '$\Gamma_{max} x10^{10}$', 5, 105, valinit=20, valstep=10, color='turquoise')


def update(val):
    k=sol(sg.val, 10**sk0.val, sC.val*1e-6, sRs.val, sGa.val*1e-10)
    L1.set_ydata(k[1])
    L1.set_xdata(k[0])
    ax.set_ylim(k[1].min(), k[1].max())
    ax.set_xlim(k[0].min(), k[0].max())
    fig.canvas.draw_idle()

sg.on_changed(update)
sk0.on_changed(update)
sC.on_changed(update)
sRs.on_changed(update)
sGa.on_changed(update)

resetax = plt.axes([0.8, 0.9, 0.1, 0.04])
button = Button(resetax, 'Reset', color='lightgoldenrodyellow', hovercolor='0.975')


def reset(event):
    sg.reset()
    sk0.reset()
    sC.reset()
    sRs.reset()
    sGa.reset()
button.on_clicked(reset)


plt.show()
