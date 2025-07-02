import numpy as np
from cmath import log10, rect
import math

def Celsius_to_Kelvin(C):
    return (C + 273.15)

def P2R(amplitude, angles):
    nprect = np.vectorize(rect)
    return nprect(amplitude, np.deg2rad(angles))

def R2P(x):
    return abs(x),np.angle(x,deg=True)

def DB10(x):
  # to power DB
  x=abs(x)
  return 10*log10(x)

def FDB10(x):  
  #from power DB
  return round(10**(x/10),3)

def DBm10(x):
    #from power DB
    return DB10(x)+30

def DB20(x):
  # to voltage DB
  x=abs(x)
  return 20*log10(x)

def FDB20(x):  
  #from voltage DB
  return round(10**(x/20),3)

def calc_circle(c, r):
  theta = np.linspace(0, 2*np.pi, 1000)
  return c + r*np.exp(1.0j*theta)




def gamas(zs, z0):
    # zs is unnormalized
    a=(zs-z0)/(zs+z0)
    print("gamas read from smith chart :",R2P(a))
    return a

def gamal(gs, s11, s12, s21, s22):
    return np.conj(s22+((s12*s21*gs)/(1-(s11*gs))))

def findzin(zs):
    zin=-zs
    if zin.real<0:
        return zin
    print("zin real part most be negative")
    return False

def findzout(zl):
    zout=-zl
    if zout.real<0:
        return zout
    print("zout real part most be negative")
    return False

def findzfromgama(g, z0):
    # find z from gama
    return [z0*((g+1)/(1-g)),z0*((-g+1)/(1+g))]




def swr(g):
    return (1+abs(g))/(1-abs(g))

def Bandwidth(s21_list,freq_step):
    sum=0
    for i in s21_list:
        sum=sum+(abs(i)**2)*freq_step
    return sum/(abs(np.amax(s21_list))**2)    