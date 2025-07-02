import numpy as np
import math
from .basic_calculations import DB10

def rollet(s11,s12,s21,s22):
  delta=(s11*s22)-(s12*s21)
  print('delta:',delta)
  amps11=abs(s11)**2
  amps22=abs(s22)**2
  ampdelta=abs(delta)**2
  k=(1-amps11-amps22+ampdelta)/(2*abs(s12*s21))
  print('k:',k,'|delta|:',abs(delta))
  if amps11<1 and amps22<1 and abs(delta)<1 and k>1:
    print('unconditionally stable')
  else:
    print('potentially unstable check stable circules')  

def mou(s11,s12,s21,s22,verbus=True):
    delta=(s11*s22)-(s12*s21)
    amps11=abs(s11)**2
    mou=(1-amps11)/(abs(s22-(delta*s11.conjugate()))+abs(s12*s21))
    if verbus:
        print('delta:',delta)
        print('mou:',mou)
        if mou>1:
            print('unconditionally stable') 
    return mou

def MAG(s11,s12,s21,s22):
    # Maximum Available Gain for unconditionaly stable
    delta=(s11*s22)-(s12*s21)
    amps11=abs(s11)**2
    amps22=abs(s22)**2
    ampdelta=abs(delta)**2
    k=(1-amps11-amps22+ampdelta)/(2*abs(s12*s21))
    if k<1:
        raise("it is not unconditionaly stable")
    B1=1+amps11-amps22-ampdelta
    if B1<0:
        MAG=DB10(abs(s21/s12))+DB10(abs(k+math.sqrt(k**2-1)))
    else:                                
        MAG=DB10(abs(s21/s12))+DB10(abs(k-math.sqrt(k**2-1)))
    return MAG                    