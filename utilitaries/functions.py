import skrf as rf
import numpy as np
import schemdraw
import schemdraw.elements as elm
from skrf import Network,Frequency
from scipy.optimize import minimize
from skrf.media import DistributedCircuit
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from matplotlib import style
from skrf.data import ring_slot
from cmath import log10,rect
from pylab import *
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





def sensitivity(Temp,bandwidth,noisefig,snr):
    # normaly negative
    Temp=Celsius_to_Kelvin(Temp)
    return DBm10(KTB(Temp,bandwidth))+noisefig+snr

def Bandwidth(s21_list,freq_step):
    sum=0
    for i in s21_list:
        sum=sum+(abs(i)**2)*freq_step
    return sum/(abs(np.amax(s21_list))**2)    
    
def KTB(Temp,bandwidth):
    Temp=Celsius_to_Kelvin(Temp)
    return 1.38065 * 10**-23*Temp*bandwidth

def Thermal_noise_floor(Temp,bandwidth):
    return DBm10(KTB(Temp,bandwidth))

def swr(g):
    return (1+abs(g))/(1-abs(g))
    
def GT(s11,s12,s21,s22,gs,gl):
    # Transducer Gain : 
    # the actual gain of an amplifier stage including the effects of input and output matching and device gain
    qabs = lambda x: np.square(np.absolute(x))
    return (qabs(s21)*(1-qabs(gs))*(1-qabs(gl)))/qabs((1-s11*gs)*(1-s22*gl)-s12*s21*gl*gs)

def Rout(zout,zl):
    "series resistance"
    return -zout.real-zl.real

def Rin(zin,zs):
    "series resistance"
    return -zin.real-zs.real

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

def gamas(zs, z0):
    # zs is unnormalized
    a=(zs-z0)/(zs+z0)
    print("gamas read from smith chart :",R2P(a))
    return a

def gamal(gs, s11, s12, s21, s22):
    return np.conj(s22+((s12*s21*gs)/(1-(s11*gs))))

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

def smithcircule(R):
    fig = go.Figure(go.Scattersmith(imag=np.imag(R), real=np.real(R),marker_color="red",showlegend=True,name='guess'))
    fig.update_layout(height=800,width=1000)
    fig.show()
    
def calc_circle(c, r):
  theta = np.linspace(0, 2*np.pi, 1000)
  return c + r*np.exp(1.0j*theta)

def plotstablecircules(s11,s12,s21,s22,R):
    sqabs = lambda x: np.square(np.absolute(x))
    #plotlist=[]
    delta = (s11*s22) - (s12*s21)
    rl = np.absolute((s12 * s21)/(sqabs(s22) - sqabs(delta)))
    cl = np.conj(s22 - delta*np.conj(s11))/(sqabs(s22) - sqabs(delta))
    print("out (","R:",rl,"C:",R2P(cl),")")
    rs = np.absolute((s12 * s21)/(sqabs(s11) - sqabs(delta)))
    cs = np.conj(s11 - delta*np.conj(s22))/(sqabs(s11) - sqabs(delta))
    print("in (","R:",rs,"C:",R2P(cs),")")
    plt.figure(figsize=(10, 8), dpi=80)
    n = rf.Network(name="out", s=calc_circle(cl, rl))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(cl, rl)), real=np.real(calc_circle(cl, rl)),marker_color="green",showlegend=True,name='out'))
    n = rf.Network(name="in", s=calc_circle(cs, rs))
    n.plot_s_smith(lw=2,draw_labels=True)  
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(cs, rs)), real=np.real(calc_circle(cs, rs)),marker_color="red",showlegend=True,name='in'))
    for i in R:
        n = rf.Network(name=str(i), s=calc_circle(0, i))
        n.plot_s_smith(lw=1,draw_labels=True)      
        #plotlist.append(go.Scattersmith(imag=[0], real=[i],showlegend=True,name=str(i)))
    #fig = go.Figure(plotlist)
    #fig.update_layout(height=800,width=1000)
    #fig.show()      
    if abs(s11)<1 and abs(s22)<1:
        print('hole smith except circles')
    else:
        print('inside of circles')   

def gain_ciecles(s11,gama_s,s22,gama_l,guess, gamma_opt):
    sqabs = lambda x: np.square(np.absolute(x))
    #plotlist=[]
    gs=(( 1-sqabs(gama_s) / sqabs(1-s11*gama_s)))*(1-sqabs(s11))
    gl=(( 1-sqabs(gama_l) / sqabs(1-s22*gama_l)))*(1-sqabs(s22))

    Cs = (gs*np.conjugate(s11))/(1-(1-gs)*sqabs(s11))
    Rs = (np.sqrt(1-gs)*(1-sqabs(s11)))/(1-(1-gs)*sqabs(s11))

    Cl = (gl*np.conjugate(s22))/(1-(1-gl)*sqabs(s22))
    Rl = (np.sqrt(1-gl)*(1-sqabs(s22)))/(1-(1-gl)*sqabs(s22))

    plt.figure(figsize=(10, 8), dpi=80)
    n = rf.Network(name="out", s=calc_circle(Cl, Rl))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(Cl, Rl)), real=np.real(calc_circle(Cl, Rl)),marker_color="green",showlegend=True,name='out'))
    n = rf.Network(name="in", s=calc_circle(Cs, Rs))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(Cs, Rs)), real=np.real(calc_circle(Cs, Rs)),marker_color="red",showlegend=True,name='in'))
    n = rf.Network(name="guess", s=calc_circle(guess,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='black')
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(cl, rl)), real=np.real(calc_circle(cl, rl)),marker_color="green",showlegend=True,name='guess'))
    n = rf.Network(name="gama_opt", s=calc_circle(gamma_opt,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='brown')
    i=rf.Network(name="in", s=calc_circle(Cs, Rs))
    o=rf.Network(name="out", s=calc_circle(Cl, Rl))
    crash,O,I,index=[],np.round(o.s_db,4),np.round(i.s_db,4),0
    i,o=i.s,o.s
    for point1 in I:
        for point2 in O:
            if point2==point1:
                #if abs(R2P(i[index][0][0])[0]-R2P(o[index][0][0])[0])<0.07:
                    crash.append(i[index][0][0])
        index=index+1 
    index=0    
    for i in crash:
        n = rf.Network(name="c"+str(index), s=calc_circle(i,0.02))
        n.plot_s_smith(lw=1,draw_labels=True)
        index=index+1
    return crash    
    #plotlist.append(go.Scattersmith(imag=[gamma_opt.imag], real=[gamma_opt.real],marker_color="brown",showlegend=True,name="gama_opt"))
    #fig = go.Figure(plotlist)
    #fig.update_layout(height=800,width=1000)
    #fig.show()   
# we need the normalized equivalent noise and optimum source coefficient to calculate the constant noise circles
def noisecircle(rn,gamma_opt,fmin,noise,guess, z0):
    #plotlist=[]
    N = ((FDB10(noise)  - FDB10(fmin))/(4*rn/z0))*(abs(1+gamma_opt)**2)
    c_n = gamma_opt/(1+N)
    r_n = (1/(1-N))*np.sqrt(N**2 + N-N*(abs(gamma_opt)**2))

    n = rf.Network(name=str(noise), s=calc_circle(c_n, r_n))
    plt.figure(figsize=(10, 8), dpi=80)
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(c_n, r_n)), real=np.real(calc_circle(c_n, r_n)),marker_color="purple",showlegend=True,name='Noise:'+str(noise)))
    n = rf.Network(name="guess", s=calc_circle(guess,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='black')
    n = rf.Network(name="gama_opt", s=calc_circle(gamma_opt,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='brown')
    #plotlist.append(go.Scattersmith(imag=[gamma_opt.imag], real=[gamma_opt.real],marker_color="brown",showlegend=True,name="gama_opt"))
    print("the optimum source reflection coefficient is ", gamma_opt)    
    #fig = go.Figure(plotlist)
    #fig.update_layout(height=800,width=1000)
    #fig.show()
def gain_noise(gama_s,noise, gamma_opt, fmin, rn, z0, s11):
    sqabs = lambda x: np.square(np.absolute(x))
    #plotlist=[]
    N = ((FDB10(noise)  - FDB10(fmin))/(4*rn/z0))*(abs(1+gamma_opt)**2)
    c_n = gamma_opt/(1+N)
    r_n = (1/(1-N))*np.sqrt(N**2 + N-N*(abs(gamma_opt)**2))
    gs=(( 1-sqabs(gama_s) / sqabs(1-s11*gama_s)))*(1-sqabs(s11))
    Cs = (gs*np.conjugate(s11))/(1-(1-gs)*sqabs(s11))
    Rs = (np.sqrt(1-gs)*(1-sqabs(s11)))/(1-(1-gs)*sqabs(s11))
    plt.figure(figsize=(10, 8), dpi=80)
    n = rf.Network(name=str(noise), s=calc_circle(c_n, r_n))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(c_n, r_n)), real=np.real(calc_circle(c_n, r_n)),marker_color="purple",showlegend=True,name='Noise:'+str(noise)))
    n = rf.Network(name="in", s=calc_circle(Cs, Rs))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(Cs, Rs)), real=np.real(calc_circle(Cs, Rs)),marker_color="red",showlegend=True,name='in'))
    n = rf.Network(name="gama_opt", s=calc_circle(gamma_opt,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='brown')
    i=rf.Network(name="in", s=calc_circle(Cs, Rs))
    o=rf.Network(name=str(noise), s=calc_circle(c_n, r_n))
    crash,O,I,index=[],np.round(o.s_db,4),np.round(i.s_db,4),0
    i,o=i.s,o.s
    for point1 in I:
        for point2 in O:
            if point2==point1:
                #if abs(R2P(i[index][0][0])[0]-R2P(o[index][0][0])[0])<0.05:
                    crash.append(i[index][0][0])
        index=index+1 
    index=0    
    for i in crash:
        n = rf.Network(name="c"+str(index), s=calc_circle(i,0.02))
        n.plot_s_smith(lw=1,draw_labels=True)
        index=index+1
    return crash  
    #plotlist.append(go.Scattersmith(imag=[gamma_opt.imag], real=[gamma_opt.real],marker_color="brown",showlegend=True,name="gama_opt"))
    #fig = go.Figure(plotlist)
    #fig.update_layout(height=800,width=1000)
    #fig.show()  
def summary(gama_s,gama_l,noise,zs,zl,s11,s12,s21,s22,guess,gamma_opt,rn,z0,fmin):
    sqabs = lambda x: np.square(np.absolute(x))
    plotlist=[]
    delta = (s11*s22) - (s12*s21)
    rl = np.absolute((s12 * s21)/(sqabs(s22) - sqabs(delta)))
    cl = np.conj(s22 - delta*np.conj(s11))/(sqabs(s22) - sqabs(delta))

    rs = np.absolute((s12 * s21)/(sqabs(s11) - sqabs(delta)))
    cs = np.conj(s11 - delta*np.conj(s22))/(sqabs(s11) - sqabs(delta))

    plt.figure(figsize=(10, 8), dpi=80)
    n = rf.Network(name="out", s=calc_circle(cl, rl))
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(cl, rl)), real=np.real(calc_circle(cl, rl)),marker_color="green",showlegend=True,name='out'))
    n.plot_s_smith(lw=2,draw_labels=True)
    n = rf.Network(name="in", s=calc_circle(cs, rs))
    n.plot_s_smith(lw=2,draw_labels=True)  
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(cs, rs)), real=np.real(calc_circle(cs, rs)),marker_color="red",showlegend=True,name='in'))
    n = rf.Network(name="guess", s=calc_circle(guess,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='black')
    gs=(( 1-sqabs(gama_s) / sqabs(1-s11*gama_s)))*(1-sqabs(s11))
    gl=(( 1-sqabs(gama_l) / sqabs(1-s22*gama_l)))*(1-sqabs(s22))

    Cs = (gs*np.conjugate(s11))/(1-(1-gs)*sqabs(s11))
    Rs = (np.sqrt(1-gs)*(1-sqabs(s11)))/(1-(1-gs)*sqabs(s11))

    Cl = (gl*np.conjugate(s22))/(1-(1-gl)*sqabs(s22))
    Rl = (np.sqrt(1-gl)*(1-sqabs(s22)))/(1-(1-gl)*sqabs(s22))

    n = rf.Network(name="gain out", s=calc_circle(Cl, Rl))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(Cl, Rl)), real=np.real(calc_circle(Cl, Rl)),showlegend=True,name="gain out"))
    n = rf.Network(name="gain in", s=calc_circle(Cs, Rs))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(Cs, Rs)), real=np.real(calc_circle(Cs, Rs)),showlegend=True,name='gain in'))
    N = ((FDB10(noise)  - FDB10(fmin))/(4*rn/z0))*(abs(1+gamma_opt)**2)
    c_n = gamma_opt/(1+N)
    r_n = (1/(1-N))*np.sqrt(N**2 + N-N*(abs(gamma_opt)**2))
    n = rf.Network(name="noise "+str(noise), s=calc_circle(c_n, r_n))
    n.plot_s_smith(lw=2,draw_labels=True)
    #plotlist.append(go.Scattersmith(imag=np.imag(calc_circle(c_n, r_n)), real=np.real(calc_circle(c_n, r_n)),marker_color="purple",showlegend=True,name='Noise:'+str(noise)))
    n = rf.Network(name="gama_opt", s=calc_circle(gamma_opt,0.01))
    n.plot_s_smith(lw=1,draw_labels=True,color='brown')  
    n = rf.Network(name="Load", s=calc_circle(zl/z0,0.01))
    n.plot_s_smith(lw=1,draw_labels=True)  
    n = rf.Network(name="Source", s=calc_circle(zs/z0,0.01))
    n.plot_s_smith(lw=1,draw_labels=True)  
    #plotlist.append(go.Scattersmith(imag=[gamma_opt.imag], real=[gamma_opt.real],marker_color="brown",showlegend=True,name="gama_opt"))
    #fig = go.Figure(plotlist)
    #fig.update_layout(height=800,width=1000)
    #fig.show()  
def matching_network_LC_1(L, C,f360772a,z0):
    ' L and C in nH and pF'
    line = rf.media.DefinedGammaZ0(frequency=f360772a.frequency, z0=z0)
    return line.shunt_inductor(L*1e-9)**line.capacitor(C*1e-12)

def matching_network_CL_1(C,L,f360772a,z0):
    ' C and L in nH and pF'
    line = rf.media.DefinedGammaZ0(frequency=f360772a.frequency, z0=z0)
    return line.shunt_capacitor(C*1e-12)**line.inductor(L*1e-9)

def matching_network_LC_2(L, C,f360772a,z0):
    ' C and L in nH and pF'
    line = rf.media.DefinedGammaZ0(frequency=f360772a.frequency, z0=z0)
    return line.inductor(L*1e-9)**line.shunt_capacitor(C*1e-12)

def matching_network_CL_2(C,L,f360772a,z0):
    ' L and C in nH and pF'
    line = rf.media.DefinedGammaZ0(frequency=f360772a.frequency, z0=z0)
    return line.capacitor(C*1e-12)**line.shunt_inductor(L*1e-9)

def plotscheme(zs,zl,rin,rout,inputpkg,outputpkg,z0):
    rs,rl=zs.real,zl.real
    _,res,jx1,jx2,xskind,xpkind,restype,absorb=inputpkg
    res=res*1e+9
    d  = schemdraw.Drawing()
    d += elm.Ground()
    d += (V1 := elm.SourceV())
    d += elm.Resistor().right().label(str(z0)+'Ohm')
    if absorb=='s':
        if restype=='inductor':
            d.push()
            d += elm.Inductor().down().label(str(round(res,3))+'μH')
            d += elm.Ground()
            d.pop()
        else:
            d.push()
            d += elm.Capacitor().down().label(str(round(res,3))+'pF')
            d += elm.Ground()
            d.pop()     
    if rs>rl:
        if xskind=='Capacitor':
            jx1,jx2=jx1*1e+12,jx2*1e+9
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(jx1,3))+'pF')
            d += elm.Ground()
            d.pop()
            d += elm.Inductor().right().label(str(round(jx2,3))+'μH')
        else:
            jx1,jx2=jx1*1e+9,jx2*1e+12
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(jx1,3))+'μH')
            d += elm.Ground()
            d.pop()
            d += elm.Capacitor().right().label(str(round(jx2,3))+'pF')         
    else:
        if xskind=='Capacitor':
            jx1,jx2=jx1*1e+12,jx2*1e+9
            d += elm.Capacitor().label(str(round(jx1,3))+'pF')
            d.push()
            d += elm.Inductor().down().label(str(round(jx2,3))+'μH')
            d += elm.Ground()
            d.pop()
        else:
            jx1,jx2=jx1*1e+9,jx2*1e+12
            d += elm.Inductor().label(str(round(jx1,3))+'μH')
            d.push()
            d += elm.Capacitor().down().label(str(round(jx2,3))+'pF')
            d += elm.Ground()
            d.pop()
    if absorb=='l':
        if restype=='inductor':
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(res,3))+'μH')
            d += elm.Ground()
            d.pop()
        else:
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(res,3))+'pF')
            d += elm.Ground()
            d.pop()          
    d += elm.Resistor().right().label(str(rin)+'Ohm')
    d += (Q1 := elm.Bjt())
    d += elm.Ground().at((Q1, 'emitter'))
    d += elm.Resistor().right().at((Q1, 'collector')).label(str(round(rout,3))+'Ohm')
    _,res,jx1,jx2,xskind,xpkind,restype,absorb=outputpkg
    res=res*1e+9
    if absorb=='s':
        if restype=='inductor':
            d.push()
            d += elm.Inductor().down().label(str(round(res,3))+'μH')
            d += elm.Ground()
            d.pop()
        else:
            d.push()
            d += elm.Capacitor().down().label(str(round(res,3))+'pF')
            d += elm.Ground()
            d.pop()
    if rs>rl:
        if xskind=='Capacitor':
            jx1,jx2=jx1*1e+12,jx2*1e+9
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(jx1,3))+'pF')
            d += elm.Ground()
            d.pop()
            d += elm.Inductor().right().label(str(round(jx2,3))+'μH')
        else:
            jx1,jx2=jx1*1e+9,jx2*1e+12
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(jx1,3))+'μH')
            d += elm.Ground()
            d.pop()
            d += elm.Capacitor().right().label(str(round(jx2,3))+'pF')            
    else:
        if xskind=='Capacitor':
            jx1,jx2=jx1*1e+12,jx2*1e+9
            d += elm.Capacitor().label(str(round(jx1,3))+'pF')
            d.push()
            d += elm.Inductor().down().label(str(round(jx2,3))+'μH')
            d += elm.Ground()
            d.pop()
        else:
            jx1,jx2=jx1*1e+9,jx2*1e+12
            d += elm.Inductor().label(str(round(jx1,3))+'μH')
            d.push()
            d += elm.Capacitor().down().label(str(round(jx2,3))+'pF')
            d += elm.Ground()
            d.pop()
    if absorb=='l':
        if restype=='inductor':
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(res,3))+'μH')
            d += elm.Ground()
            d.pop()
        else:
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(res,3))+'pF')
            d += elm.Ground()
            d.pop()        
    d+=elm.Line()
    d += elm.Resistor().down().label(str(z0)+'Ohm')
    d += elm.Ground()
    return d
    
"""
ω =2πf, 
X =the reactance as read from the chart, 
B =the susceptance as read from the chart, 
N=the number used to normalize the original impedances that are to be matched.
"""
def series_C_component(freq,x,z0):
  w=(2*np.pi*freq)
  return 1/(w*x*z0)

def series_L_component(freq,x,z0):
  w=(2*np.pi*freq)
  return (x*z0)/w

def shunt_C_component(freq,b,z0):
  w=(2*np.pi*freq)
  return b/(w*z0)

def shunt_L_component(freq,b,z0):
  w=(2*np.pi*freq)
  return z0/(w*b) 

def L_match(zs,zl,freq,f360772a,z0):
    rs,rl=zs.real,zl.real
    line = rf.media.DefinedGammaZ0(frequency=f360772a.frequency, z0=z0)
    if zs.imag==0:
        absorb='l'
        if zl.imag<0:
            restype='inductor'
            Res=shunt_L_component(freq,-zl.imag,z0)
            res=line.shunt_inductor(Res*1e-9)
        else:    
            restype='capacitor'
            Res=shunt_C_component(freq,-zl.imag,z0)
            res=line.shunt_capacitor(Res*1e-12)
    else:   
        absorb='s'
        if zs.imag<0:
            restype='inductor'
            Res=shunt_L_component(freq,-zs.imag,z0)
            res=line.shunt_inductor(Res*1e-9)
        else:    
            restype='capacitor'
            Res=shunt_C_component(freq,-zs.imag,z0)
            res=line.shunt_capacitor(Res*1e-12)        
    if rs>rl:
        m=rs/rl
        q=np.sqrt(m-1)
        xs=q*rl
        xss=xs*(1+1/(q**2))
        xp=-xss
    else:
        m=rl/rs
        q=np.sqrt(m-1)
        xp=rl/q
        xpp=xp/(1+1/(q**2))
        xs=-xpp
    if xs>0:
        xskind='Inductor'
    else:
        xskind='Capacitor'
    if xp>0:
        xpkind='Inductor'
    else:
        xpkind='Capacitor'
    xs,xp,Res=abs(xs),abs(xp),abs(Res)    
    if rs>rl:
        if xskind=='Inductor':
            jx1=shunt_L_component(freq,xs,z0)
            jx2=series_C_component(freq,xp,z0)
            m=matching_network_LC_1(jx1,jx2,f360772a,z0)
        else:    
            jx1=shunt_C_component(freq,xs,z0)  
            jx2=series_L_component(freq,xp,z0)
            m=matching_network_CL_1(jx1,jx2,f360772a,z0)
    else:
        if xskind=='Inductor':
            jx1=series_L_component(freq,xs,z0)
            jx2=shunt_C_component(freq,xp,z0)
            m=matching_network_LC_2(jx1,jx2,f360772a,z0)
        else:  
            jx1=series_C_component(freq,xs,z0)
            jx2=shunt_L_component(freq,xp,z0)  
            m=matching_network_CL_2(jx1,jx2,f360772a,z0)
    if absorb=='s':
        m=res**m
    else:
        m=m**res
    print('jx1',jx1,xskind,'jx2',jx2,xpkind,'restype',restype,'res',Res,'Q:',q)   
    return [m,Res,jx1,jx2,xskind,xpkind,restype,absorb]
