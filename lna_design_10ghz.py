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


from utilitaries.basic_calculations import *
from utilitaries.circuit_drawing import *
from utilitaries.gain_calculations import *
from utilitaries.matching_network import *
from utilitaries.noise_calculations import *
from utilitaries.smith_chart_utils import *
from utilitaries.stability_analysis import *

#from utilitaries.functions import *

f360772a = rf.Network('f360772a.s2p')
plt.figure(figsize=(10, 10), dpi=80)
f360772a.plot_it_all()
begin=f360772a

freq=10 #GHz
noisefig=0.5 #dB
gamma_opt = P2R(0.60,129)
fmin = 0.44 #dB
z0=f360772a.z0[freq][0].real
rn = 0.05*z0 #noise resistance
s11,s12=f360772a.s[freq].tolist()[0]
s21,s22=f360772a.s[freq].tolist()[1]
print(s11,s12,s21,s22)


rollet(s11,s12,s21,s22)
_=mou(s11,s12,s21,s22)

R=[0.83,0.67]
plotstablecircules(s11,s12,s21,s22,R)
#read it from chart
rin=(z0*0.4).real
rout=(z0*0.22).real
print("for unconditionally stable we need rin:",rin,"Ohm","rout:",rout,"Ohm")


line = rf.media.DefinedGammaZ0(frequency=f360772a.frequency, z0=50)
f360772a=line.resistor(rin) ** f360772a ** line.resistor(rout)
plt.figure(figsize=(10, 10), dpi=80)
f360772a.plot_it_all()


s11,s12=f360772a.s[freq].tolist()[0]
s21,s22=f360772a.s[freq].tolist()[1]
rollet(s11,s12,s21,s22)
mou(s11,s12,s21,s22)


print("Maximum Available Gain :",MAG(s11,s12,s21,s22),"dB")
R=[0]
plotstablecircules(s11,s12,s21,s22,R)




fewquencylist=begin.frequency.f
fewquencylist=fewquencylist/(10**9)
moulist,moulist2=[],[]

s11list=begin.s11.s
s12list=begin.s12.s
s21list=begin.s21.s
s22list=begin.s22.s
for i in range(len(fewquencylist)):   
     moulist.append(mou(s11list[i],s12list[i],s21list[i],s22list[i],verbus=False))
moulist=np.concatenate( moulist, axis=0 )     

s11list=f360772a.s11.s
s12list=f360772a.s12.s
s21list=f360772a.s21.s
s22list=f360772a.s22.s
for i in range(len(fewquencylist)):   
     moulist2.append(mou(s11list[i],s12list[i],s21list[i],s22list[i],verbus=False))
moulist2=np.concatenate( moulist2, axis=0 )  

plt.figure(figsize=(10, 8), dpi=80)
plt.grid(visible=True, which='major', axis='both')
plt.plot(fewquencylist,np.ones(len(fewquencylist)),label=' Stable line')
plt.plot(fewquencylist, moulist,label=' Transistor only')
plt.plot(fewquencylist, moulist2,label=' Uncondiotionaly stable',color='red')
plt.title(str(begin.frequency.start/10**9)+' To '+str(begin.frequency.stop/10**9)+' Ghz')
plt.xlabel('Frequency GHz')
plt.ylabel('MOU')
plt.legend(loc=4)
#specify axis tick step sizes
plt.xticks(np.arange(min(fewquencylist), max(fewquencylist)+1, 1))
plt.show()



guess=(0.3309+0.5188j)
noisecircle(rn,gamma_opt,fmin,noisefig,guess,z0)


#zs=55-53j
zs=55+5j
gs=gamas(zs, z0)
Volunteer=gain_noise(gs,noisefig, gamma_opt, fmin, rn, z0, s11)
guess=Volunteer[3]




#gl=0.01+0.02j
gl=gamal(guess, s11, s12, s21, s22)
zl=findzfromgama(gl, z0)[0]
print('gs:',gs,'gl:',gl,'\nzs:',zs,'zl:',zl)
C=gain_ciecles(s11,gs,s22,gl,guess/z0,gamma_opt)



z11=rf.network.s2z(f360772a.s11.s,z0)
z22=rf.network.s2z(f360772a.s22.s,z0)
if (z11!=z22).all() :
    print("network is reciprocal")    
else:
    print("network is symmetrical")
zl=rf.network.s2z(gl.reshape(1,1,1),z0)[0][0][0]   
zs=rf.network.s2z(np.array(guess).reshape(1,1,1),z0)[0][0][0] 
print('zs:',zs,'zl:',zl)



zout=-zl
zin=-zs
zin=findzin(zs)
zout=findzout(zl)
print("Zin is:",zin,"Zout is:",zout)



summary(gs,gl,noisefig,zs,zl,s11,s12,s21,s22,guess,gamma_opt,rn,z0,fmin)
print('Transducer Gain (Total):',GT(s11,s12,s21,s22,guess,gl),'dB','lossed gain:',MAG(s11,s12,s21,s22)-GT(s11,s12,s21,s22,guess,gl),"dB")
print('input SWR:',swr(gs),'output SWR:',swr(gl))
print('gama S:',R2P(guess),'gama L:',R2P(gl))


#https://www.mathworks.com/help/rf/ug/designing-matching-networks-part-1-networks-with-an-lna-and-lumped-elements.html
#https://scikit-rf.readthedocs.io/en/latest/examples/circuit/Lumped%20Element%20Circuits.html
#https://www.allaboutcircuits.com/tools/l-match-impedance-matching-circuits/
gs=guess


outputpkg =L_match(np.conj(zl),z0,freq*10**9,f360772a,z0)
plt.figure(figsize=(5, 4), dpi=80)
outputpkg[0].plot_s_mag(lw=2)
plt.ylim(bottom=0)
plt.show()



inputpkg =L_match(z0,np.conj(zs),freq*10**9,f360772a,z0)
plt.figure(figsize=(5, 4), dpi=80)
inputpkg[0].plot_s_mag(lw=2)
plt.ylim(bottom=0)
plt.show()


d=plotscheme(zs,zl,rin,rout,inputpkg,outputpkg,z0)
d.draw()



plt.figure(figsize=(12, 8), dpi=80)
plt.grid(visible=True, which='major', axis='both')
f360772a.s21.plot_s_db(label=' s21')
f360772a.s11.plot_s_db(label=' s11')
#specify axis tick step sizes
fewquencylist=f360772a.frequency.f
#fewquencylist=fewquencylist/(10**9)
plt.xticks(np.arange(min(fewquencylist), max(fewquencylist)+1, 10**9/2))
plt.yticks(np.arange(min(f360772a.s11.s_db), max(f360772a.s21.s_db)+1, 1))
plt.title('At '+str(freq)+' GHz s11:'+str(DB20(f360772a.s[freq].tolist()[0][0]))+' s21: '+str(DB20(f360772a.s[freq].tolist()[1][0])))
plt.show()



plt.figure(figsize=(12, 10), dpi=80)
f360772a.plot_it_all()


plt.figure(figsize=(10, 10), dpi=80)
amplifier = f360772a ** outputpkg[0]
amplifier.plot_it_all()


plt.figure(figsize=(10, 10), dpi=80) 
amplifier = inputpkg[0] ** f360772a 
amplifier.plot_it_all()



amplifier = inputpkg[0] ** f360772a ** outputpkg[0]
plt.figure(figsize=(12, 10), dpi=80)
plt.grid(visible=True, which='major', axis='both')
#amplifier.s21.plot_s_db(label=' s21')
#amplifier.s11.plot_s_db(label=' s11')
amplifier.plot_it_all()


s11,s12=f360772a.s[freq].tolist()[0]
s21,s22=f360772a.s[freq].tolist()[1]
rollet(s11,s12,s21,s22)
_=mou(s11,s12,s21,s22)



print('Transducer Gain (Total):',GT(s11,s12,s21,s22,guess,gl),'dB','lossed gain:',(MAG(s11,s12,s21,s22)-GT(s11,s12,s21,s22,guess,gl)).real,"dB")



import pandas as pd
Temp=[15,20,25,30,35,40,45,50,55,60] # celcious
snr=4 # dB C/N
freq_step=f360772a.frequency.f[1]-f360772a.frequency.f[0]
s21_list=f360772a.s21.s
bandwidthl,ktbl,noisefloorl,nfl,sensl=[],[],[],[],[]
print('input SNR {} dB output SNR {} db'.format(snr,snr/FDB10(noisefig)))
for T in Temp:
    nfl.append(noisefig)
    bandwidth=Bandwidth(s21_list,freq_step)[0]
    ktb=KTB(T,bandwidth)[0]
    noisefloor=Thermal_noise_floor(T,bandwidth).real
    sens=sensitivity(T,bandwidth,noisefig,snr).real
    bandwidthl.append(bandwidth[0]/(10**9))
    ktbl.append(ktb)
    noisefloorl.append(noisefloor)
    sensl.append(sens)
data={
    'Temp C':Temp,
    'Bandwidth GHz':bandwidthl,
    'KTB':ktbl,
    'Thermal_noise_floor dBm':noisefloorl,
    'NF dB':nfl,
    'Sensitivity dBm':sensl
}
pd.DataFrame.from_dict(data)


freq = f360772a .frequency
port1 = rf.Circuit.Port(frequency=freq, name='input_port', z0=50)
port2 = rf.Circuit.Port(frequency=freq, name='output_port', z0=50)
line = rf.media.DefinedGammaZ0(frequency=freq, z0=50)
cap = line.capacitor(70e-12, name='cap')
ind = line.inductor(24e-9, name='ind')
res_series = line.resistor(1e-2, name='res_series')
res_parallel = line.resistor(20e6, name='res_parallel')
cap_shunt = line.capacitor(50e-12, name='cap_shunt')
ground = rf.Circuit.Ground(frequency=freq, name='ground', z0=50)

connections = [
    [(port1, 0), (res_series, 0), (res_parallel, 0)],
    [(res_series, 1), (cap, 0)],
    [(cap, 1), (ind, 0)],
    [(ind, 1), (cap_shunt, 0), (res_parallel, 1), (port2, 0)],
    [(cap_shunt, 1), (ground, 0)],
]
circuit = rf.Circuit(connections)
LCC_from_circuit = circuit.network
circuit.plot_graph(network_labels=True, edge_labels=True, port_labels=True,inter_labels=True,network_fontsize=20,edges_fontsize=20,port_fontsize=20)


power = [1, 0]  # 1 Watt at port1 and 0 at port2
phase = [0, 0]  # 0 radians
V_at_ports = circuit.voltages_external(power, phase)
print(V_at_ports)
I_at_ports = circuit.currents_external(power, phase)
print(I_at_ports)