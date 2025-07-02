import numpy as np
import skrf as rf
from matplotlib import pyplot as plt
from .basic_calculations import FDB10, DBm10, Celsius_to_Kelvin, calc_circle

def KTB(Temp,bandwidth):
    Temp=Celsius_to_Kelvin(Temp)
    return 1.38065 * 10**-23*Temp*bandwidth

def Thermal_noise_floor(Temp,bandwidth):
    return DBm10(KTB(Temp,bandwidth))

def sensitivity(Temp,bandwidth,noisefig,snr):
    # normaly negative
    Temp=Celsius_to_Kelvin(Temp)
    return DBm10(KTB(Temp,bandwidth))+noisefig+snr

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