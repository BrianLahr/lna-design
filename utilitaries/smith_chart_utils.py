import numpy as np
import skrf as rf
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from .basic_calculations import R2P, FDB10

def smithcircule(R):
    fig = go.Figure(go.Scattersmith(imag=np.imag(R), real=np.real(R),marker_color="red",showlegend=True,name='guess'))
    fig.update_layout(height=800,width=1000)
    fig.show()
    
def calc_circle(c, r):
    theta = np.linspace(0, 2 * np.pi, 1000)
    return c + r * np.exp(1.0j * theta)

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