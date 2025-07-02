import numpy as np
from .basic_calculations import calc_circle
from matplotlib import pyplot as plt
import skrf as rf

def GT(s11,s12,s21,s22,gs,gl):
    # Transducer Gain : 
    # the actual gain of an amplifier stage including the effects of input and output matching and device gain
    qabs = lambda x: np.square(np.absolute(x))
    return (qabs(s21)*(1-qabs(gs))*(1-qabs(gl)))/qabs((1-s11*gs)*(1-s22*gl)-s12*s21*gl*gs)

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