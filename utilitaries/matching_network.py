import numpy as np
import skrf as rf

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
