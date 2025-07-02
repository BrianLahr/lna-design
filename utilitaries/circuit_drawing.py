import schemdraw
import schemdraw.elements as elm

def plotscheme(zs, zl, rin, rout, inputpkg, outputpkg, z0):
    rs, rl = zs.real, zl.real
    _, res, jx1, jx2, xskind, xpkind, restype, absorb = inputpkg
    res = res * 1e+9
    d = schemdraw.Drawing()
    d += elm.Ground()
    d += (V1 := elm.SourceV())
    d += elm.Resistor().right().label(str(z0) + 'Ohm')
    if absorb == 's':
        if restype == 'inductor':
            d.push()
            d += elm.Inductor().down().label(str(round(res, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
        else:
            d.push()
            d += elm.Capacitor().down().label(str(round(res, 3)) + 'pF')
            d += elm.Ground()
            d.pop()     
    if rs > rl:
        if xskind == 'Capacitor':
            jx1, jx2 = jx1 * 1e+12, jx2 * 1e+9
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(jx1, 3)) + 'pF')
            d += elm.Ground()
            d.pop()
            d += elm.Inductor().right().label(str(round(jx2, 3)) + 'μH')
        else:
            jx1, jx2 = jx1 * 1e+9, jx2 * 1e+12
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(jx1, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
            d += elm.Capacitor().right().label(str(round(jx2, 3)) + 'pF')         
    else:
        if xskind == 'Capacitor':
            jx1, jx2 = jx1 * 1e+12, jx2 * 1e+9
            d += elm.Capacitor().label(str(round(jx1, 3)) + 'pF')
            d.push()
            d += elm.Inductor().down().label(str(round(jx2, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
        else:
            jx1, jx2 = jx1 * 1e+9, jx2 * 1e+12
            d += elm.Inductor().label(str(round(jx1, 3)) + 'μH')
            d.push()
            d += elm.Capacitor().down().label(str(round(jx2, 3)) + 'pF')
            d += elm.Ground()
            d.pop()
    if absorb == 'l':
        if restype == 'inductor':
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(res, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
        else:
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(res, 3)) + 'pF')
            d += elm.Ground()
            d.pop()          
    d += elm.Resistor().right().label(str(rin) + 'Ohm')
    d += (Q1 := elm.Bjt())
    d += elm.Ground().at((Q1, 'emitter'))
    d += elm.Resistor().right().at((Q1, 'collector')).label(str(round(rout, 3)) + 'Ohm')
    _, res, jx1, jx2, xskind, xpkind, restype, absorb = outputpkg
    res = res * 1e+9
    if absorb == 's':
        if restype == 'inductor':
            d.push()
            d += elm.Inductor().down().label(str(round(res, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
        else:
            d.push()
            d += elm.Capacitor().down().label(str(round(res, 3)) + 'pF')
            d += elm.Ground()
            d.pop()
    if rs > rl:
        if xskind == 'Capacitor':
            jx1, jx2 = jx1 * 1e+12, jx2 * 1e+9
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(jx1, 3)) + 'pF')
            d += elm.Ground()
            d.pop()
            d += elm.Inductor().right().label(str(round(jx2, 3)) + 'μH')
        else:
            jx1, jx2 = jx1 * 1e+9, jx2 * 1e+12
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(jx1, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
            d += elm.Capacitor().right().label(str(round(jx2, 3)) + 'pF')            
    else:
        if xskind == 'Capacitor':
            jx1, jx2 = jx1 * 1e+12, jx2 * 1e+9
            d += elm.Capacitor().label(str(round(jx1, 3)) + 'pF')
            d.push()
            d += elm.Inductor().down().label(str(round(jx2, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
        else:
            jx1, jx2 = jx1 * 1e+9, jx2 * 1e+12
            d += elm.Inductor().label(str(round(jx1, 3)) + 'μH')
            d.push()
            d += elm.Capacitor().down().label(str(round(jx2, 3)) + 'pF')
            d += elm.Ground()
            d.pop()
    if absorb == 'l':
        if restype == 'inductor':
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(str(round(res, 3)) + 'μH')
            d += elm.Ground()
            d.pop()
        else:
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(str(round(res, 3)) + 'pF')
            d += elm.Ground()
            d.pop()        
    d += elm.Line()
    d += elm.Resistor().down().label(str(z0) + 'Ohm')
    d += elm.Ground()
    return d