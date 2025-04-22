import numpy as np

def get_phase(u,pu,betu,alfu,dispu=0,dispuprime=0,delta=0):
    # get phase space amplitude and angle
    u = u-dispu*delta
    pu = pu-dispuprime*delta
    un = u/np.sqrt(betu)
    pun = u*alfu/np.sqrt(betu) + pu*np.sqrt(betu)
    Ju = un**2 + pun**2
    angle = np.angle(un-1j*pun)
    return Ju, angle

def get_particle_tunes(u0, pu0, u1, pu1, 
                       betu, alfu, 
                       dispu=0,dispuprime=0,
                       delta0=0, delta1=0):
    # get particle tunes from phase space coordinates
    Ju0, angle0 = get_phase(u0, pu0, betu, alfu, dispu, dispuprime, delta0)
    Ju1, angle1 = get_phase(u1, pu1, betu, alfu, dispu, dispuprime, delta1)
    qu = (angle1 - angle0)/(2*np.pi)
    #qu[np.where(qu<0)] += 1
    return qu