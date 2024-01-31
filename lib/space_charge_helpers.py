import numpy as np

def get_phase(u,pu,betu,alfu):
    # get phase space amplitude and angle
    un = u/np.sqrt(betu)
    pun = u*alfu/np.sqrt(betu) + pu*np.sqrt(betu)
    Ju = un**2 + pun**2
    angle = np.angle(un-1j*pun)
    return Ju, angle

def get_particle_tunes(u0, pu0, u1, pu1, betu, alfu):
    # get particle tunes from phase space coordinates
    Ju0, angle0 = get_phase(u0, pu0, betu, alfu)
    Ju1, angle1 = get_phase(u1, pu1, betu, alfu)
    qu = (angle1 - angle0)/(2*np.pi)
    #qu[np.where(qu<0)] += 1
    return qu