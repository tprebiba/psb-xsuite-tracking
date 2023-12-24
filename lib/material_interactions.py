# Based on https://github.com/PyORBIT-Collaboration/py-orbit/blob/master/src/orbit/MaterialInteractions/MaterialInteractions.cc
# C++ module written in Python
# Adapted for Xsuite 

import math
import random
import scipy.constants as cnst
from scipy.special import j0,j1


def mcsJackson(stepsize, z, a, rho, beta, pfac, particles, ip):
    """
    Multiple coulomb scatter a particle through a set distance using a simple formulation from JD Jackson.
    Units are cgs, esu.
    
    PARAMETERS:
    - stepsize: length (thickness) of scattering material {mm}.
    - z: atomic number of scattering material.
    - a: atomic weight of scattering material.
    - rho: mass density of scattering material.
    - beta: v/c.
    - pfac: 1+dp/p0.
    - particles: xpart.Particles() object
    - ip: particle index
    
    RETURNS: 
        Nothing
    """
    
    # Convert to mm and mrad for this routine. And density to g/cm3.
    x = particles.x[ip]*1000.0
    y = particles.y[ip]*1000.0
    px = particles.px[ip]*1000.0
    py = particles.py[ip]*1000.0
    stepsize *= 1000.0
    rho /= 1000.0
    
    pi = math.pi
    cvel = 2.99792458e+10 # cm/s
    eesu = 4.8032068e-10 # esu
    hbar = 1.05457e-34 # J s
    AMU = 1.66053906660e-27 # kg
    aBohr = 5.2917721067e-11 # m
    AP = 1.00727647 # atomic mass of proton
    gamma = 1.0 / math.sqrt(1.0 - beta * beta)
    
    aMax = 1.4 * aBohr / pow(z, 0.333333)
    aMin = 1.4e-13 * pow(a, 0.333333)
    
    N = rho / (a * AMU)
    
    T = stepsize / 10.0
    
    pmom = gamma * beta * cvel * AP * AMU
    
    thMin = hbar / (pmom * aMax)
    thMax = hbar / (pmom * aMin)
    
    thMin2i = 1.0 / (thMin * thMin)
    thMax2i = 1.0 / (thMax * thMax)
    th2iDiff = thMin2i - thMax2i
    
    th2s = 2.0 * math.log(thMax / thMin) / th2iDiff
    coeff = 2.0 * z * eesu * eesu / (pmom * beta * cvel)
    
    sigTot = pi * coeff * coeff * th2iDiff
    nColl = N * T * sigTot
    
    th2Tot = nColl * th2s
    
    probrp = random.random()
    probxy = 2.0 * pi * random.random()
    
    angle = math.sqrt(-th2Tot * math.log(probrp))
    anglexMCS = angle * math.cos(probxy)
    angleyMCS = angle * math.sin(probxy)
    
    xpfac = px / (1 * pfac) # 1000.0 -> 1
    ypfac = py / (1 * pfac) # 1000.0 -> 1
    
    anglex = math.atan(xpfac) + anglexMCS
    angley = math.atan(ypfac) + angleyMCS

    tanglex = math.tan(anglex)
    tangley = math.tan(angley)
    
    # to be checked again
    px = tanglex * (1 * pfac) # 1000.0 -> 1
    py = tangley * (1 * pfac) # 1000.0 -> 1

    directionfac = math.sqrt(1.0 + xpfac * xpfac + ypfac * ypfac)
    zstep = stepsize / directionfac
    x += zstep * (xpfac + tanglex) / 2.0
    y += zstep * (ypfac + tangley) / 2.0

    # Convert back to m and rad.
    x /= 1000.0
    y /= 1000.0
    px /= 1000.0
    py /= 1000.0
    particles.x[ip] = x
    particles.y[ip] = y
    particles.px[ip] = px
    particles.py[ip] = py


def ruthScattJackson(stepsize, z, a, rho, beta, trackit, pfac, thx, thy):
    """
    Rutherford Scattering from an atomic nucleus.
    Simple formulation from JD Jackson:
    Classical Electrodynamics, Chapter 13.
    Units are cgs, esu.
    
    PARAMETERS:
        stepsize: length (thickness) of scattering material {mm}.
        Z: atomic number of scattering material.
        A: atomic weight of scattering material.
        rho: mass density of scattering material.
        beta: v/c.
        trackit: 
        pfac: 1+dp/p0.
    
    RETURNS: 
        A double array with thetax and thetay. It also calculates and assigns the 
        Rutherford cross section rcross. 
    """
    stepsize *= 1000.0 # Convert to mm
    rho /= 1000.0 # Convert to g/cm3

    pi = math.pi
    cvel = 2.99792458e+08 * 100.0
    eesu = 4.80320451e-10
    hbar = 1.05443e-34 / (1.602176634e-19 * 1e+9)
    AMU = 1.66053904e-27
    aBohr = 5.29177210903e-11
    AP = 1.00727647
    gamma = 1.0 / math.sqrt(1.0 - beta**2)

    aMax = 1.4 * aBohr / math.pow(z, 0.333333)
    aMin = 1.4e-13 * math.pow(a, 0.333333)

    N = rho / (a * AMU)
    T = stepsize / 10.0

    pmom = gamma * beta * cvel * AP * AMU

    thMin = hbar / (pmom * aMax)
    thMax = hbar / (pmom * aMin)

    thMin2i = 1.0 / (thMin * thMin)
    thMax2i = 1.0 / (thMax * thMax)
    th2iDiff = thMin2i - thMax2i

    th2s = 2.0 * math.log(thMax / thMin) / th2iDiff
    coeff = 2.0 * z * eesu * eesu / (pmom * beta * cvel)

    sigTot = pi * coeff * coeff * th2iDiff

    nColl = N * T * sigTot

    th2Tot = nColl * th2s

    if(thMin < (2.0 * math.sqrt(th2Tot))): 
        thMin = 2.0 * math.sqrt(th2Tot)

    thMin2i = 1.0 / (thMin * thMin)
    th2iDiff = thMin2i - thMax2i
    rcross = 1.e+24 * pi * coeff * coeff * th2iDiff

    if(rcross < 0.0): 
        rcross = 0.0

    thMin *= 1.4
    thMin2i = 1.0 / (thMin * thMin)
    th2iDiff = thMin2i - thMax2i

    if(trackit != 0):
        thx = 0.0
        thy = 0.0

        if(thMin < thMax):
            probrp = random.uniform(0.0, 1.0)
            probxy = 2.0 * pi * random.uniform(0.0, 1.0)

            denom2 = probrp * th2iDiff + thMax2i
            th = math.sqrt(1.0 / denom2)

            thx = th * math.cos(probxy)
            thy = th * math.sin(probxy)

    return rcross, thx, thy


def momentumKick(t, p, dpx, dpy):
    """
    Generates random, 2D projection of momentum kick.  See Monte 
    Carlo Methods in PDG, or the K2 code description in the PhD
    thesis of Nuria Catalan-Lasheras.
    
    PARAMETERS:
        t: magnitude of momentum transfer
        p: particle momentum
    
    RETURNS: 
        dp, the momentum kick generated in each plane.
    """
    dp = [0.0, 0.0]
    r2 = 10.0
    theta = math.acos(1 - t / (2*p**2))

    while r2 > 1.0:
        randomvalue = random.random()
        va = 2 * randomvalue - 1
        vb = randomvalue - 1
        va2 = va * va
        vb2 = vb * vb
        r2 = va2 + vb2

    dp[0] = theta * (2*va*vb)/r2
    dp[1] = theta * (va2 - vb2)/r2
    dpx = dp[0]
    dpy = dp[1]
    return dpx, dpy


def elastic_t(p, a):
    """
    Generates random momentum angle for low energy elastic scattering.
    
    PARAMETERS:
        p: particle momentum
        a: nuclear mass number
    
    RETURNS: 
        double. 
    """
    c = cnst.c
    PI = math.pi

    m = .938 # proton, GeV/c^2
    h = 4.135e-15/1e9
    M_nuc = a*.931494 # GeV/c^2
    E = math.sqrt((p*1e-9)**2+m**2)
    E_cm = math.sqrt(m**2 + M_nuc**2 + 2*E*M_nuc) # See PDG.
    p_cm = (p*1e-9)*M_nuc/E_cm
    lamda = h*c/p_cm
    R = 1.4e-15 * pow(a, 1. / 3.) + lamda
    u = 2 * R /lamda * math.sin(PI/2)

    cnorm = 1. / 2. * (math.sqrt(PI) - math.sqrt(PI) * pow(j0(u), 2)
                        - math.sqrt(PI) * pow(j1(u), 2)) / math.sqrt(PI)

    random_value = random.random()

    theta = 0.0001
    angle_cm = 0
    t = 0
    found = 0
    while theta <= PI:
        x = 2. * R * math.sin(theta / 2.) / lamda
        f_x = 1. / 2. / cnorm * (math.sqrt(PI) - math.sqrt(PI) * pow(j0(x), 2)
                                    - math.sqrt(PI) * pow(j1(x), 2)) / math.sqrt(PI) - random_value
        if -0.01 < f_x < 0.01:
            angle_cm = theta
            theta = PI + 1.
            found = 1
        theta += 1.768e-3

    if found == 0:
        print("Warning, never found elastic t.\n")
    t = 2. * p * a * .931494 / M_nuc * pow(math.sin(angle_cm / 2.), 2)
    return t


def ionEnergyLoss(beta, z, a):
    """
    Compute ionization energy loss for a proton.  See Bethe-Bloch
    equation in the Particle Data Book (eq. 34.5)
    
    PARAMETERS:
        beta: relativistic beta.
        z: material charge number.
        a: atomic number.
    
    RETURNS: 
        double. Rate of energy (GeV) loss per stepsize where stepsize is a measure of 
        mass per unit area (kg/m^2). See PDG handbook. Units are GeV/(kg/m^2)
    """
    m_e = 0.511 # in MeV
    M = 938.27 # in MeV
    I = z*10.0
    
    gamma = 1.0/math.sqrt(1.0-beta**2)
    T_max = 1e6 * (2*m_e*beta**2*gamma**2)/(1.0 + 2.0*gamma*m_e/M + pow(m_e/M, 2.0)) # in eV
    arg = 2.0*(m_e*1e6)*beta**2*gamma**2*T_max/(I**2) # unitless
    K = 0.307075 # MeV*cm^2/mol
    dE = K*z/a/(beta**2)*(math.log(arg)/2.0 - beta**2) # MeV/(g/cm^2)
    
    # Unit conversion from MeV/(g/cm^2) to GeV/(kg/m^2)
    dE /= 10000.0
    
    return dE