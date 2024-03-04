# Based on https://github.com/PyORBIT-Collaboration/py-orbit/blob/master/src/orbit/MaterialInteractions/Foil.cc
# C++ module written in Python
# Adapted for Xsuite

import math
import random
import scipy.constants as cnst
#import cross_sections as cs
#import material_interactions as mi
import lib.cross_sections as cs
import lib.material_interactions as mi


class Foil():

    def __init__(self, xmin: float, xmax: float, ymin: float, ymax: float, thick: float, scatterChoice=0, activate_foil=0,
                 context = 'CPU'):
        """
        Foil class constructor.

        Parameters:
            xmin (float): Left edge of the foil in m
            xmax (float): Right edge of the foil in m
            ymin (float): Bottom edge of the foil in m
            ymax (float): Top edge of the foil in m
            thick (float): Thickness of the foil in ug/cm^2
            scatterChoice (int): The user choice of scattering routine, 
                                 0 for full scatter (default), 1 for simple scatter (no losses)
            activate_foil (int): 1 to activate foil, 0 to deactivate foil
            context (string): Context to use for the calculations (default = 'CPU')
        
        Returns:
            Nothing
        """
        self.xmin_ = xmin
        self.xmax_ = xmax
        self.ymin_ = ymin
        self.ymax_ = ymax
        self.thick_ = thick
        self.scatterChoice_ = scatterChoice
        self.activate_foil_ = activate_foil

        self.length_ = 0.0 # length of accelerator element (thin)
        self.ma_ = 0 # foil material index
        self.nHits = 0 # number of hits on the foil
        self.nLost = 0 # number of lost particles  

        if context == 'GPU':
            try:
                import cupy as cp
            except ImportError:
                print("# Foil : cupy module is required to run in GPU. ")
            self.np = cp
        elif context == 'CPU':
            try:
                import numpy as np
            except ImportError:
                print("# Foil : numpy module is required to run in CPU. ")
            self.np = np
        print('Foil object created in %s context.'%context)

    def setScatterChoice(self, choice):
        self.scatterChoice_ = choice

    def setActivateFoil(self, activate_foil):
        self.activate_foil_ = activate_foil

    def track(self, particles):
        """
		Tracks the particles 
		"""
        if self.activate_foil_ == 1:
            if(self.scatterChoice_ == 0):
                self.transverseFoilFullScatter(particles)
                particles.reorganize() # to put lost particles from the foil at the end of the array
            else:
                self.transverseFoilSimpleScatter(particles)
        else:
            #print('Foil is not activated.')
            pass

    def checkFoilFlag(self, x, y):
        """
        Filters the particles that fall within the foil. 

        Parameters:
            x (array): x coordinates of particles
            y (array): y coordinates of particles
        
        Returns:
            Returnes a mask of all particles that fall transversely within the foil.
        """
        mask_foil = (x >= self.xmin_) & (x <= self.xmax_) & (y >= self.ymin_) & (y <= self.ymax_)
        return mask_foil


    def transverseFoilSimpleScatter(self, particles):
        """
        Material scatters particles through a foil with simplified MCS scatter.

        Parameters:
            particles: particles object from xsuite
        
        Returns:
            Nothing
        """
        BohrRadius = 0.52917706e-8  # hydrogenic Bohr radius in cm
        hBar = 1.054571817e-27  # Planck's constant in erg*sec (1e-7 J*sec)
        echarge = 4.80320425e-10  # in esu or statcoulombs (cm^(3/2)*g^(1/2) in CGS or 3.3356e-10 C in SI)
        nAvogadro = 6.0221408e23 # Avogadro's number (1/mol)
        muScatter = 1.35 # ?
        emass = 9.1093837e-28 # electron mass in grams
        mass_proton = 0.938272013 # mass proton in GeV
        mass_proton_g = 1.67262e-27*1e3 # in gramms
        GeV2gramm = mass_proton_g/mass_proton # to convert GeV to gramms
        TFRadius = muScatter*BohrRadius*pow(cs.get_z(self.ma_), -0.33333) # Thomas-Fermi atom radius (cm):

        # only CPU for now.
        #context = particles._context
        #ctx2np = context.nparray_from_context_array

        # mask particles
        mask_alive = particles.state>=1 # alive particles
        mask_foil = self.checkFoilFlag(particles.x, particles.y) # within the foil particles
        mask = mask_alive & mask_foil
        nparticles = len(particles.x[mask])
        
        # Momentum in g*cm/sec
        pInj0 = particles.gamma0[mask] * (particles.mass0*1e-9*GeV2gramm) * particles.beta0[mask]*cnst.c*1e2
        
        # Minimum scattering angle:
        thetaScatMin = hBar/(pInj0*TFRadius)

        # Theta max as per Jackson (13.102)
        # Original script raises it to the 1/3 power, why not to the 1/2 ? Re-check Jackson
        thetaScatMax = 274.0*emass*100.0*cnst.c/(pInj0*pow(cs.get_a(self.ma_), 0.33333))
        
        # momentum times particle velocity (in units of g * cm^2 / sec^2)
        pv = pInj0 * particles.beta0[mask]*cnst.c*1.e2

        # Factor in parenthesis of (13.104)
        # z=1 (protons) while Z=cs.get_z(self.ma_)
        term = cs.get_z(self.ma_)*echarge**2/pv # in units statC^2 / (g*cm^2/sec^2) = cm
        
        # total scattering cross section as per Jackson (10.104) in units of cm^2
        sigmacoul = 4*math.pi*term**2/(thetaScatMin**2)

        # Scattering area per area (no units)
        nscatters = nAvogadro*(cs.get_rho(self.ma_)/1000.0)/cs.get_a(self.ma_)*self.thick_/(1.0e3*cs.get_rho(self.ma_))*sigmacoul

        # Mean free path (units of cm)
        lscatter = (self.thick_/(1.0e3*cs.get_rho(self.ma_)))/nscatters 

        # Distance remaining in foil in cm 
        zrl = self.np.full(nparticles, self.thick_ / (1.0e3 * cs.get_rho(self.ma_)))
        thetaX = self.np.zeros(nparticles)
        thetaY = self.np.zeros(nparticles)
        zrl_mask = self.np.where(zrl >= 0.0)[0]
        zrl_mask_len = len(zrl_mask)
        # Generate interaction points until particle exits foil
        while zrl_mask_len > 0:
            
            zrl[zrl_mask] += lscatter[zrl_mask] * self.np.log(self.np.random.uniform(0,1, size=zrl_mask_len))
            
            # Generate random angles
            phi = 2 * math.pi * self.np.random.uniform(0,1,size=zrl_mask_len)
            random1 = self.np.random.uniform(0,1,size=zrl_mask_len)
            theta = thetaScatMin[zrl_mask] * self.np.sqrt(random1 / (1.0 - random1))
            theta_mask = self.np.where(theta > (2.0 * thetaScatMax[zrl_mask]))[0]
            if len(theta_mask) > 0:
                theta[theta_mask] = 2.0 * thetaScatMax[zrl_mask] * self.np.ones(len(theta_mask))
            else:
                #print('No theta_mask')
                pass
            thetaX[zrl_mask] += theta * self.np.cos(phi)
            thetaY[zrl_mask] += theta * self.np.sin(phi)

            zrl_mask = self.np.where(zrl >= 0.0)[0]
            zrl_mask_len = len(zrl_mask)
        
        particles.px[mask] += thetaX
        particles.py[mask] += thetaY
        print('Number of particles hitting foil: ', nparticles)
        self.nHits += nparticles


    def transverseFoilFullScatter(self, particles):
        """
        Material scatters particles through a foil. MCS and nuclear scattering.

        Parameters:
            particles: particles object from xsuite
        
        Returns:
            Nothing
        """
        nAvogadro = 6.0221408e23 # Avogadro's number (1/mol)
        radlengthfac = 0.0 #?
        t = 0.0
        dp_x = 0.0 # for mi.momentumKick()
        dp_y = 0.0 # for mi.momentumKick()
        thx = 0.0
        thy = 0.0
        stepsize=0

        z = cs.get_z(self.ma_)
        a = cs.get_a(self.ma_)
        density = cs.get_rho(self.ma_) # in kg/m^3
        b_pN = 14.5 * pow(a, 0.6666667) # ?
        radlength = cs.get_radlength(self.ma_) # in m
        length = self.thick_ / (1.0e5 * density) # Foil length (in m)
        dlength = length * 1.0e-4

        nParts = len(particles.x)
        for ip in range(nParts):
        #for ip in range(1):
            #print(ip)
            zrl = length # remaining Foil length in the z direction (in m)
            foil_flag = self.checkFoilFlag(particles.x[ip], particles.y[ip])

            while zrl > 0:
                # If in the foil, tally the hit and start tracking
                if foil_flag == 1:
                    self.nHits += 1
                    directionfac = self.getDirection(particles, ip)
                    rl = zrl * directionfac # remaining Foil length in the direction of particle momentum (in m)

                    beta = self.getBeta(particles, ip)
                    gamma = 1/math.sqrt(1-beta**2)
                    p = self.getP(particles, ip) # in ev, to be checked!
                    theta = 0.0136 / (beta * p) / math.sqrt(radlengthfac)
                    pfac = self.getPFactor(particles, ip)
                    KE = (gamma-1)*particles.mass0*1e-9 # kinetic energy of particle in GeV
                    ecross = cs.get_elastic_crosssection(KE, self.ma_)
                    icross = cs.get_inelastic_crosssection(KE, self.ma_)
                    rcross, thx, thy = mi.ruthScattJackson(stepsize, z, a, density, beta, 0, pfac, thx, thy)
                    totcross = ecross + icross + rcross
                    meanfreepath = cs.get_a(self.ma_) / ((nAvogadro) * density*1e3 * (totcross*1.0e-28)) # in m?
                    stepsize = -meanfreepath * math.log(random.uniform(0,1)) # stepsize to be taken (in m)?

                    # Take the step but no nuclear scattering event
                    # Only MCS and energy loss
                    if stepsize > rl:
                        stepsize = rl + dlength
                        foil_flag, zrl, rl = self.takeStep(particles, ip, z, a, density, stepsize, zrl, rl, foil_flag)
                    # Take the step and allow nuclear scatter
                    if stepsize <= rl:
                        foil_flag, zrl, rl = self.takeStep(particles, ip, z, a, density, stepsize, zrl, rl, foil_flag)

                        # If it still exists after MCS and energy loss, nuclear scatter
                        if foil_flag == 1 and zrl > 0:
                            print('foil_flag still 1 and zrl>0')
                            beta = self.getBeta(particles, ip)
                            p = self.getP(particles, ip)
                            theta = 0.0136 / (beta * p) / math.sqrt(radlengthfac)
                            pfac = self.getPFactor(particles, ip)
                            ecross = cs.get_elastic_crosssection(KE, self.ma_)
                            icross = cs.get_inelastic_crosssection(KE, self.ma_)
                            rcross, thx, thy = mi.ruthScattJackson(stepsize, z, a, density, beta, 0, pfac, thx, thy)
                            totcross = ecross + icross + rcross
                            e_frac = ecross / totcross
                            i_frac = icross / totcross
                            r_frac = rcross / totcross
                            choice = random.uniform(0,1)

                            # Nuclear Elastic Scattering
                            if choice >= 0. and choice <= e_frac:
                                print('nuclear elastic scattering!!!')
                                if KE <= 0.4:
                                    t = mi.elastic_t(p, a)
                                if KE > 0.4:
                                    t = -math.log(random.uniform(0,1)) / b_pN
                                dp_x, dp_y = mi.momentumKick(t, p, dp_x, dp_y)
                                particles.px[ip] += dp_x * pfac
                                particles.py[ip] += dp_y * pfac

                            # Rutherford Coulomb scattering
                            if (choice > e_frac) and (choice <= (1 - i_frac)):
                                print('ratherford coulomb scattering')
                                rcross, thx, thy = mi.ruthScattJackson(stepsize, z, a, density, beta, 1, pfac, thx, thy)
                                xpfac = particles.px[ip] / pfac
                                ypfac = particles.py[ip] / pfac
                                anglex = math.atan(xpfac) + thx
                                angley = math.atan(ypfac) + thy                                
                                particles.px[ip] = math.tan(anglex) * pfac
                                particles.py[ip] = math.tan(angley) * pfac

                            # Nuclear Inelastic absorption
                            if ((choice > (1 - i_frac)) and (choice <= 1)):
                                print('Nucmlear inelastic absorption')
                                foil_flag, zrl = self.loseParticle(particles, ip, foil_flag, zrl)
                
                # If outside of the foil, drift particle
                else:
                    #self.driftParticle(particles, ip, length)
                    self.driftParticle(particles, ip, self.length_)
                    zrl = 0
                    print('Drifting particle.')


    def driftParticle(self, particles, ip, length):
        """
        Drifts a single particle. 

        Parameters:
            particles: xpart.Particles() object
            ip: particle index
            length = length of the drift

        Returns:
            nothing 
        """
        if length>0:
            # needs to be adapted for xsuite
            # not working (for now)
            beta = particles.beta0[ip]
            gamma = 1/math.sqrt(1-beta**2)
            gamma2i = 1.0/(gamma**2)
            dp_p_coeff = 1.0 / (particles.p0c[ip] * beta)
            dp_p = particles.delta[ip] #arr[5] * dp_p_coeff # ?
            KNL = 1.0 / (1.0 + dp_p)
            particles.x[ip] += KNL * length * particles.px[ip]
            particles.y[ip] += KNL * length * particles.py[ip]
            phifac = (particles.px[ip] ** 2 + particles.py[ip] ** 2 + dp_p ** 2 * gamma2i) / 2.0
            phifac = (phifac * KNL - dp_p * gamma2i) * KNL
            #arr[4] -= length * phifac
        else:
            pass


    def getDirection(self, particles, ip) -> float:
        """
        Gets the normalized unit vector direction of particle momentum

        Parameters:
            particles: xpart.Particles() object
            ip: particle index

        Returns:
            float 
        """
        pfac = self.getPFactor(particles, ip)
        xpfac = particles.px[ip] / pfac
        ypfac = particles.py[ip] / pfac
        return math.sqrt(1.0 + xpfac**2 + ypfac**2)


    def getPFactor(self, particles, ip):
        """
        Gets the fractional deviation from sync part momentum (1+dp/p).

        Parameters:
            particles: xpart.Particles() object
            ip: particle index

        Returns:
            float
        """
        # particles.delta := ( Pc m0/m - p0c )/p0c
        return (1.0+particles.delta[ip])
    

    def getBeta(self, particles, ip) -> float:
        """
        Gets current beta for particle. 

        Parameters:
            particles: xpart.Particles() object
            ip: particle index

        Returns:
            float
        """
        return particles.rvv[ip]*particles.beta0[ip]
    

    def getP(self, particles, ip):
        """
        Gets current momentum of particle in eV

        Parameters:
            particles: xpart.Particles() object
            ip: particle index

        Returns:
            float
        """
        # Note: valid only for m0/m = 1 (same species of particles)
        return (particles.delta[ip]*particles.p0c[ip]+particles.p0c[ip])
    

    def takeStep(self, particles, ip, z: float, a: float, density: float, stepsize: float, zrl: float, rl: float, foil_flag: int):
        """
        Takes a step in the Foil, allows MCS and ionization energy loss.
        Loses the particle if the energy is too low and updates the tracking params.

        Parameters:
            particles: xpart.Particles() object
            ip: particle index
            z: z number of material
            a: a number of material
            density: density of material
            stepsize: the stepsize to be taken
            zrl: remaining Foil length in the z direction
            rl:	remaining Foil length in the direction of particle momentum
            foil_flag: flag for in or out of the Foil.

        Returns:
            Nothing
        """
        beta = self.getBeta(particles, ip)
        p = self.getP(particles, ip)
        pfac = self.getPFactor(particles, ip)
        
        mi.mcsJackson(stepsize, z, a, density, beta, pfac, particles, ip)
        dE = mi.ionEnergyLoss(beta, z, a) # GeV/(kg/m^2)
        dE = -dE * density * stepsize * 1e9 # energy loss in eV 
        Enew = math.sqrt(p**2+particles.mass0**2)+dE
        pnew = math.sqrt(Enew**2-particles.mass0**2)
        dp = pnew-p
        
        # update particles object
        particles.ptau[ip] += dE/particles.p0c[ip] # this assumes m0/m=1
        particles.delta[ip] += dp/particles.p0c[ip]
        #particles.pzeta[ip] = particles.ptau[ip]/particles.beta0[ip] # read only!
        
        energy_threshold = 0.02*1e9 # in eV
        if (particles.ptau[ip]*particles.p0c[ip]+particles.energy0[ip]) < energy_threshold:
            foil_flag, zrl = self.loseParticle(particles, ip, foil_flag, zrl)
        else:
            # particle coordinates have slightly changed due to MCS
            directionfac = self.getDirection(particles, ip)
            zrl -= stepsize / directionfac
            rl = zrl * directionfac
            foil_flag = self.checkFoilFlag(particles.x[ip], particles.y[ip])
        
        return foil_flag, zrl, rl

        
    def loseParticle(self, particles, ip, foil_flag, zrl):
        """
        Lose the particle from the alive bunch and add it to the lost bunch

        Parameters:
            particles: xpart.Particles() object
            ip: particle index
            foil_flag: flag for if particle is still in Foil (and alive)
            zrt: remaining length in z the particle has in Foil

        Returns:
            Nothing
        """
        particles.state[ip]=-350 # lost by foil absorption
        self.nLost += 1
        foil_flag = 0
        zrl = -1.0

        return foil_flag, zrl

# %%
