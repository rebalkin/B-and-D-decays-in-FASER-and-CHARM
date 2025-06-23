### HEADER
import numpy as np
from . import kinematics as kin
### Constants ### 

mB0, mK0, mD0, mPi0 = np.array([5279.65, 497.611, 1864.84, 134.976])* 10**(-3) # Meson masses in GeV;

#import spectra
# Load the spectra data
b_spectrum_LHC = np.genfromtxt('spectra/bbbar_spectrum_normalized_LHC.txt', delimiter=' ', skip_header=0)
c_spectrum_LHC = np.genfromtxt('spectra/ccbar_spectrum_normalized_LHC.txt', delimiter=' ', skip_header=0)
b_spectrum_CHARM = np.genfromtxt('spectra/bbbar_spectrum_normalized_CHARM.txt', delimiter=' ', skip_header=0)
c_spectrum_CHARM = np.genfromtxt('spectra/ccbar_spectrum_normalized_CHARM.txt', delimiter=' ', skip_header=0)

def decay_spectrum_1_bin(m1,m2,m3,theta1_lab,p1_lab,binWeight,Ntheta,Nphi,eps,exp):
    ## values for 4pi solid angle spread for decaying particle in the rest frame of (1)
    costheta_R = np.linspace(-1+eps,1-eps,Ntheta) #probing the edge values usually generates undesired artifacts (particle is produced colinear with parent), we avoid it by defining a threshold eps 
    phi_R = np.linspace(0,2*np.pi,Nphi)
    ## spread the total bin weight across that 4pi solid angle
    weight = binWeight/(Ntheta*Nphi)
    if exp == 'LHC':
        result = np.empty([int(Ntheta),int(Nphi),3],dtype=float)
    elif exp == 'CHARM':
        result = np.empty([int(Ntheta),int(Nphi),4],dtype=float)
    else:
        raise ValueError("Experiment must be either 'LHC' or 'CHARM'")
    ## creating the spectrum of the decay particle {theta_x,p_x,weight} in the lab frame
    for itheta, costheta in enumerate(costheta_R):
        for iphi, phi in enumerate(phi_R):
            theta = np.arccos(costheta)
            p3Lab = kin.p3_lab(m1,m2,m3,theta1_lab,p1_lab,theta,phi)
            if exp == 'LHC':
                result[itheta,iphi,0] = np.log10(kin.getTheta(p3Lab))
                result[itheta,iphi,1] = np.log10(kin.getMom(p3Lab))
                result[itheta,iphi,2] = weight
            elif exp == 'CHARM':
                result[itheta,iphi,0] = kin.getTheta(p3Lab)
                result[itheta,iphi,1] = kin.getPhi(p3Lab)
                result[itheta,iphi,2] = kin.getMom(p3Lab)
                result[itheta,iphi,3] = weight
    if exp == 'LHC':
        return result.reshape(-1,3)
    elif exp == 'CHARM':
        return result.reshape(-1,4)
    

    #Given the input spectrum (spec), generated the full production spectrum of 3 from 1->2+3 decay.
    #spec is a numpy array of shape (N, 3) where each row is [theta1_lab, p1_lab, binWeight]
    #m1, m2, m3 are the masses of the particles in GeV.
    #Ntheta and Nphi are the number of bins in theta and phi respectively.
    #eps is a small value to avoid numerical issues at the edges of the theta range.
    #Returns a numpy array of shape (N * Ntheta * Nphi, 3) where each row is [log10(theta), log10(p), weight].
def spectrum_all_bins(spec,m1,m2,m3,Ntheta,Nphi,eps,exp):

    spec_length = spec.shape[0]
    spectrum = []
    
    for i in range(spec_length):
        
        theta1_lab, p1_lab, binWeight = spec[i]
        spectrum.append(decay_spectrum_1_bin(m1,m2,m3,theta1_lab,p1_lab,binWeight,Ntheta,Nphi,eps,exp))
    if exp == 'LHC':
        return np.array(spectrum).reshape(-1, 3)
    elif exp == 'CHARM':
        return np.array(spectrum).reshape(-1, 4)
    

def Xspectrum(mX,decaying_meson,exp,Ntheta=100,Nphi=50,eps=10**(-6)):
    if decaying_meson == 'B':
        m1 = mB0
        m2 = mK0
    elif decaying_meson == 'D':
        m1 = mD0
        m2 = mPi0
    else:
        raise ValueError("decaying_meson must be either 'B' or 'D'")
    
    if exp == 'LHC' and decaying_meson == 'B':
        spec = b_spectrum_LHC
    elif exp == 'LHC' and decaying_meson == 'D':
        spec = c_spectrum_LHC
    elif exp == 'CHARM' and decaying_meson == 'B':
        spec = b_spectrum_CHARM
    elif exp == 'CHARM' and decaying_meson == 'D':
        spec = c_spectrum_CHARM
    else:
        raise ValueError("Invalid combination of experiment and decaying meson")
    return spectrum_all_bins(spec,m1,m2,mX,Ntheta,Nphi,eps,exp)
