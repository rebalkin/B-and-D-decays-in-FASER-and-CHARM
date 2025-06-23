import numpy as np
from . import kinematics as kin

#### FASER/2 Geometry ####
# in meters
# Radius
RFASER = 0.1
RFASER2 =  1
#Length of the decay volume
DFASER = 1.5
DFASER2  = 10
# Distance from the beam line to the front of the decay volume
LFASER = 480
LFASER2 = 650

thetaMax_FASER = np.arctan(RFASER/(LFASER+DFASER))
thetaMax_FASER2 = np.arctan(RFASER2/(LFASER2+DFASER2))

log_thetaMax_FASER = np.log10(thetaMax_FASER)
log_thetaMax_FASER2 = np.log10(thetaMax_FASER2)


# Total yields in LHC - B meson
NB_total_LHC = 6.68619*10**13
NB_total_HL_LHC = 6.68619*10**14

# Total yields in LHC - D meson
ND_total_LHC = 1.02234*10**15
ND_total_HL_LHC = 1.02234*10**16

# Decay factors FASER/2
# tau should have same units as L (ie meters)

def p_volume_decay_FASER(p, m, tau):
    return kin.p_volume_decay(p, m, tau, LFASER, DFASER)

def p_volume_decay_FASER2(p, m, tau):
    return kin.p_volume_decay(p, m, tau, LFASER2, DFASER2)


def get_events_FASER(Xspec, mX, tau,decaying_meson):
    # Xspec has the form (log10theta, log10p, weight)
    logtheta = Xspec[:,0]
    p = 10**Xspec[:,1]
    weights = Xspec[:,2]
    
    mask_FASER  = logtheta < log_thetaMax_FASER
    mask_FASER2 = logtheta < log_thetaMax_FASER2
    
    prob_FASER  = p_volume_decay_FASER(p, mX, tau)
    prob_FASER2 = p_volume_decay_FASER2(p, mX, tau)
    if decaying_meson == 'B':
        N_total_FASER = NB_total_LHC
        N_total_FASER2 = NB_total_HL_LHC
    elif decaying_meson == 'D':
        N_total_FASER = ND_total_LHC
        N_total_FASER2 = ND_total_HL_LHC
    else:
        raise ValueError("Decaying meson must be either 'B' or 'D'")
    NFASER = N_total_FASER*np.dot(weights[mask_FASER], prob_FASER[mask_FASER])
    NFASER2 = N_total_FASER2*np.dot(weights[mask_FASER2], prob_FASER2[mask_FASER2])
    
   
    return np.array([mX, tau, NFASER, NFASER2])

import numpy as np
from . import kinematics as kin



########## CHARM GEOMETRY #############

# Everything is in meters
#CHARM decay volume is modeled as 3m x 3m x 35m  rectangle 

RCHARM = 3 

DCHARM = 35

LCHARM = 480

CHARM_Offset = 5 # There was a 5m offset between the center of the decay volume and the beam line

CHARM_delta = 3.5 # Offset between the bottom of the decay volume and the beam line

def phiMin_CHARM(theta):
    return np.arctan2(np.tan(theta)*(LCHARM+DCHARM),+RCHARM/2)

def phiMax_CHARM(theta):
    return np.arctan2(np.tan(theta)*(LCHARM+DCHARM),-RCHARM/2)


# Total B and D mesons produced in CHARM
NB_total_CHARM = 7.9*10**11
ND_total_CHARM = 1.07*10**15


# Decay factors
# tau should have same units as L (ie meters)

def p_volume_decay_box(p,theta, phi, L, D, R, delta, m, tau):
    
    dMax = np.abs((L+D)/np.cos(theta)) #if it's negative, it won't pass the geometrical acceptance 
    
    dMin1 = np.abs(L/np.cos(theta))
    
    dMin2 = np.abs(delta/(np.sin(theta)*np.sin(phi)+10**(-10))) #added a small number to avoid deviding by 0
    
    dMin = np.maximum(dMin1,dMin2);
 
    return kin.p_decay(p, m, tau, dMin)-kin.p_decay(p, m, tau, dMax)


#Angular acceptance which requires that the decay particle trajectory passes through the back panel of the decay volume
def angular_acc(theta, phi, L, D, R, delta):
    
    x = np.tan(theta)*np.cos(phi) # in units of L+D
    y = np.tan(theta)*np.sin(phi) # in units of L+D
    
    maxX = 0.5 * R / (L+D)
    
    minY = (delta)/(L+D)
    maxY = (delta+R)/(L+D)
    
    inXrange = (-maxX<x) & (x<maxX) 
    inYrange = (minY<y) & (y<maxY)
    isForward = (theta <= np.pi/2)
    
    return inXrange & inYrange & isForward

# Function to get the number of events in a box-shaped detector with given the parameters:
# L is the length of the decay volume, D is the depth, R is the radius, delta is the offset from the beam line

def get_events_box(Xspec, mX, tau,L, D, R, delta):
    # spec has the form (theteX_lab,phiX_lab,pX_lab,weight)
    theta = Xspec[:,0]
    phi = Xspec[:,1]
    p = Xspec[:,2]
    weights = Xspec[:,3]
   
    acc_mask = angular_acc(theta, phi, L, D, R, delta)
    decay_prob  = p_volume_decay_box(p,theta, phi, L, D, R, delta, mX, tau)

    N_events = (np.dot(weights[acc_mask], decay_prob[acc_mask]))
    #Still need to normalize the number of events by the total number of produced mesons
    return np.array([mX, tau, N_events])



def get_events_CHARM(Xspec, mX, tau,decaying_meson):
    # Now with the parameters for CHARM
    res =  get_events_box(Xspec, mX, tau,LCHARM, DCHARM, RCHARM, CHARM_delta)
    N_unnormalized = res[2]
    if decaying_meson == 'B':
        N_total_CHARM = NB_total_CHARM
    elif decaying_meson == 'D':
        N_total_CHARM = ND_total_CHARM
    else:
        raise ValueError("Decaying meson must be either 'B' or 'D'")
    return np.array([mX, tau,N_total_CHARM*N_unnormalized])

def get_events(Xspec, mX, tau,decaying_meson,exp):
    if exp == 'FASER':
        return get_events_FASER(Xspec, mX, tau,decaying_meson)
    elif exp == 'CHARM':
        return get_events_CHARM(Xspec,  mX, tau,decaying_meson)
    else:
        raise ValueError("Experiment must be either 'FASER' or 'CHARM'")


         