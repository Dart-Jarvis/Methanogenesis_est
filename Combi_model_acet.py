#--------------------------------------------------------------------
# Modified from Ed's model described in Taenzer et al., 2020, GCA.
# Created by Jiawen Li in November, 2023.
#--------------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from matplotlib.ticker import MultipleLocator
import statistics

# Font dictionary
font = {'family': 'sans serif',
        'color':  'black',
        'weight': 'normal',
        'size': 36,
        }
font_labels = {'family': 'sans serif',
               'color':  'black',
               'weight': 'normal',
               'size': 36,
                }
font_ticks = {'family': 'sans serif',
               'color':  'black',
               'weight': 'normal',
               'size': 32,
                }

# define functions for calculating isotopologue abundances from delta values
def CH3_abundance(dD, d13C, D13CD, DD2):
    RD=0.00015576*(dD/1000+1)
    RC=0.01118*(d13C/1000+1)
    xH=1/(1+RD)
    xD=1-xH
    x12C=1/(1+RC)
    x13C=1-x12C
    # Here start calculating the abundances of isotopologues
    CH3=x12C*xH**3 # Stochastic. Add a multiplier for the clumping of each species (default 1), same below
    QH3=x13C*xH**3
    CH2D=3*x12C*(xH**2)*xD
    QH2D=(D13CD/1000+1)*3*x13C*(xH**2)*xD # Clumping value measured by Daniel Stolper
    CHD2=(DD2/1000+1)*3*x12C*xH*xD**2
    return CH3,QH3,CH2D,QH2D,CHD2

def H2O_abundance(dD):
    RD=0.00015576*(dD/1000+1)
    xH=1/(1+RD)
    xD=1-xH
    # The abundances of isotopologues
    H2O=xH**2
    HDO=2*xH*xD
    D2O=xD**2
    return H2O,HDO,D2O

def CH4_abundance(RD, RC):
    xH=1/(1+RD)
    xD=1-xH
    x12C=1/(1+RC)
    x13C=1-x12C
    # Here start calculating the abundances of isotopologues
    CH4=x12C*xH**4 # Stochastic.
    QH4=x13C*xH**4
    CH3D=4*x12C*(xH**3)*xD
    QH3D=4*x13C*(xH**3)*xD
    CH2D2=6*x12C*(xH**2)*xD**2
    QH2D2=6*x13C*(xH**2)*xD**2
    CHD3=4*x12C*(xD**3)*xH
    QHD3=4*x13C*(xD**3)*xH
    CD4=x12C*xD**4
    QD4=x13C*xD**4
    return CH4,QH4,CH3D,QH3D,CH2D2,QH2D2,CHD3,QHD3,CD4,QD4

# Define a function to calculate collison frequency and 
# resulting CH4 isotopologue abundances (this is to make MCMC more convinient)
def col_freq(aDp, aDs, dD_MeOH, d13C_MeOH, dD_H2O, a13, g13CDp, gDDp, g13CDs, gDDs, D13CD, DD2):
# Calculate the isotopologue abundances for CH3 and H2O
    xCH3,xQH3,xCH2D,xQH2D,xCHD2=CH3_abundance(dD_MeOH,d13C_MeOH, D13CD, DD2)
    xH2O,xHDO,xD2O=H2O_abundance(dD_H2O)
    # Calculate collision frequency
    # A total of 10*2 equations for the isotopologues of interest (including both 12C and 13C versions)
    rxn1=xCH3*xH2O # H2O+12CH3 -> 12CH4
    rxn2=aDs*xCH2D*xH2O # H2O + 12CH2D -> 12CH3D 
    rxn3=(aDs**2)*gDDs*xCHD2*xH2O # H2O + 12CHD2 -> 12CH2D2
    rxn4=1/2*xHDO*xCH3 # HDO + 12CH3 -> 12CH4
    rxn5=1/2*aDp*xCH3*xHDO # HDO + 12CH3 -> 12CH3D
    rxn6=1/2*aDs*xCH2D*xHDO # HDO + 12CH2D -> 12CH3D
    rxn7=1/2*aDs*aDp*gDDp*xCH2D*xHDO # HDO + 12CH2D -> 12CH2D2
    rxn8=1/2*(aDs**2)*gDDs*xHDO*xCHD2 # HDO + 12CHD2 -> 12CH2D2 
    # HDO + 12CHD2 -> 12CHD3 
    # not included since 12CHD3 is not the isotopologue of interest, and much lower in abundance
    rxn9=aDp*xD2O*xCH3 # D2O + 12CH3 -> 12CH3D
    rxn10=aDp*aDs*gDDp*xCH2D*xD2O # D2O + 12CH2D -> 12CH2D2
    # For 13C isotopologues
    rxn11=a13*xQH3*xH2O # H2O+13CH3 -> 13CH4
    rxn12=aDs*a13*g13CDs*xQH2D*xH2O # H2O + 13CH2D -> 13CH3D 
    # rxn13=(aDs**2)*gDD*xCHD2*xH2O # H2O + 13CHD2 -> 13CH2D2
    rxn13=1/2*a13*xHDO*xQH3 # HDO + 13CH3 -> 13CH4
    rxn14=1/2*a13*aDp*g13CDp*xQH3*xHDO # HDO + 13CH3 -> 13CH3D
    rxn15=1/2*a13*aDs*g13CDs*xQH2D*xHDO # HDO + 13CH2D -> 13CH3D
    # rxn17=1/2*aDs*aDp*gDD*xCH2D*xHDO # HDO + 13CH2D -> 13CH2D2
    # rxn18=1/2*(aDs**2)*gDD*xHDO*xCHD2 # HDO + 13CHD2 -> 13CH2D2 
    # HDO + 13CHD2 -> 13CHD3 
    # not included since 12CHD3 is not the isotopologue of interest, and much lower in abundance
    rxn16=a13*aDp*g13CDp*xD2O*xQH3 # D2O + 13CH3 -> 13CH3D
    # rxn20=aDp*aDs*gDD*xCH2D*xD2O # D2O + 13CH2D -> 13CH2D2
    # Abundance of each isotopologue
    tCH4=rxn1+rxn4
    tCH3D=rxn2+rxn5+rxn6+rxn9
    tCH2D2=rxn3+rxn7+rxn8+rxn10
    tQH4=rxn11+rxn13
    tQH3D=rxn12+rxn14+rxn15+rxn16
    RD_CH4=(tCH3D+2*tCH2D2+tQH3D)/(4*tCH4+3*tCH3D+2*tCH2D2+4*tQH4+3*tQH3D) # total D/H ratio in CH4
    RC_CH4= (tQH4+tQH3D)/(tCH4+tCH3D+tCH2D2) # total 13C/12C ratio in CH4
    return RD_CH4, RC_CH4, tCH4, tCH3D, tCH2D2, tQH4, tQH3D

# Define a function to calculate log likelihood
def normpdf_ll(mu,sigma,x):
    log_likelihood=-(x-mu)*(x-mu)/(2*sigma*sigma)
    return log_likelihood

# Define a function to perform Monte Carlo Error Propagation
# Propagate the errors in the measurements to the final clumped isotope values
def MCEP(param, eparam, nsims): 
    # Unpack the needed values
    slope, intercept, dD_MeOH, d13C_MeOH, dD_H2O, f, a13, g13CDp, gDDp, g13CDs, gDDs, D13CD, DD2 = param
    eslope, eintercept, edD_MeOH, ed13C_MeOH, edD_H2O, ef, ea13, eg13CDp, egDDp, eg13CDs, egDDs, eD13CD, eDD2 = eparam
    a=1 # Controls how many sigma's
    isotope=np.zeros((nsims,4))
    alpha=np.zeros((nsims,2))
    gamma=np.zeros((nsims,4))
    for i in range(nsims):
        # Randomly picking values, assuming all variables are normally distributed
        slope_prop=np.random.normal(slope,a*eslope)
        intercept_prop=np.random.normal(intercept,a*eintercept)
        dD_MeOH_prop=np.random.normal(dD_MeOH,a*edD_MeOH)
        d13C_MeOH_prop=np.random.normal(d13C_MeOH,a*ed13C_MeOH)
        dD_H2O_prop=np.random.normal(dD_H2O,a*edD_H2O)
        f_prop=np.random.normal(f,a*ef)
        a13_prop=np.random.normal(a13,a*ea13)
        aDp_prop=slope_prop/f_prop
        aDs_prop=intercept_prop/((1-f_prop)*(dD_MeOH_prop+1000))
        g13CDp_prop=np.random.normal(g13CDp,a*eg13CDp)
        gDDp_prop=np.random.normal(gDDp,a*egDDp)
        g13CDs_prop=np.random.normal(g13CDs,a*eg13CDs)
        gDDs_prop=np.random.normal(gDDs,a*egDDs)
        D13CD_prop=np.random.normal(D13CD,a*eD13CD)
        DD2_prop=np.random.normal(DD2,a*eDD2)
        RD_CH4, RC_CH4, tCH4, tCH3D, tCH2D2, tQH4, tQH3D=col_freq(aDp_prop, aDs_prop, dD_MeOH_prop, d13C_MeOH_prop, dD_H2O_prop, a13_prop, 
                                                                  g13CDp_prop, gDDp_prop, g13CDs_prop, gDDs_prop, D13CD_prop,DD2_prop)
        stCH4,stQH4,stCH3D,stQH3D,stCH2D2,stQH2D2,stCHD3,stQHD3,stCD4,stQD4=CH4_abundance(RD_CH4,RC_CH4)
        # Convert to isotope values
        dD_CH4_prop=1000*(RD_CH4/0.00015576-1)
        d13C_CH4_prop=1000*(RC_CH4/0.01118-1)
        D13CH3D_CH4_prop=1000*((tQH3D/tCH4)/(stQH3D/stCH4)-1)
        D12CH2D2_CH4_prop=1000*((tCH2D2/tCH4)/(stCH2D2/stCH4)-1)
        # Store isotope values
        isotope[i,:]=[d13C_CH4_prop,dD_CH4_prop,D13CH3D_CH4_prop,D12CH2D2_CH4_prop]
        alpha[i,:]=[aDp_prop,aDs_prop]
        gamma[i,:]=[g13CDp_prop,gDDp_prop,g13CDs_prop,gDDs_prop]
        # Calculate the net gamma and alpha 13CD, DD
    return isotope, alpha

# Define a function to perform MCMC
# param: Input parameters (see below); eparam: errors of input parameters; 
# isotope: isotope values of the samples; eisotope: errors of isotope values; nsims: number of simulations.
# All inputs are lists except for nsims (integer)
def MCMC(param, eparam, isotope, eisotope, nsims): 
    # Unpack the needed values
    slope, intercept, dD_MeOH, d13C_MeOH, dD_H2O, f_0, a13, g13CDp_0, gDDp_0, g13CDs_0,gDDs_0, D13CD, DD2 = param
    eslope, eintercept, edD_MeOH, ed13C_MeOH, edD_H2O, ef, ea13, eg13CDp, egDDp, eg13CDs,egDDs, eD13CD, eDD2 = eparam
    d13C_CH4_0,dD_CH4_0,D13CH3D_CH4_0,D12CH2D2_CH4_0 = isotope
    ed13C_CH4,edD_CH4,eD13CH3D_CH4,eD12CH2D2_CH4 = eisotope
    k=1 # Calculating likelihood at 1 sigma
    a=1 # Random walk step
    # Randomly choose initial values for parameters, assuming the variables are normally distributed
    # 3 variables to randomize: f(and corresponding aDp, aDs), g13CD, gDD
    f=np.random.normal(f_0,a*ef)
    aDp=slope/f
    aDs=intercept/((1-f)*(dD_MeOH+1000))
    g13CDp=np.random.normal(g13CDp_0,a*eg13CDp)
    gDDp=np.random.normal(gDDp_0,a*egDDp)
    g13CDs=np.random.normal(g13CDs_0,a*eg13CDs)
    gDDs=np.random.normal(gDDs_0,a*egDDs)
    # Calculate collision frequency
    RD_CH4, RC_CH4, tCH4, tCH3D, tCH2D2, tQH4, tQH3D=col_freq(aDp, aDs, dD_MeOH, d13C_MeOH, dD_H2O, a13, g13CDp, gDDp, g13CDs, gDDs, D13CD, DD2)
    stCH4,stQH4,stCH3D,stQH3D,stCH2D2,stQH2D2,stCHD3,stQHD3,stCD4,stQD4=CH4_abundance(RD_CH4,RC_CH4)
    # Convert to isotope values
    dD_CH4=1000*(RD_CH4/0.00015576-1)
    d13C_CH4=1000*(RC_CH4/0.01118-1)
    D13CH3D_CH4=1000*((tQH3D/tCH4)/(stQH3D/stCH4)-1)
    D12CH2D2_CH4=1000*((tCH2D2/tCH4)/(stCH2D2/stCH4)-1)
    # Calculate the log-likelihood
    ll=normpdf_ll(f_0,k*ef,f)+normpdf_ll(g13CDp_0,k*eg13CDp,g13CDp)+normpdf_ll(gDDp_0,k*egDDp,gDDp)
    +normpdf_ll(g13CDs_0,k*eg13CDs,g13CDs)+normpdf_ll(gDDs_0,k*egDDs,gDDs)
    +normpdf_ll(d13C_CH4_0,k*ed13C_CH4,d13C_CH4)+normpdf_ll(dD_CH4_0,k*edD_CH4,dD_CH4)
    +normpdf_ll(D13CH3D_CH4_0,k*eD13CH3D_CH4,D13CH3D_CH4)+normpdf_ll(D12CH2D2_CH4_0,k*eD12CH2D2_CH4,D12CH2D2_CH4)
    # Initialize an array to store the results of each iteration
    param_result=np.zeros((nsims+1,5)) # Result of parameters
    alpha_result=np.zeros((nsims+1,2))
    isotope_result=np.zeros((nsims+1,4))
    # Store the initial values
    param_result[0,:]=[f,g13CDp,gDDp,g13CDs,gDDs]
    alpha_result[0,:]=[aDp,aDs]
    isotope_result[0,:]=[d13C_CH4,dD_CH4,D13CH3D_CH4,D12CH2D2_CH4]
    c=0 # Counter for the acceptance of iterations
    for i in range(nsims):
        # Random walk
        f_prop=np.random.normal(f,a*ef)
        aDp_prop=slope/f_prop
        aDs_prop=intercept/((1-f_prop)*(dD_MeOH+1000))
        g13CDp_prop=np.random.normal(g13CDp,a*eg13CDp)
        gDDp_prop=np.random.normal(gDDp,a*egDDp)
        g13CDs_prop=np.random.normal(g13CDs,a*eg13CDs)
        gDDs_prop=np.random.normal(gDDs,a*egDDs)
        RD_CH4, RC_CH4, tCH4, tCH3D, tCH2D2, tQH4, tQH3D=col_freq(aDp_prop, aDs_prop, dD_MeOH, d13C_MeOH, dD_H2O, a13, 
                                                                  g13CDp_prop, gDDp_prop, g13CDs_prop, gDDs_prop, D13CD, DD2)
        stCH4,stQH4,stCH3D,stQH3D,stCH2D2,stQH2D2,stCHD3,stQHD3,stCD4,stQD4=CH4_abundance(RD_CH4,RC_CH4)
        # Convert to isotope values
        dD_CH4_prop=1000*(RD_CH4/0.00015576-1)
        d13C_CH4_prop=1000*(RC_CH4/0.01118-1)
        D13CH3D_CH4_prop=1000*((tQH3D/tCH4)/(stQH3D/stCH4)-1)
        D12CH2D2_CH4_prop=1000*((tCH2D2/tCH4)/(stCH2D2/stCH4)-1)
        # Calculate the log-likelihood
        ll_prop=normpdf_ll(f_0,k*ef,f_prop)+normpdf_ll(g13CDp_0,k*eg13CDp,g13CDp_prop)+normpdf_ll(gDDp_0,k*egDDp,gDDp_prop)
        +normpdf_ll(g13CDs_0,k*eg13CDs,g13CDs_prop)+normpdf_ll(gDDs_0,k*egDDs,gDDs_prop)
        +normpdf_ll(d13C_CH4_0,k*ed13C_CH4,d13C_CH4_prop)+normpdf_ll(dD_CH4_0,k*edD_CH4,dD_CH4_prop)
        +normpdf_ll(D13CH3D_CH4_0,k*eD13CH3D_CH4,D13CH3D_CH4_prop)+normpdf_ll(D12CH2D2_CH4_0,k*eD12CH2D2_CH4,D12CH2D2_CH4_prop)
        # Accept Condition
        if np.log(np.random.uniform(0,1))<(ll_prop-ll):
            ll=ll_prop
            f=f_prop
            aDp=aDp_prop
            aDs=aDs_prop
            g13CDp=g13CDp_prop
            gDDp=gDDp_prop
            g13CDs=g13CDs_prop
            gDDs=gDDs_prop
            dD_CH4=dD_CH4_prop
            d13C_CH4=d13C_CH4_prop
            D13CH3D_CH4=D13CH3D_CH4_prop
            D12CH2D2_CH4=D12CH2D2_CH4_prop
            c+=1
        # Record the values of this iteration
        # Even if the new values are not accepted, this still counts as one iteration
        param_result[i+1,:]=[f, g13CDp, gDDp, g13CDs, gDDs]
        alpha_result[i+1,:]=[aDp,aDs]
        isotope_result[i+1,:]=[d13C_CH4,dD_CH4,D13CH3D_CH4,D12CH2D2_CH4]
    return param_result,alpha_result,isotope_result,c

# Load the equilibrium clumped isotope values
equib = pd.read_csv('equilibrium.csv')

# -------------------------- MCEP -----------------------------------------
# --------------------- Input parameters (initial) --------------------------------
# Use these values in the MCMC
# Slope and intercept of the regression lines
def Plot_MCEP_acet(dD_H2O):
    slope=0.087909 # 0.087909 for 0.92 ;0.0913157 for 0.925 ;0.0922087 for 0.93
    intercept=616.925 #  616.925 for 0.92; 610.325 for 0.93; 612.676 for 0.925
    # Measured isotope values
    dD_MeOH=-159.7 # +/- 0.246
    d13C_MeOH=-32.2 # +/- 0.017
    D13CD=0.0
    DD2=0.0
    # Calculate the fractionation factors here
    f=0.25 # +/- ef pick your favourite
    a13= 0.9817 # 0.98 for no mixing; 0.9817 for r=0.92; alpha13C, averaged across all spikes
    aDp= slope/f # Primary fractionation factor alphaD, presumably 0.84 for mcr from Scheller et al., 2013
    aDs= intercept/((1-f)*(dD_MeOH+1000)) # Secondary fractionation factor alphaD
    g13CDp=0.996 # 0.975 # Primary clumping factor gamma 13CD +/- e
    gDDp=0.999 # 0.997 # Primary clumping factor gamma DD +/- e
    g13CDs=0.9975 # 0.996 # Secondary clumping factor gamma 13CD +/- e
    gDDs=0.999 # 0.997 # Secondary clumping factor gamma DD +/- e
    # Errors
    eslope=7.716e-05 # 7.716e-05 for 0.92; 8.057e-5 for 0.93; 6.571e-05 for 0.925
    eintercept=0.0974266 # 0.0974266 for 0.92; 0.1007 for 0.93; 0.0955374 for 0.925
    edD_MeOH=1.9 
    ed13C_MeOH=0.05 
    edD_H2O=0
    eD13CD=0.3
    eDD2=3.0
    ef=0.01 # Pick your favourites
    ea13=0.0035 # 0.0035 for r=0.92; 0.003 for no mixing; Stdev of all alpha13C
    eg13CDp=0.002 # Pick your favourite
    egDDp=0.002 # Pick your favourite
    eg13CDs=0.002 # Pick your favourite
    egDDs=0.002 # Pick your favourite
    n=1000 # number of simulations
    o=np.zeros([len(dD_H2O),4])
    oe=np.zeros([len(dD_H2O),4])
    ao=np.array([])
    # Pack up the parameters and errors
    for i in range (len(dD_H2O)):
        param_dict=[slope, intercept, dD_MeOH, d13C_MeOH, dD_H2O[i], f, a13, g13CDp, gDDp, g13CDs, gDDs, D13CD, DD2]
        eparam_dict=[eslope, eintercept, edD_MeOH, ed13C_MeOH, edD_H2O, ef, ea13, eg13CDp, egDDp, eg13CDs, egDDs, eD13CD,eDD2]
        iso_result,alpha_result=MCEP(param_dict, eparam_dict,n) # A matrix with each variable as a column
        o[i,:]=np.mean(iso_result,axis=0)
        oe[i,:]=np.std(iso_result,axis=0)
        ao=np.append(ao,alpha_result)
    return o,oe,ao

dD_H2O=list(range(-150,9050,50))
# Target values for MCEP
# Acetate
# Lab water
d13CH3D_target_50=[-285.448, -286.947]
d12CH2D2_target_50=[-495.783, -497.874]
# +3000
d13CH3D_target_3000=[ 45.083 ,40.063]
d12CH2D2_target_3000=[192.227 , 176.13]
# +8000
d13CH3D_target_8000=[ 726.605 ,680.184]
d12CH2D2_target_8000=[2444.197 ,2212.134]
# H2/CO2
# Lab water
d13CH3D_H2CO2_50=[-527.245, -526.959]
d12CH2D2_H2CO2_50=[-759.778, -759.909]
# +3000
d13CH3D_H2CO2_3000=[833.068, 821.151]
d12CH2D2_H2CO2_3000=[2627.761, 2576.11]
# +8000
d13CH3D_H2CO2_8000=[3405.96,3395.159]
d12CH2D2_H2CO2_8000=[19918.83,19952.965]

data = pd.read_csv('isotope_data.csv')
H2CO2=data[data['Label']=='4L'] # Low temp H2/CO2 (35-39)
Acet=data[data['Label']=='5L'] # Low temp acetotrophic (low temp only)
H2CO2_3000=data[data['Label']=='4LS1'] # Low temp H2/CO2, spiked +3000
H2CO2_8000=data[data['Label']=='4LS2'] # Low temp H2/CO2, spiked +8000
Acet_3000=data[data['Label']=='5LS1'] # Low temp H2/CO2, spiked +3000
Acet_8000=data[data['Label']=='5LS2'] # Low temp H2/CO2, spiked +8000
MeOH_LT=data[data['Label']=='1L'] # Low temp MeOH
MeOH=MeOH_LT.iloc[2:5]
MeOH_3000=data[data['Label']=='1LS1'] # Low temp MeOH, spiked +3000
MeOH_8000=data[data['Label']=='1LS2'] # Low temp MeOH, spiked +8000
# define a function to calculate mixing endmember of pure acetoclastic methanogenesis
# f is the input from acetoclastic methanogenesis
# 1 is the H2/CO2 endmember; 2 is the measured acetate data; 3 is the acetoclastic endmember
def calc_acet(d13C1,d13C2,dD1,dD2,d13CD1,d13CD2,dDD1,dDD2,fmix):
    # Convert delta to x
    x13C1=0.0112372*(d13C1/1000+1)
    x13C1=x13C1/(1+x13C1)
    x13C2=0.0112372*(d13C2/1000+1)
    x13C2=x13C2/(1+x13C2)
    xD1=0.00015576*(dD1/1000+1)
    xD1=xD1/(1+xD1)
    xD2=0.00015576*(dD2/1000+1)
    xD2=xD2/(1+xD2)
    r13CD1=6.34600e-6*(d13CD1/1000+1)
    r13CD2=6.34600e-6*(d13CD2/1000+1)
    rDD1=1.28866e-7*(dDD1/1000+1)
    rDD2=1.28866e-7*(dDD2/1000+1)
    # Calculate mixing, use F
    d13C3=(d13C2-d13C1*(1-fmix))/fmix
    dD3=(dD2-dD1*(1-fmix))/fmix
    x13C3=(x13C2-x13C1*(1-fmix))/fmix
    xD3=(xD2-xD1*(1-fmix))/fmix
    r13CD3=(r13CD2-r13CD1*(1-fmix))/fmix
    rDD3=(rDD2-rDD1*(1-fmix))/fmix
    # Calculate stochastic state from new x values
    sto_13CD=4*x13C3*xD3*(1-xD3)**3
    sto_DD=6*(1-x13C3)*(xD3**2)*((1-xD3)**2)
    sto_CH4=(1-x13C3)*(1-xD3)**4
    r13CD_sto=sto_13CD/sto_CH4
    rDD_sto=sto_DD/sto_CH4
    # Convert d13CD and dDD to big Delta
    D13CD3=(r13CD3/r13CD_sto-1)*1000
    DDD3=(rDD3/rDD_sto-1)*1000
    return np.array([d13C3,dD3,D13CD3,DDD3])


# Model the combinatorial effect
dD_H2O_m=np.array(range(-150,9050,50))
p1,ep1,ap1=Plot_MCEP_acet(dD_H2O_m)
aDp=ap1[::2]
aDs=ap1[1::2]
eaDp=np.std(aDp)
eaDs=np.std(aDs)
aDp=np.average(aDp)
aDs=np.average(aDs)
# ----------------------------- Plotting -------------------------------------
def quick_plot2(ax,y):
    ax.errorbar(Acet['dD_H2O'].iloc[0], iso50[y],xerr=Acet['dse_H2O'].iloc[0],yerr=0,
                markersize=16,label='Acetoclastic (' +r'$r$ ='+str(mixing_ratio[i])+')', fmt='D', 
                markerfacecolor=a, markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_3000['dD_H2O'].iloc[0],iso3000[y],xerr=Acet_3000['dse_H2O'].iloc[0],yerr=0,
                    markersize=16,fmt='D', 
                    markerfacecolor=a, markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_8000['dD_H2O'].iloc[0],iso8000[y],xerr=Acet_8000['dse_H2O'].iloc[0],yerr=0,
                    markersize=16,fmt='D', 
                    markerfacecolor=a, markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.set_xlim([-450,9000])
    ax.xaxis.set_minor_locator(MultipleLocator(400))   
    ax.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
    ax.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
    ax.set_xlabel(r'$\delta$D$_{\rm H_2O}$' ' (\u2030)',fontdict=font_labels)

def quick_plot_star(ax,y):
    ax.errorbar(Acet['dD_H2O'].iloc[0], iso50[y],xerr=Acet['dse_H2O'].iloc[0],yerr=0.0,
                markersize=40,label='Acetoclastic (' +r'$r$ ='+str(mixing_ratio[i])+')', fmt='*', 
                markerfacecolor=a, markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=10)
    ax.errorbar(Acet_3000['dD_H2O'].iloc[0],iso3000[y],xerr=Acet_3000['dse_H2O'].iloc[0],yerr=0,
                    markersize=40,fmt='*', 
                    markerfacecolor=a, markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5,  zorder=10)
    ax.errorbar(Acet_8000['dD_H2O'].iloc[0],iso8000[y],xerr=Acet_8000['dse_H2O'].iloc[0],yerr=0,
                    markersize=40,fmt='*', 
                    markerfacecolor=a, markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=10)
    ax.set_xlim([-450,9000])
    ax.xaxis.set_minor_locator(MultipleLocator(400))   
    ax.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
    ax.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
    ax.set_xlabel(r'$\delta$D$_{\rm H_2O}$' ' (\u2030)',fontdict=font_labels)

def plot_original(ax,y):
    ax.errorbar(Acet['dD_H2O'],Acet[y], xerr=0.0, yerr=0.0,
                markersize=16,label='Acetoclastic (original)', fmt='D', 
                markerfacecolor='0.0', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_3000['dD_H2O'],Acet_3000[y], xerr=0.0, yerr=0.0,
                markersize=16, fmt='D', 
                markerfacecolor='0.0', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_8000['dD_H2O'],Acet_8000[y], xerr=0.0, yerr=0.0,
                markersize=16, fmt='D', 
                markerfacecolor='0.0', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)

def quick_plot_reg(ax):
    ax.errorbar(Acet['dD_H2O'].iloc[0]+1000, iso50[1]+1000,xerr=Acet['dse_H2O'].iloc[0],yerr=0,
                markersize=20,label='Acetoclastic (' +r'$r$ ='+str(mixing_ratio[i])+')', fmt='s', 
                markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_3000['dD_H2O'].iloc[0]+1000,iso3000[1]+1000,xerr=Acet_3000['dse_H2O'].iloc[0],yerr=0,
                    markersize=20,fmt='s', 
                    markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_8000['dD_H2O'].iloc[0]+1000,iso8000[1]+1000,xerr=Acet_8000['dse_H2O'].iloc[0],yerr=0,
                    markersize=20,fmt='s', 
                    markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.set_xlim([550,10000])
    ax.xaxis.set_minor_locator(MultipleLocator(400))   
    ax.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
    ax.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
    ax.set_xlabel(r'$\delta$D$_{\rm H_2O}$' ' + 1000 (\u2030)',fontdict=font_labels)

frac, ax_frac=plt.subplots(figsize=(12,12))
C_dD, ax_C_dD = plt.subplots(figsize=(12,12))
CD_dD, ax_CD_dD = plt.subplots(figsize=(12,12))
DD_dD, ax_DD_dD = plt.subplots(figsize=(12,12))
regAcet,ax_regAcet=plt.subplots(figsize=(12,12))
# plot_original(ax_frac,'dD_CH4')
# plot_original(ax_C_dD,'d13C_CH4')
# plot_original(ax_CD_dD,'D13CH3D')
# plot_original(ax_DD_dD,'D12CH2D2')

mixing_ratio=[1.0, 0.97,0.94,0.92,0.9]
c=['0.0', '0.5', '0.75', '0.85', '1.0']
# c=c[::-1]

for i in range(len(mixing_ratio)):
    iso50=calc_acet(np.mean(H2CO2['d13C_CH4']),np.mean(Acet['d13C_CH4']),
                    np.mean(H2CO2['dD_CH4']),np.mean(Acet['dD_CH4']),
                    np.mean(d13CH3D_H2CO2_50),np.mean(d13CH3D_target_50),
                    np.mean(d12CH2D2_H2CO2_50),np.mean(d12CH2D2_target_50),mixing_ratio[i])
    iso3000=calc_acet(np.mean(H2CO2_3000['d13C_CH4']),np.mean(Acet_3000['d13C_CH4']),
                    np.mean(H2CO2_3000['dD_CH4']),np.mean(Acet_3000['dD_CH4']),
                    np.mean(d13CH3D_H2CO2_3000),np.mean(d13CH3D_target_3000),
                    np.mean(d12CH2D2_H2CO2_3000),np.mean(d12CH2D2_target_3000),mixing_ratio[i])
    iso8000=calc_acet(np.mean(H2CO2_8000['d13C_CH4']),np.mean(Acet_8000['d13C_CH4']),
                    np.mean(H2CO2_8000['dD_CH4']),np.mean(Acet_8000['dD_CH4']),
                    np.mean(d13CH3D_H2CO2_8000),np.mean(d13CH3D_target_8000),
                    np.mean(d12CH2D2_H2CO2_8000),np.mean(d12CH2D2_target_8000),mixing_ratio[i])
    a=c[i]
    if mixing_ratio[i]!=0.92:
        quick_plot2(ax_frac,1)
        quick_plot2(ax_C_dD,0)
        quick_plot2(ax_CD_dD,2)
        quick_plot2(ax_DD_dD,3)
    if mixing_ratio[i]==0.92:
        quick_plot_star(ax_frac,1)
        quick_plot_star(ax_C_dD,0)
        quick_plot_star(ax_CD_dD,2)
        quick_plot_star(ax_DD_dD,3)
        quick_plot_reg(ax_regAcet)
        print(iso50)
        print(iso3000)
        print(iso8000) 


ax_frac.plot(dD_H2O_m,p1[:,1],linewidth=3.5,color='black')
ax_frac.fill_between(dD_H2O_m,p1[:,1]-ep1[:,1],p1[:,1]+ep1[:,1],color='red',alpha=0.3,zorder=-2)
ax_frac.fill_between(dD_H2O_m,p1[:,1]-2*ep1[:,1],p1[:,1]+2*ep1[:,1],color='red',alpha=0.3,zorder=-2)
ax_frac.set_ylabel(r'$\delta D_{\rm CH_4}$'+' (\u2030)',fontdict=font_labels)
ax_frac.yaxis.set_minor_locator(MultipleLocator(40))
# ax_frac.set_ylim([-350,-200])

ax_C_dD.plot(dD_H2O_m,p1[:,0],linewidth=3.5,color='black')
ax_C_dD.fill_between(dD_H2O_m,p1[:,0]-ep1[:,0],p1[:,0]+ep1[:,0],color='red',alpha=0.3,zorder=-2)
ax_C_dD.fill_between(dD_H2O_m,p1[:,0]-2*ep1[:,0],p1[:,0]+2*ep1[:,0],color='red',alpha=0.3,zorder=-2)
ax_C_dD.set_ylabel(r'$\delta^{13}$C$_{\rm CH_4}$'+' (\u2030)',fontdict=font_labels)
ax_C_dD.set_ylim([-100, -30])
ax_C_dD.yaxis.set_minor_locator(MultipleLocator(2))
ax_C_dD.legend(fontsize=25)

ax_CD_dD.plot(dD_H2O_m,p1[:,2],linewidth=3.5,color='black')
ax_CD_dD.fill_between(dD_H2O_m,p1[:,2]-ep1[:,2],p1[:,2]+ep1[:,2],color='red',alpha=0.3,zorder=-2)
ax_CD_dD.fill_between(dD_H2O_m,p1[:,2]-2*ep1[:,2],p1[:,2]+2*ep1[:,2],color='red',alpha=0.3,zorder=-2)
ax_CD_dD.set_ylabel(r'$\Delta^{13}$CH$_3$D'+' (\u2030)',fontdict=font_labels)
ax_CD_dD.set_ylim([-8,3])
ax_CD_dD.yaxis.set_minor_locator(MultipleLocator(0.4))

ax_DD_dD.plot(dD_H2O_m,p1[:,3],linewidth=3.5,color='black')
ax_DD_dD.fill_between(dD_H2O_m,p1[:,3]-ep1[:,3],p1[:,3]+ep1[:,3],color='red',alpha=0.3,zorder=-2)
ax_DD_dD.fill_between(dD_H2O_m,p1[:,3]-2*ep1[:,3],p1[:,3]+2*ep1[:,3],color='red',alpha=0.3,zorder=-2)
ax_DD_dD.set_ylabel(r'$\Delta^{12}$CH$_2$D$_2$'+' (\u2030)',fontdict=font_labels)
ax_DD_dD.yaxis.set_minor_locator(MultipleLocator(20))
#ax_DD_dD.set_ylim([-60,0])

ref=np.linspace(500,9500,100)
Acet_ref=0.0879*ref+616.925
ax_regAcet.plot(ref, Acet_ref, '--', linewidth=5.0, color='red', label=r'$\delta$D$_{\rm CH_4}$+1000=0.0879($\delta$D$_{\rm H_2O}$+1000)+617')
ax_regAcet.set_ylabel(r'$\delta D_{\rm CH_4}$'+' + 1000 (\u2030)',fontdict=font_labels)
ax_regAcet.set_ylim([550,1490])
ax_regAcet.yaxis.set_minor_locator(MultipleLocator(40))
MeOH_ref=0.0798*ref+640.373
ax_regAcet.plot(ref, MeOH_ref, '-.', linewidth=5.0, color='blue', label=r'$\delta$D$_{\rm CH_4}$+1000=0.0798($\delta$D$_{\rm H_2O}$+1000)+640')
CH4_MeOH=np.array([-284.57303723, -72.19256345, 343.77])
H2O_MeOH=np.array([MeOH['dD_H2O'].iloc[0], MeOH_3000['dD_H2O'].iloc[0], MeOH_8000['dD_H2O'].iloc[0]])
eH2O_MeOH=np.array([MeOH['dse_H2O'].iloc[0], MeOH_3000['dse_H2O'].iloc[0], MeOH_8000['dse_H2O'].iloc[0]])
ax_regAcet.errorbar(H2O_MeOH+1000, CH4_MeOH+1000, xerr=eH2O_MeOH,yerr=[0]*3, markersize=20,fmt='o', 
                    markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5,
                    label='Methylotrophic (' +r'$r$ = 0.985'+')')
ax_regAcet.legend(fontsize=22)
plt.show()

# frac.savefig('Acet_model_H.pdf', bbox_inches='tight')
# C_dD.savefig('Acet_model_C.pdf', bbox_inches='tight')
# CD_dD.savefig('Acet_model_CD.pdf', bbox_inches='tight')
# DD_dD.savefig('Acet_model_DD.pdf', bbox_inches='tight')
# regAcet.savefig('reg.pdf', bbox_inches='tight')

