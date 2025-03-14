import math
import numpy
# Define a function to calculate net deltaG, from Gropp et al., 2022
def dGr(Tc,h2,co2,ch4):

# This function takes as imput temperature, CO2, H2 and CH4 concentrations
# in M, and produces a vector of transformed Gibbs free energies (dGr')
# 
# Output:
#       dGr - transformed Gibbs free energy (kJ/mol)
# Input:
#       Tc  - temperature in degree Celsius
#       H2  - H2 concentration (M)
#       CO2 - CO2 concentration (M)
#       CH4 - CH4 concentration (M)
########################################################

# Calculate the dG'0 of methanogenensis at specified temperature using
# Van't Hoffs equation. Values from Alberty, 2003.
    T1  = 298.15
    T2  = Tc + 273.15
    Ha  = [-90.68,-286.65,-413.8,-5.02]
    R   = 8.31446
    dG  = -193.13e3
    dH  = ((Ha[0] + 2*Ha[1]) - (Ha[2] + 4*Ha[3]))*1000
    vh  = (-dH/R)*((1/T2)-(1/T1))
    dg0 = -R*T2*math.log(math.exp((-dG/R/T1)+vh))
    dgr = (dg0 + R*(Tc+273.15)*math.log(ch4/(co2*(h2**4))))/1000
    return dgr

# Henry's Law constant from 
# https://henrys-law.org/henry/index.html
R=8.31446
Hs_hydrogen=7.5e-6 # in mol/m^3/Pa
Hs_ch4=1.5e-5
Hs_co2=3.3e-4

#Define a function for conditional dGr calculations using Henry's law
# T: Temperature/degreeC ; pH2, pCO2: partial pressure in Pa ; 
# nCH4_f: a list containing the final CH4 amounts in umol
# V_headspace: headspace volume in mL

def dGr_con(T,pH2,pCO2,nCH4_f, V_headspace):
    pCH4=pCO2/1000 # Almost no methane at the start
    cH2=pH2*Hs_hydrogen/1000 # Converted to M (mol/L)
    cCO2=Hs_co2*pCO2/1000
    cCH4=Hs_ch4*pCH4/1000
    # Starting CO2/H2 amount, assuming 16 mL of headspace
    nH2=pH2*V_headspace*1e-6/R/(T+273.15)*1e6 
    nCO2=pCO2*V_headspace*1e-6/R/(T+273.15)*1e6 
    dGr_ini=dGr(T,cH2,cCO2,cCH4)
    print('Initial H2, CO2, CH4 concentrations:', cH2, cCO2, cCH4)
    print('inital dG_net:', dGr_ini)
    print('--'*30)

    for i in range(len(nCH4_f)):
        eat=nCH4_f[i]*1 
        # Not considering biomass building, but the dGr is not affected significantly
        nCO2_f=nCO2-eat
        nH2_f=nH2-4*eat
        pCO2_f=nCO2_f*R*(T+273.15)/1e6/(V_headspace*1e-6)
        pH2_f=nH2_f*R*(T+273.15)/1e6/(V_headspace*1e-6)
        pCH4_f=nCH4_f[i]*R*(T+273.15)/1e6/(V_headspace*1e-6)
        cH2_f=pH2_f*Hs_hydrogen/1000 # Converted to M (mol/L)
        cCO2_f=Hs_co2*pCO2_f/1000
        cCH4_f=Hs_ch4*pCH4_f/1000
        dGr_f=dGr(T,cH2_f,cCO2_f,cCH4_f)
        print('final H2, CO2, CH4 concentrations:', cH2_f, cCO2_f, cCH4_f)
        print('fraction of substrate:',cH2_f/cH2*100,cCO2_f/cCO2*100)
        print('final dG_net:', dGr_f)
        print('Average dGr:', (dGr_ini+dGr_f)/2)
        print('--'*30)

# conditions for Dartmouth samples
print('--'*30)
print('--'*30)
print("Dartmouth")
# Starting condition: 25 psi (172369 Pa), 80:20 v/v H2/CO2
T=35 # degree C
pH2=0.8*172369
pCO2=0.2*172369
nCH4_f=[93.2,87.8] # For dartmouth samples
V_headspace=16
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)
# conditions for UMass samples
# Full H2
print("UMass Full H2, 80C")
# Starting condition: 2 atm (2.01e5 Pa), 80:20 v/v H2/CO2
T=80
pH2=0.8*2.01e5
pCO2=0.2*2.01e5
nCH4_f=[322.7, 310.2, 321]
V_headspace=35
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)

# 65 C, M.thermo, full H2
print("UMass Full H2, 65C")
T=65
pH2=0.8*2.01e5
pCO2=0.2*2.01e5
nCH4_f=[376.3, 398.9,376.3,344.0]
V_headspace=35
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)

# 30 mL H2, 80C
print("UMass 30 mL H2, 80C")
# Starting condition: 2 atm (2.01e5 Pa), 30 mL 80:20 v/v H2/CO2
# added at 1 atm, then pressurized with N2/CO2 to 2 atm
T=80
pH2=0.8*1.01e5*30/35 
pCO2=0.2*2.01e5
nCH4_f=[136.6, 193.5]
V_headspace=35
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)

# 30 mL H2, 65C
print("UMass 30 mL H2, 65C")
# Starting condition: 2 atm (2.01e5 Pa), 30 mL 80:20 v/v H2/CO2
# added at 1 atm, then pressurized with N2/CO2 to 2 atm
T=65
pH2=0.8*1.01e5*30/35 
pCO2=0.2*2.01e5
nCH4_f=[202.9]
V_headspace=35
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)

# 10 mL H2
print("UMass 10 mL H2, 80C")
# Starting condition: 2 atm (2.01e5 Pa), 80:20 v/v H2/CO2
T=80
pH2=0.8*1.01e5*10/35
pCO2=0.2*2.01e5
nCH4_f=[62.8,64.0]
V_headspace=35
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)

# 10 mL H2
print("UMass 10 mL H2, 65C")
# Starting condition: 2 atm (2.01e5 Pa), 80:20 v/v H2/CO2
T=65
pH2=0.8*1.01e5*10/35
pCO2=0.2*2.01e5
nCH4_f=[76.5]
V_headspace=35
dGr_con(T,pH2,pCO2,nCH4_f,V_headspace)
print('--'*30)

# pCH4=pCO2/1000 # Almost no methane at the start
# cH2=Hs_hydrogen*pH2/1000 # Converted to M (mol/L)
# cCO2=Hs_co2*pCO2/1000
# cCH4=Hs_ch4*pCH4/1000
# # Starting CO2/H2 amount, assuming 60-25=35 mL of headspace
# nH2=pH2*35e-6/R/(T+273.15)*1e6 
# nCO2=pCO2*35e-6/R/(T+273.15)*1e6 
# print('inital dG_net (80 C, 10 mL):', dGr(T,cH2,cCO2,cCH4))
# print(cH2,cCO2,cCH4)

# # ~340 umol of CH4 produced in each tube 
# # accounting for the CO2 and H2 used for biomass building
# nCH4_f=62.8 # umol
# eat=63*1.05
# nCO2_f=nCO2-eat
# nH2_f=nH2-4*eat
# pCO2_f=nCO2_f*R*(T+273.15)/1e6/35e-6
# pH2_f=nH2_f*R*(T+273.15)/1e6/35e-6
# pCH4_f=nCH4_f*R*(T+273.15)/1e6/35e-6
# cH2_f=pH2_f*Hs_hydrogen/1000 # Converted to M (mol/L)
# cCO2_f=Hs_co2*pCO2_f/1000
# cCH4_f=Hs_ch4*pCH4_f/1000
# print('final dG_net (80 C, 10 mL):', dGr(T,cH2_f,cCO2_f,cCH4_f))
# print('fraction of substrate:',cH2_f/cH2*100,cCO2_f/cCO2*100)

# # 2 mL H2
# print("UMass 2 mL H2 (Thought experiment)")
# # Starting condition: 2 atm (2.01e5 Pa), 80:20 v/v H2/CO2
# T=80
# pH2=0.8*1.01e5*2/35
# pCO2=0.2*2.01e5
# pCH4=pCO2/1000 # Almost no methane at the start
# cH2=Hs_hydrogen*pH2/1000 # Converted to M (mol/L)
# cCO2=Hs_co2*pCO2/1000
# cCH4=Hs_ch4*pCH4/1000
# # Starting CO2/H2 amount, assuming 60-25=35 mL of headspace
# nH2=pH2*35e-6/R/(T+273.15)*1e6 
# nCO2=pCO2*35e-6/R/(T+273.15)*1e6 
# print('inital dG_net (80 C, 2 mL):', dGr(T,cH2,cCO2,cCH4))

# Radboud
# 50 mL medium + 75 mL headspace (1 bar 80/20 H2/CO2)
# print("Radboud H2/CO2")
# # Starting condition: 2 atm (2.01e5 Pa), 80:20 v/v H2/CO2
# T=39
# pH2=0.8*1e5*2.4
# pCO2=0.2*1e5*2.4
# pCH4=pCO2/1000 # Almost no methane at the start
# cH2=Hs_hydrogen*pH2/1000 # Converted to M (mol/L)
# cCO2=Hs_co2*pCO2/1000
# cCH4=Hs_ch4*pCH4/1000
# # Starting CO2/H2 amount, assuming 60-25=35 mL of headspace
# nH2=pH2*75e-6/R/(T+273.15)*1e6 # umol
# nCO2=pCO2*75e-6/R/(T+273.15)*1e6 
# print(nH2,nCO2)
# print('inital dG_net (39 C):', dGr(T,cH2,cCO2,cCH4))

# # accounting for the CO2 and H2 used for biomass building
# nCH4_f=1272 # umol
# eat=1272*1
# nCO2_f=nCO2-eat
# nH2_f=nH2-4*eat
# pCO2_f=nCO2_f*R*(T+273.15)/1e6/75e-6
# pH2_f=nH2_f*R*(T+273.15)/1e6/75e-6
# pCH4_f=nCH4_f*R*(T+273.15)/1e6/75e-6
# cH2_f=pH2_f*Hs_hydrogen/1000 # Converted to M (mol/L)
# cCO2_f=Hs_co2*pCO2_f/1000
# cCH4_f=Hs_ch4*pCH4_f/1000
# print('final dG_net (39 C):', dGr(T,cH2_f,cCO2_f,cCH4_f))

# Calculate dG in Gruen et al., 2018
