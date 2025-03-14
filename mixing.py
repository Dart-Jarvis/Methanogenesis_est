# This file is to calculate the 'pure' methylotrophic, acetoclastic, 
# #and methoxydotrophic endmember, eliminating the input from hydrogenotrophic methanogenesis
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import *
from matplotlib.ticker import MultipleLocator
# Mixing functions
def clumped_mixing(dD_1, d13C_1, D13CH3D_1, D12CH2D2_1, dD_2, d13C_2, D13CH3D_2, D12CH2D2_2, fraction):
    # Convert it to D/H ratio, then 13C/12C ratio, then to XH, XD, X12C and X13C
    DH_1 = 0.00015576*(dD_1/1000+1)
    DH_2 = 0.00015576*(dD_2/1000+1)
    QC_1 = 0.0112372*(d13C_1/1000+1)
    QC_2 = 0.0112372*(d13C_2/1000+1)
    XD_1 = DH_1/(DH_1+1)
    XD_2 = DH_2/(DH_2+1)
    XH_1 = 1-XD_1
    XH_2 = 1-XD_2
    XQ_1 = QC_1/(QC_1+1)
    XQ_2 = QC_2/(QC_2+1)
    XC_1 = 1-XQ_1
    XC_2 = 1-XQ_2
    # Calculate the two stochastic ratios relative to CH4
    XCH4_1 = XC_1*XH_1**4
    XCH4_2 = XC_2*XH_2**4
    X13CDsto_1 = 4*XQ_1*XD_1*XH_1**3/XCH4_1
    X13CDsto_2 = 4*XQ_2*XD_2*XH_2**3/XCH4_2
    XCDDsto_1 = 6*XC_1*XH_1**2*XD_1**2/XCH4_1
    XCDDsto_2 = 6*XC_2*XH_2**2*XD_2**2/XCH4_2
    # Calculate the real concentration of clumped isotopologues
    C13CD_real_1 = (D13CH3D_1/1000+1)*X13CDsto_1*XCH4_1
    C13CD_real_2 = (D13CH3D_2/1000+1)*X13CDsto_2*XCH4_2
    C12CDD_real_1 = (D12CH2D2_1/1000+1)*XCDDsto_1*XCH4_1
    C12CDD_real_2 = (D12CH2D2_2/1000+1)*XCDDsto_2*XCH4_2
    # Calculate the bulk isotope after mixing, assuming the fraction of gas 1 is f
    D13CH3D = []
    D12CH2D2 = []
    for f in fraction:
        XD = XD_1*f + XD_2*(1-f)
        XH = XH_1*f + XH_2*(1-f)
        XQ = XQ_1*f + XQ_2*(1-f)
        XC = XC_1*f + XC_2*(1-f)
        # Calculate real 13CD concentration after mixing
        C13CD_real = C13CD_real_1*f + C13CD_real_2*(1-f)
        C12CDD_real = C12CDD_real_1*f + C12CDD_real_2*(1-f)
        # Calculate stochastic ratios after mixing
        XCH4 = XCH4_1*f + XCH4_2*(1-f)
        X13CDsto = 4*XQ*XD*XH**3/XCH4
        XCDDsto = 6*XC*XH**2*XD**2/XCH4
        # Calculate Deltas after mixing
        D13CH3D_mix = (C13CD_real/XCH4/X13CDsto-1)*1000
        D12CH2D2_mix = (C12CDD_real/XCH4/XCDDsto-1)*1000
        D13CH3D.append(D13CH3D_mix)
        D12CH2D2.append(D12CH2D2_mix)
    return D13CH3D, D12CH2D2
def bulk_mixing(dD_1,d13C_1,dD_2, d13C_2,f):
    d13C=np.zeros(len(f))
    dD=np.zeros(len(f))
    for i in range(len(f)):
        d13C[i]=d13C_1*f[i]+d13C_2*(1-f[i])
        dD[i]=dD_1*f[i]+dD_2*(1-f[i])
    return d13C,dD
def mix_delta(d13C1,d13C2,dD1,dD2,d13CD1,d13CD2,dDD1,dDD2,fmix):
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
# Calculate small detla 13CD and delta DD from big Deltas
def calc_delta(d13C,dD,D13CD,DDD):
    # Convert it to D/H ratio, then 13C/12C ratio, then to XH, XD, X12C and X13C
    DH= 0.00015576*(dD/1000+1)
    QC= 0.0112372*(d13C/1000+1)
    XD = DH/(DH+1)
    XH = 1-XD
    XQ = QC/(QC+1)
    XC = 1-XQ
    # Calculate the two stochastic ratios relative to CH4
    XCH4 = XC*XH**4
    X13CDsto = 4*XQ*XD*XH**3/XCH4
    XCDDsto = 6*XC*XH**2*XD**2/XCH4
    # Calculate the real concentration of clumped isotopologues
    d13CD = ((D13CD/1000+1)*X13CDsto/6.34600e-6-1)*1000
    dDD = ((DDD/1000+1)*XCDDsto/1.28866e-7-1)*1000
    return d13C, dD, d13CD, dDD
def convert_delta(dataset):
    temp=np.zeros([len(dataset),4])
    for i in range(len(dataset)):
        temp[i,:]=calc_delta(dataset['d13C_CH4'].iloc[i], dataset['dD_CH4'].iloc[i], 
                             dataset['D13CH3D'].iloc[i],dataset['D12CH2D2'].iloc[i])
    return temp
def calc_Delta(d13C,dD,d13CD,dDD):
    # Convert it to D/H ratio, then 13C/12C ratio, then to XH, XD, X12C and X13C
    DH= 0.00015576*(dD/1000+1)
    QC= 0.0112372*(d13C/1000+1)
    XD = DH/(DH+1)
    XH = 1-XD
    XQ = QC/(QC+1)
    XC = 1-XQ
    # Calculate the two stochastic ratios relative to CH4
    XCH4 = XC*XH**4
    X13CDsto = 4*XQ*XD*XH**3/XCH4
    XCDDsto = 6*XC*XH**2*XD**2/XCH4
    # Calculate big Delta
    D13CD=((d13CD/1000+1)*6.34600e-6/X13CDsto-1)*1000
    DDD=((dDD/1000+1)*1.28866e-7/XCDDsto-1)*1000
    return np.array([d13C,dD,D13CD,DDD])

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

plt.rcParams['savefig.facecolor']='white'
data = pd.read_csv('isotope_data.csv')
equib = pd.read_csv('equilibrium.csv')
# Data grouping
MeOH_H2=data[data['Label']=='1HL'] # MeOH + H2
MeOH_LT=data[data['Label']=='1L'] # Low temp MeOH
MeOH=MeOH_LT.iloc[2:5]
TMA=data[data['Label']=='2H'] # Low temp TMA + H2
TMA_NH=data[data['Label']=='2NH'] # Low temp TMA, no H2
TMB=data[data['Label']=='3H'] # High temp TMB (high temp only)
H2CO2=data[data['Label']=='4L'] # Low temp H2/CO2 (35-39)
H2CO2_HT=data[data['Label']=='4H'] # High temp H2/CO2 (65-80)
Acet=data[data['Label']=='5L'] # Low temp acetotrophic (low temp only)
H2CO2_3000=data[data['Label']=='4LS1'] # Low temp H2/CO2, spiked +3000
H2CO2_8000=data[data['Label']=='4LS2'] # Low temp H2/CO2, spiked +8000
MeOH_3000=data[data['Label']=='1LS1'] # Low temp MeOH, spiked +3000
MeOH_8000=data[data['Label']=='1LS2'] # Low temp MeOH, spiked +8000
Acet_3000=data[data['Label']=='5LS1'] # Low temp H2/CO2, spiked +3000
Acet_8000=data[data['Label']=='5LS2'] # Low temp H2/CO2, spiked +8000
Env_P=data[data['Label']=='7P'] # Previous environmental 
Env_coalbed=data[data['Label']=='8P']
# Methylotrophic
# Endmember 1: H2CO2
H2CO2_delta=convert_delta(H2CO2)
H2CO2_avg=np.average(H2CO2_delta,axis=0)
H2CO2_end=calc_Delta(H2CO2_avg[0],H2CO2_avg[1],H2CO2_avg[2],H2CO2_avg[3])

H2CO2_delta2=convert_delta(H2CO2_HT)
H2CO2_avg2=np.average(H2CO2_delta2,axis=0)
H2CO2_end2=calc_Delta(H2CO2_avg2[0],H2CO2_avg2[1],H2CO2_avg2[2],H2CO2_avg2[3])
# Mixture: Methanol
MeOH_delta=convert_delta(MeOH)
MeOH_avg=np.average(MeOH_delta,axis=0)
mix_MeOH_end=calc_Delta(MeOH_avg[0],MeOH_avg[1],MeOH_avg[2],MeOH_avg[3])
# Calculate 'pure' methylotrophic endmember
MeOH_pure=mix_delta(H2CO2_avg[0],MeOH_avg[0],H2CO2_avg[1],MeOH_avg[1],
                    H2CO2_avg[2],MeOH_avg[2],H2CO2_avg[3],MeOH_avg[3],0.985)
# Mixing curves
frac=np.linspace(0,1,21)
D13CD_MeOH, DDD_MeOH=clumped_mixing(MeOH_pure[1],MeOH_pure[0],MeOH_pure[2],MeOH_pure[3],
                                    H2CO2_end[1],H2CO2_end[0],H2CO2_end[2],H2CO2_end[3],frac)
D13CD_MeOH2, DDD_MeOH2=clumped_mixing(MeOH_pure[1],MeOH_pure[0],MeOH_pure[2],MeOH_pure[3],
                                    H2CO2_end2[1],H2CO2_end2[0],H2CO2_end2[2],H2CO2_end2[3],frac)
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Acetoclastic
# Endmember 1: H2CO2
# Mixture: Methanol
Acet_delta=convert_delta(Acet)
Acet_avg=np.average(Acet_delta,axis=0)
Acet_end=calc_Delta(Acet_avg[0],Acet_avg[1],Acet_avg[2],Acet_avg[3])
# Calculate 'pure' methylotrophic endmember
Acet_pure=mix_delta(H2CO2_avg[0],Acet_avg[0],H2CO2_avg[1],Acet_avg[1],
                    H2CO2_avg[2],Acet_avg[2],H2CO2_avg[3],Acet_avg[3],0.92)
D13CD_Acet, DDD_Acet=clumped_mixing(Acet_pure[1],Acet_pure[0],Acet_pure[2],Acet_pure[3],
                                    H2CO2_end[1],H2CO2_end[0],H2CO2_end[2],H2CO2_end[3],frac)
D13CD_Acet2, DDD_Acet2=clumped_mixing(Acet_pure[1],Acet_pure[0],Acet_pure[2],Acet_pure[3],
                                    H2CO2_end2[1],H2CO2_end2[0],H2CO2_end2[2],H2CO2_end2[3],frac)
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Methoxydotrophic
# Endmember 1: High-temp H2/CO2
# Mixture: TMB
TMB_delta=convert_delta(TMB)
TMB_avg=np.average(TMB_delta,axis=0)
# Pure TMB endmember
TMB_pure=mix_delta(H2CO2_avg2[0],TMB_avg[0],H2CO2_avg2[1],TMB_avg[1],
                    H2CO2_avg2[2],TMB_avg[2],H2CO2_avg2[3],TMB_avg[3],0.66)
print(TMB_pure)
# Use M. barkeri as an endmember
TMB_pure2=mix_delta(H2CO2_avg[0],TMB_avg[0],H2CO2_avg[1],TMB_avg[1],
                    H2CO2_avg[2],TMB_avg[2],H2CO2_avg[3],TMB_avg[3],0.66)
# Mixing
D13CD_TMB, DDD_TMB=clumped_mixing(TMB_pure[1],TMB_pure[0],TMB_pure[2],TMB_pure[3],
                                    H2CO2_end2[1],H2CO2_end2[0],H2CO2_end2[2],H2CO2_end2[3],frac)
D13CD_TMB2, DDD_TMB2=clumped_mixing(TMB_pure[1],TMB_pure[0],TMB_pure[2],TMB_pure[3],
                                    H2CO2_end[1],H2CO2_end[0],H2CO2_end[2],H2CO2_end[3],frac)
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Plot
fig,ax=plt.subplots(figsize=(12,12))
ax.plot(equib['D13CH3D'],equib['D12CH2D2'],'-k', label = 'Equilibrium', linewidth = 2.5, markersize = 15)
for i in range(len(equib)):
    if equib['p'].iloc[i]==1:
        ax.scatter(equib['D13CH3D'].iloc[i], equib['D12CH2D2'].iloc[i],color='black',s=60)

# ax.scatter(H2CO2['D13CH3D'],H2CO2['D12CH2D2'])
# ax.scatter(MeOH['D13CH3D'],MeOH['D12CH2D2'])
# ax.scatter(H2CO2_HT['D13CH3D'],H2CO2_HT['D12CH2D2'])
ax.scatter(H2CO2_end[2],H2CO2_end[3], color='yellow',s=240, marker='s', edgecolors='black',
           label=r'Hydrogenotrophic ($M. barkeri$)')
ax.scatter(H2CO2_end2[2], H2CO2_end2[3], color='red', s=240, marker='s', edgecolors='black',
           label=r'Hydrogenotrophic (hyperthermopilic)')
ax.scatter(mix_MeOH_end[2],mix_MeOH_end[3], color='white', s=300, marker='^', edgecolors='black',
           label=r'Methylotrophic (measured)')
ax.scatter(MeOH_pure[2],MeOH_pure[3], color='blue', s=240, marker='^', edgecolors='black',
           label=r'Methylotrophic (pure)')
ax.scatter(Acet_end[2],Acet_end[3], color='white', s=300, marker='D', edgecolors='black',
           label=r'Acetoclastic (measured)')
ax.scatter(Acet_pure[2],Acet_pure[3], color='cyan', s=240, marker='D', edgecolors='black',
           label=r'Acetoclastic (pure)')
ax.scatter(TMB['D13CH3D'],TMB['D12CH2D2'], color='white', s=300, marker='v', edgecolors='black',
           label=r'Methoxydotrophic (measured)')
ax.scatter(TMB_pure[2],TMB_pure[3], color='purple', s=240, marker='v', edgecolors='black',
           label=r'Methoxydotrophic (pure)')
ax.errorbar(Env_P['D13CH3D'],Env_P['D12CH2D2'],xerr=Env_P['13CDse'],yerr=Env_P['DDse'], markersize=10,label=r'Environmental samples (Previous)', fmt='D', 
            markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-3)
ax.errorbar(Env_coalbed['D13CH3D'],Env_coalbed['D12CH2D2'],xerr=Env_coalbed['13CDse'],yerr=Env_coalbed['DDse'], markersize=10,label=r'Incubation with coal (Previous)', fmt='v', 
            markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-3)
# ax.scatter(TMB_pure2[2],TMB_pure2[3], color='purple', s=240, marker='v', edgecolors='black',
#            label=r'Methoxydotrophic 2 (pure)')
ax.plot(D13CD_MeOH, DDD_MeOH, '--',color='blue', zorder=-2, 
        linewidth=2.5, label='Hydrogenotrophic + Methylotrophic')
ax.scatter(D13CD_MeOH[::2], DDD_MeOH[::2],color='blue',s=90,zorder=-2)
ax.plot(D13CD_MeOH2, DDD_MeOH2, '--',color='blue', zorder=-2)
ax.scatter(D13CD_MeOH2[::2], DDD_MeOH2[::2],color='blue',s=90,zorder=-2)
ax.plot(D13CD_Acet, DDD_Acet, '--',color='green', zorder=-2,
        linewidth=2.5, label='Hydrogenotrophic + Acetoclastic')
ax.scatter(D13CD_Acet[::2], DDD_Acet[::2],color='green', s=90,zorder=-2)
ax.plot(D13CD_Acet2, DDD_Acet2, '--',color='green', zorder=-2)
ax.scatter(D13CD_Acet2[::2], DDD_Acet2[::2],color='green', s=90,zorder=-2)
ax.plot(D13CD_TMB, DDD_TMB, '--',color='purple', zorder=-2,
        linewidth=2.5, label='Hydrogenotrophic + Methoxydotrophic')
ax.scatter(D13CD_TMB[::2], DDD_TMB[::2],color='purple', s=90,zorder=-2)
# ax.plot(D13CD_TMB2, DDD_TMB2, '--',color='purple', zorder=-2,
#         linewidth=2.5)
# ax.scatter(D13CD_TMB2[::2], DDD_TMB2[::2],color='purple', s=90,zorder=-2)

ax.set_xlabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
ax.set_ylabel('$\Delta^{12}$CH$_2$D$_2$ (\u2030)', fontdict = font_labels)
ax.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
ax.set_ylim([-60,20])
ax.set_xlim([-9,6])
ax.yaxis.set_minor_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
# ax.legend(fontsize=25, bbox_to_anchor=(1.02,1.00))

# Plot the bulk
fig2,ax2=plt.subplots(figsize=(12,12))
ax2.scatter(H2CO2_end[0],H2CO2_end[1], color='yellow',s=240, marker='s', edgecolors='black',
           label=r'Hydrogenotrophic ($M. barkeri$)')
ax2.scatter(H2CO2_end2[0], H2CO2_end2[1], color='red', s=240, marker='s', edgecolors='black',
           label=r'Hydrogenotrophic (hyperthermopilic)')
ax2.scatter(mix_MeOH_end[0],mix_MeOH_end[1], color='white', s=300, marker='^', edgecolors='black',
           label=r'Methylotrophic (measured)')
ax2.scatter(MeOH_pure[0],MeOH_pure[1], color='blue', s=240, marker='^', edgecolors='black',
           label=r'Methylotrophic (pure)',zorder=2)
ax2.scatter(Acet_end[0],Acet_end[1], color='white', s=300, marker='D', edgecolors='black',
           label=r'Acetoclastic (measured)')
ax2.scatter(Acet_pure[0],Acet_pure[1], color='cyan', s=240, marker='D', edgecolors='black',
           label=r'Acetoclastic (pure)')
ax2.scatter(TMB['d13C_CH4'],TMB['dD_CH4'], color='white', s=300, marker='v', edgecolors='black',
           label=r'Methoxydotrophic (measured)')
ax2.scatter(TMB_pure[0],TMB_pure[1], color='purple', s=240, marker='v', edgecolors='black',
           label=r'Methoxydotrophic (pure)')
ax2.errorbar(Env_P['d13C_CH4'],Env_P['dD_CH4'],xerr=Env_P['cse'],yerr=Env_P['dse'], markersize=10,label=r'Environmental samples (Previous)', fmt='D', 
            markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)
ax2.errorbar(Env_coalbed['d13C_CH4'],Env_coalbed['dD_CH4'],xerr=Env_coalbed['cse'],yerr=Env_coalbed['dse'], markersize=10,label=r'Incubation with coal (Previous)', fmt='v', 
            markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)
# Mixing in bulk
def bulk_mixing(dD_1,d13C_1,dD_2, d13C_2):
    f=np.linspace(0,1,21)
    d13C=np.zeros(len(f))
    dD=np.zeros(len(f))
    for i in range(len(f)):
        d13C[i]=d13C_1*f[i]+d13C_2*(1-f[i])
        dD[i]=dD_1*f[i]+dD_2*(1-f[i])
    return d13C,dD
MeOH_C, MeOH_D=bulk_mixing(H2CO2_end[1],H2CO2_end[0],MeOH_pure[1],MeOH_pure[0])
Acet_C, Acet_D=bulk_mixing(H2CO2_end[1],H2CO2_end[0],Acet_pure[1],Acet_pure[0])
TMB_C, TMB_D=bulk_mixing(H2CO2_end2[1],H2CO2_end2[0],TMB_pure[1],TMB_pure[0])
MeOH_C2, MeOH_D2=bulk_mixing(H2CO2_end2[1],H2CO2_end2[0],MeOH_pure[1],MeOH_pure[0])
Acet_C2, Acet_D2=bulk_mixing(H2CO2_end2[1],H2CO2_end2[0],Acet_pure[1],Acet_pure[0])
ax2.plot(MeOH_C, MeOH_D, '--',color='blue', zorder=-2, 
        linewidth=2.5, label='Hydrogenotrophic + Methylotrophic')
ax2.scatter(MeOH_C[::2], MeOH_D[::2],color='blue',s=90,zorder=-2)
ax2.plot(Acet_C, Acet_D, '--',color='green', zorder=-2,
        linewidth=2.5, label='Hydrogenotrophic + Acetoclastic')
ax2.plot(Acet_C2, Acet_D2, '--',color='green', zorder=-2,
        linewidth=2.5)
ax2.scatter(Acet_C[::2], Acet_D[::2],color='green', s=90,zorder=-2)
ax2.scatter(Acet_C2[::2], Acet_D2[::2],color='green', s=90,zorder=-2)
ax2.plot(TMB_C, TMB_D, '--',color='purple', zorder=-2,
        linewidth=2.5, label='Hydrogenotrophic + Methoxydotrophic')
ax2.scatter(TMB_C[::2], TMB_D[::2],color='purple', s=90,zorder=-2)
ax2.set_xlabel('$\delta^{13}$C (\u2030)', fontdict = font_labels)
ax2.set_ylabel('$\delta$D (\u2030)', fontdict = font_labels)
ax2.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax2.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax2.set_ylim([-70,30])
ax2.set_xlim([-110,-35])
ax2.yaxis.set_minor_locator(MultipleLocator(10))
ax2.xaxis.set_minor_locator(MultipleLocator(4))
ax2.legend(fontsize=25, bbox_to_anchor=(1.02,1.00))

fig.savefig('Grand_unity.pdf',bbox_inches='tight')
fig2.savefig('Grand_unity_bulk.pdf',bbox_inches='tight')