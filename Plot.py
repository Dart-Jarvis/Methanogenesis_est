import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import *
from Combi_model_function import MCEP
from Combi_model_function import Plot_MCEP
from Combi_model_function import calc_acet

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
H2CO2_model=pd.read_csv('gropp_model.csv')
# Data grouping
MeOH_H2=data[data['Label']=='1HL'] # MeOH + H2
MeOH_HT=data[data['Label']=='1H'] # High temp MeOH
MeOH_LT=data[data['Label']=='1L'] # Low temp MeOH
TMA=data[data['Label']=='2H'] # Low temp TMA + H2
TMA_NH=data[data['Label']=='2NH'] # Low temp TMA, no H2
TMB=data[data['Label']=='3H'] # High temp TMB (high temp only)
H2CO2=data[data['Label']=='4L'] # Low temp H2/CO2 (35-39)
H2CO2_HT=data[data['Label']=='4H'] # High temp H2/CO2 (65-80)
H2CO2_Out=data[data['Label']=='4HO'] # The outlier from H2/CO2, 65C, M. thermo, full H2 
Acet=data[data['Label']=='5L'] # Low temp acetotrophic (low temp only)
MeOH_P=data[data['Label']=='1P'] # Previous MeOH data
H2CO2_P=data[data['Label']=='4P'] # Previous H2/CO2 data
Acet_P=data[data['Label']=='5P'] # Previous Acetate data
H2CO2_3000=data[data['Label']=='4LS1'] # Low temp H2/CO2, spiked +3000
H2CO2_8000=data[data['Label']=='4LS2'] # Low temp H2/CO2, spiked +8000
MeOH_3000=data[data['Label']=='1LS1'] # Low temp MeOH, spiked +3000
MeOH_8000=data[data['Label']=='1LS2'] # Low temp MeOH, spiked +8000
Acet_3000=data[data['Label']=='5LS1'] # Low temp H2/CO2, spiked +3000
Acet_8000=data[data['Label']=='5LS2'] # Low temp H2/CO2, spiked +8000
Mpn_P=data[data['Label']=='6P'] # Mpn data from Taenzer et al., 2020 GCA
Env_P=data[data['Label']=='7P'] # Previous environmental 
Coal=data[data['Label']=='8P'] # Previous environmental 
# Calculate average isotope values each groups - get the typical isotope value for a group
avg_MeOH_H2=[np.average(MeOH_H2["d13C_CH4"]), np.average(MeOH_H2["dD_CH4"]),
             np.average(MeOH_H2["D13CH3D"]), np.average(MeOH_H2["D12CH2D2"])]
avg_MeOH_HT=[np.average(MeOH_HT["d13C_CH4"]), np.average(MeOH_HT["dD_CH4"]),
             np.average(MeOH_HT["D13CH3D"]), np.average(MeOH_HT["D12CH2D2"])]
avg_MeOH_LT=[np.average(MeOH_LT["d13C_CH4"]), np.average(MeOH_LT["dD_CH4"]),
             np.average(MeOH_LT["D13CH3D"]), np.average(MeOH_LT["D12CH2D2"])]
avg_TMA_H2=[np.average(TMA["d13C_CH4"]), np.average(TMA["dD_CH4"]),
             np.average(TMA["D13CH3D"]), np.average(TMA["D12CH2D2"])]
avg_TMA_NH=[np.average(TMA_NH["d13C_CH4"]), np.average(TMA_NH["dD_CH4"]),
             np.average(TMA_NH["D13CH3D"]), np.average(TMA_NH["D12CH2D2"])]
avg_H2CO2=[np.average(H2CO2["d13C_CH4"]), np.average(H2CO2["dD_CH4"]),
             np.average(H2CO2["D13CH3D"]), np.average(H2CO2["D12CH2D2"])]
avg_H2CO2_HT=[np.average(H2CO2_HT["d13C_CH4"]), np.average(H2CO2_HT["dD_CH4"]),
             np.average(H2CO2_HT["D13CH3D"]), np.average(H2CO2_HT["D12CH2D2"])]
avg_H2CO2_P=[np.average(H2CO2_P["d13C_CH4"]), np.average(H2CO2_P["dD_CH4"]),
             np.average(H2CO2_P["D13CH3D"]), np.average(H2CO2_P["D12CH2D2"])]
# Define where are the spiked experiments, where are non-spiked experiments
# nsH2CO2=[0]
nsMeOH=[0,1]
sH2CO2=[0,1]
sMeOH=[2,3,4]
# Summary plot
summary_clump, ax_clump=plt.subplots(figsize=(12,12))
# Plot equilibrium
ax_clump.plot(equib['D13CH3D'],equib['D12CH2D2'],'-k', label = 'Equilibrium', linewidth = 2.5, markersize = 15)
for i in range(len(equib)):
    if equib['p'].iloc[i]==1:
        ax_clump.scatter(equib['D13CH3D'].iloc[i], equib['D12CH2D2'].iloc[i],color='black',s=60)
# Plot our data
# define a function for quick plotting
def quick_plot(ax_clump,x,y,xerr,yerr):
    ax_clump.errorbar(MeOH_H2[x],MeOH_H2[y],xerr=MeOH_H2[xerr],yerr=MeOH_H2[yerr], markersize=19,label=r'CH$_3$OH + H$_2$', fmt='o', 
                    markerfacecolor='orange', markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=10)
    ax_clump.errorbar(MeOH_HT[x],MeOH_HT[y],xerr=MeOH_HT[xerr],yerr=MeOH_HT[yerr], markersize=19,label=r'CH$_3$OH (65 $^{\rm O}$C)', fmt='o', 
                    markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5)
    ax_clump.errorbar(MeOH_LT[x],MeOH_LT[y],xerr=MeOH_LT[xerr],yerr=MeOH_LT[yerr], markersize=16,label=r'CH$_3$OH (35-39 $^{\rm O}$C)', fmt='o', 
                    markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(TMA[x],TMA[y],xerr=TMA[xerr],yerr=TMA[yerr], markersize=19,label=r'TMA + H$_2$', fmt='^', 
                    markerfacecolor='green', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(TMA_NH[x],TMA_NH[y],xerr=TMA_NH[xerr],yerr=TMA_NH[yerr], markersize=19,label=r'TMA', fmt='^', 
                    markerfacecolor='crimson', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(TMB[x],TMB[y],xerr=TMB[xerr],yerr=TMB[yerr], markersize=20,label=r'TMB', fmt='v', 
                    markerfacecolor='purple', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(H2CO2[x],H2CO2[y],xerr=H2CO2[xerr],yerr=H2CO2[yerr], markersize=19,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                    markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(H2CO2_HT[x],H2CO2_HT[y],xerr=H2CO2_HT[xerr],yerr=H2CO2_HT[yerr], markersize=19,label=r'H$_2$/CO$_2$ (65-80 $^{\rm O}$C)', fmt='s', 
                    markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(H2CO2_Out[x],H2CO2_Out[y],xerr=H2CO2_Out[xerr],yerr=H2CO2_Out[yerr], markersize=19,label=r'H$_2$/CO$_2$ (65 $^{\rm O}$C), outlier', fmt='s', 
                    markerfacecolor='white', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(Acet[x],Acet[y],xerr=Acet[xerr],yerr=Acet[yerr], markersize=19,label=r'Acetate', fmt='D', 
                    markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    # Plot previous data
    ax_clump.errorbar(MeOH_P[x],MeOH_P[y],xerr=MeOH_P[xerr],yerr=MeOH_P[yerr], markersize=12,label=r'CH$_3$OH (Previous)', fmt='o', 
                    markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)
    ax_clump.errorbar(H2CO2_P[x],H2CO2_P[y],xerr=H2CO2_P[xerr],yerr=H2CO2_P[yerr], markersize=12,label=r'H$_2$/CO$_2$ (Previous)', fmt='s', 
                    markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)
    ax_clump.errorbar(Mpn_P[x],Mpn_P[y],xerr=Mpn_P[xerr],yerr=Mpn_P[yerr], markersize=24,label=r'Methylphosphonic acid (Previous)', fmt='*', 
                    markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)
    ax_clump.errorbar(Env_P[x],Env_P[y],xerr=Env_P[xerr],yerr=Env_P[yerr], markersize=12,label=r'Environmental samples (Previous)', fmt='D', 
                    markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)
    ax_clump.errorbar(Coal[x],Coal[y],xerr=Coal[xerr],yerr=Coal[yerr], markersize=12,label=r'Incubation with coal (Previous)', fmt='v', 
                    markerfacecolor='gray', alpha=0.7, markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=-1)

# Quick plot with the necessary end members for mixing
def quick_plot_mix(ax_clump,x,y,xerr,yerr):
    ax_clump.errorbar(MeOH_H2[x],MeOH_H2[y],xerr=MeOH_H2[xerr],yerr=MeOH_H2[yerr], markersize=24,label=r'CH$_3$OH + H$_2$', fmt='*', 
                    markerfacecolor='orange', markeredgecolor='black',markeredgewidth=2.5, ecolor='black', elinewidth=2.5, zorder=10)
    ax_clump.errorbar(MeOH_LT[x],MeOH_LT[y],xerr=MeOH_LT[xerr],yerr=MeOH_LT[yerr], markersize=16,label=r'CH$_3$OH (35-39 $^{\rm O}$C)', fmt='o', 
                    markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(TMA[x],TMA[y],xerr=TMA[xerr],yerr=TMA[yerr], markersize=24,label=r'TMA + H$_2$', fmt='*', 
                    markerfacecolor='green', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(TMA_NH[x],TMA_NH[y],xerr=TMA_NH[xerr],yerr=TMA_NH[yerr], markersize=16,label=r'TMA', fmt='^', 
                    markerfacecolor='crimson', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(H2CO2[x],H2CO2[y],xerr=H2CO2[xerr],yerr=H2CO2[yerr], markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                    markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax_clump.errorbar(H2CO2_HT[x],H2CO2_HT[y],xerr=H2CO2_HT[xerr],yerr=H2CO2_HT[yerr], markersize=16,label=r'H$_2$/CO$_2$ (65-80 $^{\rm O}$C)', fmt='s', 
                    markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)

quick_plot(ax_clump,'D13CH3D','D12CH2D2','13CDse','DDse')
ax_clump.set_xlabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
ax_clump.set_ylabel('$\Delta^{12}$CH$_2$D$_2$ (\u2030)', fontdict = font_labels)
#ax_clump.yaxis.set_major_locator(MultipleLocator(5))
ax_clump.xaxis.set_major_locator(MultipleLocator(5))
ax_clump.yaxis.set_minor_locator(MultipleLocator(5))
ax_clump.xaxis.set_minor_locator(MultipleLocator(1))
ax_clump.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_clump.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
ax_clump.set_ylim([-60,10])
ax_clump.set_xlim([-10.5,5.1])
# ax_clump.set_ylim([-70,30])
# ax_clump.set_xlim([-9,8])
ax_clump.legend(fontsize=25,bbox_to_anchor=(1.02,0.999))
# ax_clump.legend(fontsize=25,bbox_to_anchor=(1.02,0.999))
# Plot bulk isotope values
summary_bulk, ax_bulk=plt.subplots(figsize=(12,12))
quick_plot(ax_bulk,'d13C_CH4','dD_CH4','cse','dse')
ax_bulk.set_xlabel('$\delta^{13}$C (\u2030)', fontdict = font_labels)
ax_bulk.set_ylabel('$\delta$D (\u2030)', fontdict = font_labels)
#ax_clump.yaxis.set_major_locator(MultipleLocator(5))
ax_bulk.xaxis.set_major_locator(MultipleLocator(20))
ax_bulk.yaxis.set_minor_locator(MultipleLocator(10))
ax_bulk.xaxis.set_minor_locator(MultipleLocator(5))
ax_bulk.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_bulk.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
ax_bulk.legend(fontsize=25,bbox_to_anchor=(1.02,0.999))
summary_clump.savefig('summary_clump.pdf',bbox_inches='tight')
summary_bulk.savefig('summary_bulk.pdf',bbox_inches='tight')
# Model the combinatorial effect
# dD_H2O_m=np.array(range(-150,9050,50))
# f=0.32 # Need to change these two when the model changes
# ef=0.03
# p1,ep1,ap1=Plot_MCEP(dD_H2O_m)
# aDp=ap1[::2]
# aDs=ap1[1::2]
# eaDp=std(aDp)
# eaDs=std(aDs)
# aDp=average(aDp)
# aDs=average(aDs)
# Fractionation curves
frac, ax_frac=plt.subplots(figsize=(12,12))
ax_frac.errorbar(H2CO2['dD_H2O'].iloc[sH2CO2]+1000,H2CO2['dD_CH4'].iloc[sH2CO2]+1000,xerr=H2CO2['dse_H2O'].iloc[sH2CO2],yerr=H2CO2['dse'].iloc[sH2CO2],
                markersize=16, fmt='s', 
                markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, label=r'H$_2$/CO$_2$ (${\it M. barkeri}$)')
ax_frac.errorbar(H2CO2_3000['dD_H2O']+1000,H2CO2_3000['dD_CH4']+1000,xerr=H2CO2_3000['dse_H2O'],yerr=H2CO2_3000['dse'],
                markersize=16,fmt='s', 
                markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_frac.errorbar(H2CO2_8000['dD_H2O']+1000,H2CO2_8000['dD_CH4']+1000,xerr=H2CO2_8000['dse_H2O'],yerr=H2CO2_8000['dse'],
                markersize=16,fmt='s', 
                markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_frac.errorbar(MeOH_LT['dD_H2O'].iloc[sMeOH]+1000,MeOH_LT['dD_CH4'].iloc[sMeOH]+1000,xerr=MeOH_LT['dse_H2O'].iloc[sMeOH],yerr=MeOH_LT['dse'].iloc[sMeOH], 
                markersize=16, fmt='o', 
                markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_frac.errorbar(MeOH_3000['dD_H2O']+1000,MeOH_3000['dD_CH4']+1000,xerr=MeOH_3000['dse_H2O'],yerr=MeOH_3000['dse'], 
                markersize=16,fmt='o', 
                markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, label= r'CH$_3$OH (${\it M. barkeri}$)')
ax_frac.errorbar(MeOH_8000['dD_H2O']+1000,MeOH_8000['dD_CH4']+1000,xerr=MeOH_8000['dse_H2O'],yerr=MeOH_8000['dse'], 
                markersize=16,fmt='o', 
                markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_frac.errorbar(MeOH_LT['dD_H2O'].iloc[nsMeOH]+1000,MeOH_LT['dD_CH4'].iloc[nsMeOH]+1000,xerr=MeOH_LT['dse_H2O'].iloc[nsMeOH],yerr=MeOH_LT['dse'].iloc[nsMeOH], 
                markersize=16,fmt='o',
                markerfacecolor='gray', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
ax_frac.errorbar(MeOH_HT['dD_H2O']+1000,MeOH_HT['dD_CH4']+1000,xerr=MeOH_HT['dse_H2O'],yerr=MeOH_HT['dse'], 
                markersize=16,fmt='^',
                markerfacecolor='gray', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
ax_frac.errorbar(Acet['dD_H2O']+1000,Acet['dD_CH4']+1000,xerr=Acet['dse_H2O'],yerr=Acet['dse'],
                markersize=16, fmt='D', 
                markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1, label=r'Acetate (${\it M. barkeri}$)')
ax_frac.errorbar(Acet_3000['dD_H2O']+1000,Acet_3000['dD_CH4']+1000,xerr=Acet_3000['dse_H2O'],yerr=Acet_3000['dse'],
                markersize=16,fmt='D', 
                markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
ax_frac.errorbar(Acet_8000['dD_H2O']+1000,Acet_8000['dD_CH4']+1000,xerr=Acet_8000['dse_H2O'],yerr=Acet_8000['dse'],
                markersize=16,fmt='D', 
                markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
# Plot the reference lines
# MeOH: y=0.0873424x+629.984
# CO2/H2 y=0.469525x+19.0921
ref=np.linspace(950,9500,100)
MeOH_ref=0.0873424*ref+629.984 # +/- 0.0000823 & +/- 0.0898757
H2_ref=0.469525*ref+19.0921 # +/- 
Acet_ref=0.119969*ref+567.767 # +/- 7.066e-05 & +/- 0.0793792 
ax_frac.plot(ref,MeOH_ref,'--b',linewidth=3.5, label=r'$\delta$D$_{\rm CH_4}$+1000=0.0873($\delta$D$_{\rm H_2O}$+1000)+630')
ax_frac.plot(ref,H2_ref,'--y',linewidth=3.5, label=r'$\delta$D$_{\rm CH_4}$+1000=0.470($\delta$D$_{\rm H_2O}$+1000)+19')
ax_frac.plot(ref,Acet_ref,'--c',linewidth=3.5, label=r'$\delta$D$_{\rm CH_4}$+1000=0.12($\delta$D$_{\rm H_2O}$+1000)+568')
# ax_frac.plot(dD_H2O_m+1000,p1[:,1]+1000,linewidth=3.5,color='black',
#              label='Modeled mean', zorder=-2)
# ax_frac.fill_between(dD_H2O_m+1000,p1[:,1]-ep1[:,1]+1000,p1[:,1]+ep1[:,1]+1000,color='red',alpha=0.3,zorder=-2)
# ax_frac.fill_between(dD_H2O_m+1000,p1[:,1]-2*ep1[:,1]+1000,p1[:,1]+2*ep1[:,1]+1000,color='red',alpha=0.3,zorder=-2)
ax_frac.set_xlabel(r'$\delta$D$_{\rm H_2O}$' '+1000 (\u2030)',fontdict=font_labels)
ax_frac.set_ylabel(r'$\delta$D$_{\rm CH_4}$' '+1000 (\u2030)',fontdict=font_labels)
ax_frac.xaxis.set_minor_locator(MultipleLocator(400))
ax_frac.yaxis.set_minor_locator(MultipleLocator(100))
ax_frac.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_frac.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
ax_frac.legend(fontsize=20)
frac.savefig('frac.pdf',bbox_inches='tight')
# define a function to plot isotope-dD_H2O ax:axis, y: the isotope signauture to plot with ys: the error
def quick_plot2(ax,y,ys):
    ax.errorbar(H2CO2['dD_H2O'].iloc[sH2CO2],H2CO2[y].iloc[sH2CO2],xerr=H2CO2['dse_H2O'].iloc[sH2CO2],yerr=H2CO2[ys].iloc[sH2CO2],
                markersize=16,label=r'H$_2$/CO$_2$ (${\it M. barkeri}$)', fmt='s', 
                markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(H2CO2_3000['dD_H2O'],H2CO2_3000[y],xerr=H2CO2_3000['dse_H2O'],yerr=H2CO2_3000[ys],
                    markersize=16,fmt='s', 
                    markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(H2CO2_8000['dD_H2O'],H2CO2_8000[y],xerr=H2CO2_8000['dse_H2O'],yerr=H2CO2_8000[ys],
                    markersize=16,fmt='s', 
                    markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    # ax.errorbar(H2CO2['dD_H2O'].iloc[nsH2CO2],H2CO2[y].iloc[nsH2CO2],xerr=H2CO2['dse_H2O'].iloc[nsH2CO2],yerr=H2CO2[ys].iloc[nsH2CO2],
    #             markersize=16,label=r'H$_2$/CO$_2$ (${\it M. mazei}$)', fmt='s', 
    #             markerfacecolor='gray', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
    ax.errorbar(MeOH_LT['dD_H2O'].iloc[sMeOH],MeOH_LT[y].iloc[sMeOH],xerr=MeOH_LT['dse_H2O'].iloc[sMeOH],yerr=MeOH_LT[ys].iloc[sMeOH], 
                    markersize=16,label=r'CH$_3$OH (${\it M. barkeri}$)', fmt='o', 
                    markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(MeOH_3000['dD_H2O'],MeOH_3000[y],xerr=MeOH_3000['dse_H2O'],yerr=MeOH_3000[ys], 
                    markersize=16,fmt='o', 
                    markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(MeOH_8000['dD_H2O'],MeOH_8000[y],xerr=MeOH_8000['dse_H2O'],yerr=MeOH_8000[ys], 
                    markersize=16,fmt='o', 
                    markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet['dD_H2O'],Acet[y],xerr=Acet['dse_H2O'],yerr=Acet[ys], 
                    markersize=16,label=r'Acetate (${\it M. barkeri}$)', fmt='D', 
                    markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_3000['dD_H2O'],Acet_3000[y],xerr=Acet_3000['dse_H2O'],yerr=Acet_3000[ys], 
                    markersize=16,fmt='D', 
                    markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar(Acet_8000['dD_H2O'],Acet_8000[y],xerr=Acet_8000['dse_H2O'],yerr=Acet_8000[ys], 
                    markersize=16,fmt='D', 
                    markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)  
    ax.errorbar(MeOH_LT['dD_H2O'].iloc[nsMeOH],MeOH_LT[y].iloc[nsMeOH],xerr=MeOH_LT['dse_H2O'].iloc[nsMeOH],yerr=MeOH_LT[ys].iloc[nsMeOH], 
                    markersize=16,label=r'CH$_3$OH (${\it M. mazei}$)', fmt='o', 
                    markerfacecolor='gray', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
    ax.errorbar(MeOH_HT['dD_H2O'],MeOH_HT[y],xerr=MeOH_HT['dse_H2O'],yerr=MeOH_HT[ys], 
                    markersize=16,label=r'CH$_3$OH (65 $^{\rm o}$C)', fmt='^', 
                    markerfacecolor='gray', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)  
    ax.errorbar(MeOH_H2['dD_H2O'],MeOH_H2[y],xerr=MeOH_H2['dse_H2O'],yerr=MeOH_H2[ys], 
                markersize=16,label=r'CH$_3$OH + H$_2$', fmt='v', 
                markerfacecolor='gray', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, zorder=-1)
  
    ax.set_xlim([-450,9000])
    ax.xaxis.set_minor_locator(MultipleLocator(400))   
    ax.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
    ax.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
    ax.set_xlabel(r'$\delta$D$_{\rm H_2O}$' ' (\u2030)',fontdict=font_labels)
D_dD, ax_D_dD = plt.subplots(figsize=(12,12))
C_dD, ax_C_dD = plt.subplots(figsize=(12,12))
CD_dD, ax_CD_dD = plt.subplots(figsize=(12,12))
DD_dD, ax_DD_dD = plt.subplots(figsize=(12,12))
quick_plot2(ax_D_dD,'dD_CH4','dse')
ax_D_dD.set_ylabel(r'$\delta D_{\rm CH_4}$'+' (\u2030)',fontdict=font_labels)
# ax_D_dD.set_ylim([-140,-40])
ax_D_dD.yaxis.set_minor_locator(MultipleLocator(100))
quick_plot2(ax_C_dD,'d13C_CH4','cse')
# ax_C_dD.plot(dD_H2O_m,p1[:,0],linewidth=3.5,color='black')
# ax_C_dD.fill_between(dD_H2O_m,p1[:,0]-ep1[:,0],p1[:,0]+ep1[:,0],color='red',alpha=0.3,zorder=-2)
# ax_C_dD.fill_between(dD_H2O_m,p1[:,0]-2*ep1[:,0],p1[:,0]+2*ep1[:,0],color='red',alpha=0.3,zorder=-2)
ax_C_dD.set_ylabel(r'$\delta^{13}$C$_{\rm CH_4}$'+' (\u2030)',fontdict=font_labels)
ax_C_dD.set_ylim([-140,-40])
ax_C_dD.yaxis.set_minor_locator(MultipleLocator(4))
ax_C_dD.legend(fontsize=25)
quick_plot2(ax_CD_dD,'D13CH3D','13CDse')
# ax_CD_dD.plot(dD_H2O_m,p1[:,2],linewidth=3.5,color='black')
# ax_CD_dD.fill_between(dD_H2O_m,p1[:,2]-ep1[:,2],p1[:,2]+ep1[:,2],color='red',alpha=0.3,zorder=-2)
# ax_CD_dD.fill_between(dD_H2O_m,p1[:,2]-2*ep1[:,2],p1[:,2]+2*ep1[:,2],color='red',alpha=0.3,zorder=-2)
ax_CD_dD.set_ylabel(r'$\Delta^{13}$CH$_3$D'+' (\u2030)',fontdict=font_labels)
ax_CD_dD.yaxis.set_minor_locator(MultipleLocator(0.4))
#ax_CD_dD.legend(fontsize=25)
quick_plot2(ax_DD_dD,'D12CH2D2','DDse')
# ax_DD_dD.plot(dD_H2O_m,p1[:,3],linewidth=3.5,color='black')
# ax_DD_dD.fill_between(dD_H2O_m,p1[:,3]-ep1[:,3],p1[:,3]+ep1[:,3],color='red',alpha=0.3,zorder=-2)
# ax_DD_dD.fill_between(dD_H2O_m,p1[:,3]-2*ep1[:,3],p1[:,3]+2*ep1[:,3],color='red',alpha=0.3,zorder=-2)
ax_DD_dD.set_ylabel(r'$\Delta^{12}$CH$_2$D$_2$'+' (\u2030)',fontdict=font_labels)
ax_DD_dD.yaxis.set_minor_locator(MultipleLocator(5))
#ax_DD_dD.legend(fontsize=25)
C_dD.savefig('C_dD.pdf',bbox_inches='tight')
D_dD.savefig('D_dD.pdf',bbox_inches='tight')
CD_dD.savefig('CD_dD.pdf',bbox_inches='tight')
DD_dD.savefig('DD_dD.pdf',bbox_inches='tight')
plt.show()

# print('Modeled slope: %.4f +/- %.4f' %(aDp*f, sqrt((aDp*ef)**2+(f*eaDp)**2)))
# print('Modeled intercept: %d +/- %d' %(aDs*(1-f)*(-54.113+1000), sqrt(((1-f)*(-54.113+1000)*eaDs)**2+(aDs*(-54.113+1000)*ef)**2+(aDs*(1-f)*0.246)**2)))
# Plot H2/CO2 data against average dGr
# D13CH3D
H2CO2_D13CD,ax_H2CO2_D13CD = plt.subplots(figsize=(12,12))
ax_H2CO2_D13CD.errorbar(H2CO2['dG_a'],H2CO2['D13CH3D'], xerr=H2CO2['dG_range'],yerr=H2CO2['13CDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_D13CD.errorbar(H2CO2_HT['dG_a'],H2CO2_HT['D13CH3D'], xerr=H2CO2_HT['dG_range'],yerr=H2CO2_HT['13CDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (65 - 80 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_D13CD.plot(H2CO2_model['dGr'],H2CO2_model['D13CH3D'], linewidth=4.0, color='blue', 
                    label='Modeled mean values (without Hmd)')
ax_H2CO2_D13CD.fill_between(H2CO2_model['dGr'], H2CO2_model['D13CH3D']-H2CO2_model['s13CD'], H2CO2_model['D13CH3D']+H2CO2_model['s13CD'],
                            color='blue', alpha=0.2)
ax_H2CO2_D13CD.fill_between(H2CO2_model['dGr'], H2CO2_model['D13CH3D']-2*H2CO2_model['s13CD'], H2CO2_model['D13CH3D']+2*H2CO2_model['s13CD'],
                            color='blue',alpha=0.3) # 2-sigma error area
ax_H2CO2_D13CD.set_ylabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
ax_H2CO2_D13CD.set_xlabel(r'$\Delta G_{r}$ (kJ/mol)', fontdict = font_labels)
ax_H2CO2_D13CD.yaxis.set_minor_locator(MultipleLocator(0.4))
ax_H2CO2_D13CD.xaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_D13CD.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_H2CO2_D13CD.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
ax_H2CO2_D13CD.legend(fontsize=18,loc='lower right')
H2CO2_D13CD.savefig('gropp_D13CD.pdf', bbox_inches='tight')
# ax_H2CO2_D13CD.set_ylim([-60,30])
# ax_H2CO2_D13CD.set_xlim([-6,6])

# DD 
H2CO2_DD2,ax_H2CO2_DD2 = plt.subplots(figsize=(12,12))
ax_H2CO2_DD2.errorbar(H2CO2['dG_a'],H2CO2['D12CH2D2'], xerr=H2CO2['dG_range'],yerr=H2CO2['DDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_DD2.errorbar(H2CO2_HT['dG_a'],H2CO2_HT['D12CH2D2'], xerr=H2CO2_HT['dG_range'],yerr=H2CO2_HT['DDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (65 - 80 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_DD2.plot(H2CO2_model['dGr'],H2CO2_model['D12CH2D2'], linewidth=4.0, color='blue', 
                    label='Modeled mean values (without Hmd)')
ax_H2CO2_DD2.fill_between(H2CO2_model['dGr'], H2CO2_model['D12CH2D2']-H2CO2_model['sDD'], H2CO2_model['D12CH2D2']+H2CO2_model['sDD'],
                            color='blue', alpha=0.2)
ax_H2CO2_DD2.fill_between(H2CO2_model['dGr'], H2CO2_model['D12CH2D2']-2*H2CO2_model['sDD'], H2CO2_model['D12CH2D2']+2*H2CO2_model['sDD'],
                            color='blue', alpha=0.3) # 2-sigma error area
ax_H2CO2_DD2.set_ylabel('$\Delta^{12}$CH$_2$D$_2$ (\u2030)', fontdict = font_labels)
ax_H2CO2_DD2.set_xlabel(r'$\Delta G_{r}$ (kJ/mol)', fontdict = font_labels)
ax_H2CO2_DD2.yaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_DD2.xaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_DD2.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_H2CO2_DD2.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax_H2CO2_DD2.set_ylim([-60,30])
# ax_H2CO2_DD2.set_xlim([-6,6])
# ax_H2CO2_DD2.legend(fontsize=25, bbox_to_anchor=(1.02,0.999))
H2CO2_DD2.savefig('gropp_DD.pdf', bbox_inches='tight')

# epsilon_CH4_H2O
H2CO2_eD,ax_H2CO2_eD = plt.subplots(figsize=(12,12))
ax_H2CO2_eD.errorbar(H2CO2['dG_a'],H2CO2['eps_CH4_H2O'], xerr=H2CO2['dG_range'],yerr=H2CO2['eps_se'],
                        markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_eD.errorbar(H2CO2_HT['dG_a'],H2CO2_HT['eps_CH4_H2O'], xerr=H2CO2_HT['dG_range'],yerr=H2CO2_HT['eps_se'],
                        markersize=16,label=r'H$_2$/CO$_2$ (65 - 80 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_eD.plot(H2CO2_model['dGr'],H2CO2_model['eD'], linewidth=4.0, color='blue', 
                    label='Modeled mean values (without Hmd)')
ax_H2CO2_eD.fill_between(H2CO2_model['dGr'], H2CO2_model['eD']-H2CO2_model['sD'], H2CO2_model['eD']+H2CO2_model['sD'],
                            color='blue', alpha=0.2)
ax_H2CO2_eD.fill_between(H2CO2_model['dGr'], H2CO2_model['eD']-2*H2CO2_model['sD'], H2CO2_model['eD']+2*H2CO2_model['sD'],
                            color='blue', alpha=0.3) # 2-sigma error area
ax_H2CO2_eD.set_ylabel(r'$^{2}\epsilon_{\rm CH_4 - H_2O}$' + '(\u2030)', fontdict = font_labels)
ax_H2CO2_eD.set_xlabel(r'$\Delta G_{r}$ (kJ/mol)', fontdict = font_labels)
ax_H2CO2_eD.yaxis.set_minor_locator(MultipleLocator(20))
ax_H2CO2_eD.xaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_eD.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_H2CO2_eD.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
H2CO2_eD.savefig('gropp_eD.pdf', bbox_inches='tight')

# Plot H2/CO2 data against average dGr (with Hmd model)
H2CO2_model=pd.read_csv('gropp_model_H.csv')
# D13CH3D
H2CO2_D13CD_H,ax_H2CO2_D13CD_H = plt.subplots(figsize=(12,12))
ax_H2CO2_D13CD_H.errorbar(H2CO2['dG_a'],H2CO2['D13CH3D'], xerr=H2CO2['dG_range'],yerr=H2CO2['13CDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_D13CD_H.errorbar(H2CO2_HT['dG_a'],H2CO2_HT['D13CH3D'], xerr=H2CO2_HT['dG_range'],yerr=H2CO2_HT['13CDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (65 - 80 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_D13CD_H.plot(H2CO2_model['dGr'],H2CO2_model['D13CH3D'], linewidth=4.0, color='black', 
                    label='Modeled mean values (with Hmd)')
ax_H2CO2_D13CD_H.fill_between(H2CO2_model['dGr'], H2CO2_model['D13CH3D']-H2CO2_model['s13CD'], H2CO2_model['D13CH3D']+H2CO2_model['s13CD'],
                            color='red', alpha=0.2)
ax_H2CO2_D13CD_H.fill_between(H2CO2_model['dGr'], H2CO2_model['D13CH3D']-2*H2CO2_model['s13CD'], H2CO2_model['D13CH3D']+2*H2CO2_model['s13CD'],
                            color='red',alpha=0.3) # 2-sigma error area
ax_H2CO2_D13CD_H.set_ylabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
ax_H2CO2_D13CD_H.set_xlabel(r'$\Delta G_{r}$ (kJ/mol)', fontdict = font_labels)
ax_H2CO2_D13CD_H.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_H2CO2_D13CD_H.xaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_D13CD_H.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_H2CO2_D13CD_H.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax_H2CO2_D13CD.set_ylim([-60,30])
# ax_H2CO2_D13CD.set_xlim([-6,6])
ax_H2CO2_D13CD_H.legend(fontsize=18, loc='lower right')
H2CO2_D13CD_H.savefig('gropp_D13CD_Hmd.pdf', bbox_inches='tight')
# DD 
H2CO2_DD2_H,ax_H2CO2_DD2_H = plt.subplots(figsize=(12,12))
ax_H2CO2_DD2_H.errorbar(H2CO2['dG_a'],H2CO2['D12CH2D2'], xerr=H2CO2['dG_range'],yerr=H2CO2['DDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_DD2_H.errorbar(H2CO2_HT['dG_a'],H2CO2_HT['D12CH2D2'], xerr=H2CO2_HT['dG_range'],yerr=H2CO2_HT['DDse'],
                        markersize=16,label=r'H$_2$/CO$_2$ (65 - 80 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_DD2_H.plot(H2CO2_model['dGr'],H2CO2_model['D12CH2D2'], linewidth=4.0, color='black', 
                    label='Modeled mean values (with Hmd)')
ax_H2CO2_DD2_H.fill_between(H2CO2_model['dGr'], H2CO2_model['D12CH2D2']-H2CO2_model['sDD'], H2CO2_model['D12CH2D2']+H2CO2_model['sDD'],
                            color='red', alpha=0.2)
ax_H2CO2_DD2_H.fill_between(H2CO2_model['dGr'], H2CO2_model['D12CH2D2']-2*H2CO2_model['sDD'], H2CO2_model['D12CH2D2']+2*H2CO2_model['sDD'],
                            color='red', alpha=0.3) # 2-sigma error area
ax_H2CO2_DD2_H.set_ylabel('$\Delta^{12}$CH$_2$D$_2$ (\u2030)', fontdict = font_labels)
ax_H2CO2_DD2_H.set_xlabel(r'$\Delta G_{r}$ (kJ/mol)', fontdict = font_labels)
ax_H2CO2_DD2_H.yaxis.set_minor_locator(MultipleLocator(4))
ax_H2CO2_DD2_H.xaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_DD2_H.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_H2CO2_DD2_H.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax_H2CO2_DD2.set_ylim([-60,30])
# ax_H2CO2_DD2.set_xlim([-6,6])
# ax_H2CO2_DD2_H.legend(fontsize=25, bbox_to_anchor=(1.02,0.999))
H2CO2_DD2_H.savefig('gropp_DD2_Hmd.pdf', bbox_inches='tight')

# epsilon_CH4_H2O
def a_H2O_CH4(T):
    a=1.0997+8456/T**2+0.9611*1e9/T**4-27.82*1e12/T**6 # in degree K, from Hribe and Craig 1995
    return a
H2CO2_eD_H,ax_H2CO2_eD_H = plt.subplots(figsize=(12,12))
ax_H2CO2_eD_H.errorbar(H2CO2['dG_a'],H2CO2['eps_CH4_H2O'], xerr=H2CO2['dG_range'],yerr=H2CO2['eps_se'],
                        markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_eD_H.errorbar(H2CO2_HT['dG_a'],H2CO2_HT['eps_CH4_H2O'], xerr=H2CO2_HT['dG_range'],yerr=H2CO2_HT['eps_se'],
                        markersize=16,label=r'H$_2$/CO$_2$ (65 - 80 $^{\rm O}$C)', fmt='s', 
                        markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
ax_H2CO2_eD_H.plot(H2CO2_model['dGr'],H2CO2_model['eD'], linewidth=4.0, color='black', 
                    label='Modeled mean values (with Hmd)')
ax_H2CO2_eD_H.fill_between(H2CO2_model['dGr'], H2CO2_model['eD']-H2CO2_model['sD'], H2CO2_model['eD']+H2CO2_model['sD'],
                            color='red', alpha=0.2)
ax_H2CO2_eD_H.fill_between(H2CO2_model['dGr'], H2CO2_model['eD']-2*H2CO2_model['sD'], H2CO2_model['eD']+2*H2CO2_model['sD'],
                            color='red', alpha=0.3) # 2-sigma error area
ax_H2CO2_eD_H.set_ylabel(r'$^{2}\epsilon_{\rm CH_4 - H_2O}$' + '(\u2030)', fontdict = font_labels)
ax_H2CO2_eD_H.set_xlabel(r'$\Delta G_{r}$ (kJ/mol)', fontdict = font_labels)
ax_H2CO2_eD_H.yaxis.set_minor_locator(MultipleLocator(20))
ax_H2CO2_eD_H.xaxis.set_minor_locator(MultipleLocator(10))
ax_H2CO2_eD_H.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
ax_H2CO2_eD_H.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
H2CO2_eD_H.savefig('gropp_eD_Hmd.pdf', bbox_inches='tight')
# Fraction factors 
fig1,ax1=plt.subplots(figsize=(12,6))
fig2,ax2=plt.subplots(figsize=(12,6))
fig3,ax3=plt.subplots(figsize=(6,6))
def plot_alpha(ax,y, yse):
    ax.errorbar([r'H$_2$/CO$_2$ (35-39 $^{\rm O}$C)']*len(H2CO2[y]),H2CO2[y], yerr=H2CO2[yse],
                markersize=16, fmt='s', label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)',
                markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'H$_2$/CO$_2$ (65-80 $^{\rm O}$C)']*len(H2CO2_HT[y]),H2CO2_HT[y], yerr=H2CO2_HT[yse],
                markersize=16, fmt='s', label=r'H$_2$/CO$_2$ (65-80 $^{\rm O}$C)',
                markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'CH$_3$OH (35-39 $^{\rm O}$C)']*len(MeOH_LT[y]), MeOH_LT[y],yerr=MeOH_LT[yse], 
                markersize=16, fmt='o', label=r'CH$_3$OH (35-39 $^{\rm O}$C)',
                markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'CH$_3$OH (65 $^{\rm O}$C)']*len(MeOH_HT[y]), MeOH_HT[y],yerr=MeOH_HT[yse], 
                markersize=16, fmt='o', label=r'CH$_3$OH (65 $^{\rm O}$C)',
                markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'CH$_3$OH + H$_2$']*len(MeOH_H2[y]), MeOH_H2[y],yerr=MeOH_H2[yse], 
                markersize=16, fmt='o', label=r'CH$_3$OH + H$_2$',
                markerfacecolor='orange', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'TMA']*len(TMA_NH[y]), TMA_NH[y],yerr=TMA_NH[yse], 
                markersize=16, fmt='^', label='TMA',
                markerfacecolor='crimson', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'TMA + H$_2$']*len(TMA[y]), TMA[y],yerr=TMA[yse], 
                markersize=16, fmt='^', label=r'TMA + H$_2$',
                markerfacecolor='green', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'TMB']*len(TMB[y]), TMB[y],yerr=TMB[yse], 
                markersize=16, fmt='v', label=r'TMB',
                markerfacecolor='purple', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
    ax.errorbar([r'Acetate']*len(Acet[y]), Acet[y],yerr=Acet[yse], 
                markersize=16, fmt='D', label=r'Acetate',
                markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)

ax3.errorbar([r'Full H$_2$']*len(H2CO2_HT['d13C_CH4'].iloc[0:5]),H2CO2_HT['d13C_CH4'].iloc[0:5], yerr=H2CO2_HT['cse'].iloc[0:5],
            markersize=16, fmt='s', label=r'Full H$_2$',
            markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, alpha=1.0)
ax3.errorbar([r'30 mL H$_2$']*len(H2CO2_HT['d13C_CH4'].iloc[5:8]),H2CO2_HT['d13C_CH4'].iloc[5:8], yerr=H2CO2_HT['cse'].iloc[5:8],
            markersize=16, fmt='s', label=r'30 mL H$_2$',
            markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, alpha=1.0)
ax3.errorbar([r'10 mL H$_2$']*len(H2CO2_HT['d13C_CH4'].iloc[8:11]),H2CO2_HT['d13C_CH4'].iloc[8:11], yerr=H2CO2_HT['cse'].iloc[8:11],
            markersize=16, fmt='s', label=r'10 mL H$_2$',
            markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5, alpha=1.0)

plot_alpha(ax1,'alphaH','alphaH_se')
plot_alpha(ax2,'alphaC','alphaC_se')
# x=[r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', r'H$_2$/CO$_2$ (65-80 $^{\rm O}$C)', r'CH$_3$OH (35-39 $^{\rm O}$C)',
#    r'CH$_3$OH (65 $^{\rm O}$C)', r'CH$_3$OH + H$_2$', r'TMA', r'TMA + H$_2$', r'TMB', r'Acetate']
x=np.linspace(0,9,9)
a35=a_H2O_CH4(35+273.15)
a65=a_H2O_CH4(65+273.15)
a80=a_H2O_CH4(80+273.15)
ax1.plot(x,[1/a35]*len(x), '--', linewidth=3.0, color='blue', label=r'Equilibrium (35 $^{\rm O}$C)')
# ax1.plot(x,[1/a65]*len(x), '--', linewidth=3.0, color='orange')
ax1.plot(x,[1/a80]*len(x), '--', linewidth=3.0, color='red', label=r'Equilibrium (80 $^{\rm O}$C)')
ax1.set_ylabel(r'$^{\rm D}\alpha_{\rm CH4-H2O}$', fontdict=font_labels)
ax2.set_ylabel(r'$^{13}\alpha_{\rm CH4-Carbon}$', fontdict=font_labels)
ax3.set_ylabel(r'$\delta^{13}{\rm C}$'+' (\u2030)', fontdict=font_labels)
ax1.set_ylim([0.45,0.9])
ax2.set_ylim([0.918,0.99])
ax1.yaxis.set_minor_locator(MultipleLocator(0.02))
ax2.yaxis.set_minor_locator(MultipleLocator(0.004))
ax3.yaxis.set_minor_locator(MultipleLocator(2))
ax1.tick_params(axis='x',which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=18)
ax1.tick_params(axis='y',which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=28)
ax1.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=18)
ax2.tick_params(axis='x',which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=18)
ax2.tick_params(axis='y',which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=28)
ax2.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=18)
ax3.tick_params(axis='x',which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=28)
ax3.tick_params(axis='y',which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=28)
ax3.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=18)
ax1.legend(fontsize=24, bbox_to_anchor=(1.02,1.10))
# ax3.legend(fontsize=24)
fig1.savefig('alphaH.pdf',bbox_inches='tight')
fig2.savefig('alphaC.pdf',bbox_inches='tight')
fig3.savefig('dC_H2.pdf',bbox_inches='tight')
plt.show()

# def plot_alpha2(ax, x, y, xse, yse):
#     ax.errorbar(H2CO2[x],H2CO2[y], xerr=H2CO2[xse], yerr=H2CO2[yse],
#                 markersize=16, fmt='s', 
#                 markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(H2CO2_HT[x],H2CO2_HT[y], xerr=H2CO2_HT[xse], yerr=H2CO2_HT[yse],
#                 markersize=16, fmt='s', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(MeOH_LT[x], MeOH_LT[y],xerr=MeOH_LT[xse],yerr=MeOH_LT[yse], 
#                 markersize=16, fmt='o', 
#                 markerfacecolor='blue', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(MeOH_HT[x], MeOH_HT[y],xerr=MeOH_HT[xse],yerr=MeOH_HT[yse], 
#                 markersize=16, fmt='o', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(MeOH_H2[x], MeOH_H2[y],xerr=MeOH_H2[xse],yerr=MeOH_H2[yse], 
#                 markersize=16, fmt='o', 
#                 markerfacecolor='orange', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(TMA_NH[x], TMA_NH[y],xerr=TMA_NH[xse], yerr=TMA_NH[yse], 
#                 markersize=16, fmt='^', 
#                 markerfacecolor='crimson', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(TMA[x], TMA[y],xerr=TMA[xse],yerr=TMA[yse], 
#                 markersize=16, fmt='^', 
#                 markerfacecolor='green', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(TMB[x], TMB[y],xerr=TMB[xse],yerr=TMB[yse], 
#                 markersize=16, fmt='v', 
#                 markerfacecolor='purple', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
#     ax.errorbar(Acet[x], Acet[y],xerr=Acet[xse], yerr=Acet[yse], 
#                 markersize=16, fmt='D', 
#                 markerfacecolor='cyan', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# plot_alpha2(ax3,'alphaC', 'alphaH','alphaC_se','alphaH_se')

# # Mixing scenario 2, mixing between MeOH and high-temp H2/CO2
# d13C_1, dD_1, D13CH3D_1, D12CH2D2_1 = avg_MeOH_LT
# d13C_2, dD_2, D13CH3D_2, D12CH2D2_2 = avg_H2CO2_HT
# D13CH3D_msc, D12CH2D2_msc = clumped_mixing(dD_1, d13C_1, D13CH3D_1, D12CH2D2_1, dD_2, d13C_2, D13CH3D_2, D12CH2D2_2)
# d13C_msc, dD_msc = bulk_mixing(dD_1, d13C_1,dD_2, d13C_2)
# ax_msc_clump.plot(D13CH3D_msc,D12CH2D2_msc, linewidth=3.0, color='red',
#                   label=r"Mixing curve: CH$_3$OH + H$_2$/CO$_2$(high-temp)")
# ax_msc_bulk.plot(d13C_msc,dD_msc, linewidth=3.0, color='red', 
#                  label=r"Mixing curve: CH$_3$OH + H$_2$/CO$_2$(high-temp)")

# # Mixing scenario 3, mixing between TMA and low-temp H2/CO2
# d13C_1, dD_1, D13CH3D_1, D12CH2D2_1 = avg_TMA_NH
# d13C_2, dD_2, D13CH3D_2, D12CH2D2_2 = avg_H2CO2
# D13CH3D_msc, D12CH2D2_msc = clumped_mixing(dD_1, d13C_1, D13CH3D_1, D12CH2D2_1, dD_2, d13C_2, D13CH3D_2, D12CH2D2_2)
# d13C_msc, dD_msc = bulk_mixing(dD_1, d13C_1,dD_2, d13C_2)
# ax_msc_clump.plot(D13CH3D_msc,D12CH2D2_msc, linewidth=3.0, ls='--', color='blue',
#                   label=r"Mixing curve: TMA + H$_2$/CO$_2$(low-temp)")
# ax_msc_bulk.plot(d13C_msc,dD_msc, linewidth=3.0, ls='--', color='blue', 
#                  label=r"Mixing curve: TMA + H$_2$/CO$_2$(low-temp)")

# # Mixing scenario 4, mixing between TMA and high-temp H2CO2
# d13C_1, dD_1, D13CH3D_1, D12CH2D2_1 = avg_TMA_NH
# d13C_2, dD_2, D13CH3D_2, D12CH2D2_2 = avg_H2CO2_HT
# D13CH3D_msc, D12CH2D2_msc = clumped_mixing(dD_1, d13C_1, D13CH3D_1, D12CH2D2_1, dD_2, d13C_2, D13CH3D_2, D12CH2D2_2)
# d13C_msc, dD_msc = bulk_mixing(dD_1, d13C_1,dD_2, d13C_2)
# ax_msc_clump.plot(D13CH3D_msc,D12CH2D2_msc, linewidth=3.0, ls='--', color='red',
#                   label=r"Mixing curve: TMA + H$_2$/CO$_2$(high-temp)")
# ax_msc_bulk.plot(d13C_msc,dD_msc, linewidth=3.0, ls='--', color='red', 
#                  label=r"Mixing curve: TMA + H$_2$/CO$_2$(high-temp)")

# ax_msc_clump.legend(fontsize=20, bbox_to_anchor=(1.02,0.999))
# ax_msc_bulk.legend(fontsize=25,bbox_to_anchor=(1.02,0.999))
# plt.show()
# msc_bulk.savefig('mixing_bulk.pdf',bbox_inches='tight')
# msc_clump.savefig('mixing_clump.pdf',bbox_inches='tight')

# Plot only the H2/CO2
# H2CO2_model,ax_H2CO2_model = plt.subplots(figsize=(12,12))
# ax_H2CO2_model.plot(equib['D13CH3D'],equib['D12CH2D2'],'-k', label = 'Equilibrium', linewidth = 2.5, markersize = 15)
# for i in range(len(equib)):
#     if equib['p'].iloc[i]==1:
#         ax_H2CO2_model.scatter(equib['D13CH3D'].iloc[i], equib['D12CH2D2'].iloc[i],color='black',s=60)
# ax_H2CO2_model.errorbar(H2CO2['D13CH3D'].iloc[1:3],H2CO2['D12CH2D2'].iloc[1:3],xerr=H2CO2['13CDse'].iloc[1:3],yerr=H2CO2['DDse'].iloc[1:3], markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
#                 markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_model.errorbar(H2CO2_HT['D13CH3D'].iloc[0:2],H2CO2_HT['D12CH2D2'].iloc[0:2],xerr=H2CO2_HT['13CDse'].iloc[0:2],yerr=H2CO2_HT['DDse'].iloc[0:2], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (full H$_2$, 80 $^{\rm O}$C)', fmt='s', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_model.errorbar(H2CO2_HT['D13CH3D'].iloc[2:4],H2CO2_HT['D12CH2D2'].iloc[2:4],xerr=H2CO2_HT['13CDse'].iloc[2:4],yerr=H2CO2_HT['DDse'].iloc[2:4], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (full H$_2$, 65 $^{\rm O}$C)', fmt='o', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_model.errorbar(H2CO2_HT['D13CH3D'].iloc[4:6],H2CO2_HT['D12CH2D2'].iloc[4:6],xerr=H2CO2_HT['13CDse'].iloc[4:6],yerr=H2CO2_HT['DDse'].iloc[4:6], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (30 mL H$_2$, 80 $^{\rm O}$C)', fmt='^', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_model.errorbar(H2CO2_HT['D13CH3D'].iloc[6],H2CO2_HT['D12CH2D2'].iloc[6],xerr=H2CO2_HT['13CDse'].iloc[6],yerr=H2CO2_HT['DDse'].iloc[6], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (30 mL H$_2$, 65 $^{\rm O}$C)', fmt='d', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_model.errorbar(H2CO2_HT['D13CH3D'].iloc[7],H2CO2_HT['D12CH2D2'].iloc[7],xerr=H2CO2_HT['13CDse'].iloc[7],yerr=H2CO2_HT['DDse'].iloc[7], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (10 mL H$_2$, 80 $^{\rm O}$C)', fmt='v', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_model.errorbar(H2CO2_HT['D13CH3D'].iloc[8],H2CO2_HT['D12CH2D2'].iloc[8],xerr=H2CO2_HT['13CDse'].iloc[8],yerr=H2CO2_HT['DDse'].iloc[8], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (10 mL H$_2$, 65 $^{\rm O}$C)', fmt='>', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)

# ax_H2CO2_model.set_xlabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
# ax_H2CO2_model.set_ylabel('$\Delta^{12}$CH$_2$D$_2$ (\u2030)', fontdict = font_labels)
# ax_H2CO2_model.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
# ax_H2CO2_model.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax_H2CO2_model.set_ylim([-60,30])
# ax_H2CO2_model.set_xlim([-6,6])
# ax_H2CO2_model.legend(fontsize=25)

# # Meaured by Andy Masterson
# # d13C_CO2 = -36.01±0.05‰ vs. VPDB
# # d2H_H2 = -261.0±5‰ vs. VSMOW

# # epsilon_CH4_H2O vs. D13CH3D
# H2CO2_dD,ax_H2CO2_dD = plt.subplots(figsize=(12,12))
# ax_H2CO2_dD.errorbar(H2CO2['epsilon_CH4_H2O'].iloc[1:3],H2CO2['D13CH3D'].iloc[1:3],xerr=H2CO2['13CDse'].iloc[1:3],yerr=H2CO2['DDse'].iloc[1:3], markersize=16,label=r'H$_2$/CO$_2$ (35 $^{\rm O}$C)', fmt='s', 
#                 markerfacecolor='yellow', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.errorbar(H2CO2_HT['epsilon_CH4_H2O'].iloc[0:2],H2CO2_HT['D13CH3D'].iloc[0:2],xerr=H2CO2_HT['13CDse'].iloc[0:2],yerr=H2CO2_HT['13CDse'].iloc[0:2], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (full H$_2$, 80 $^{\rm O}$C)', fmt='s', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.errorbar(H2CO2_HT['epsilon_CH4_H2O'].iloc[2:4],H2CO2_HT['D13CH3D'].iloc[2:4],xerr=H2CO2_HT['13CDse'].iloc[2:4],yerr=H2CO2_HT['13CDse'].iloc[2:4], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (full H$_2$, 65 $^{\rm O}$C)', fmt='o', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.errorbar(H2CO2_HT['epsilon_CH4_H2O'].iloc[4:6],H2CO2_HT['D13CH3D'].iloc[4:6],xerr=H2CO2_HT['13CDse'].iloc[4:6],yerr=H2CO2_HT['13CDse'].iloc[4:6], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (30 mL H$_2$, 80 $^{\rm O}$C)', fmt='^', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.errorbar(H2CO2_HT['epsilon_CH4_H2O'].iloc[6],H2CO2_HT['D13CH3D'].iloc[6],xerr=H2CO2_HT['13CDse'].iloc[6],yerr=H2CO2_HT['13CDse'].iloc[6], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (30 mL H$_2$, 65 $^{\rm O}$C)', fmt='d', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.errorbar(H2CO2_HT['epsilon_CH4_H2O'].iloc[7],H2CO2_HT['D13CH3D'].iloc[7],xerr=H2CO2_HT['13CDse'].iloc[7],yerr=H2CO2_HT['13CDse'].iloc[7], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (10 mL H$_2$, 80 $^{\rm O}$C)', fmt='v', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.errorbar(H2CO2_HT['epsilon_CH4_H2O'].iloc[8],H2CO2_HT['D13CH3D'].iloc[8],xerr=H2CO2_HT['13CDse'].iloc[8],yerr=H2CO2_HT['13CDse'].iloc[8], 
#                 markersize=16,label=r'H$_2$/CO$_2$ (10 mL H$_2$, 65 $^{\rm O}$C)', fmt='>', 
#                 markerfacecolor='red', markeredgecolor='black',markeredgewidth=2.5,ecolor='black',elinewidth=2.5)
# ax_H2CO2_dD.set_xlabel(r'$\epsilon_{\rm CH_4 - H_2O}$' + '(\u2030)', fontdict = font_labels)
# ax_H2CO2_dD.set_ylabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
# ax_H2CO2_dD.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
# ax_H2CO2_dD.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax_H2CO2_dD.set_ylim([-5,6])
# ax_H2CO2_dD.set_xlim([-550,-100])
# plt.show()

# Mixing curves
def clumped_mixing(dD_1, d13C_1, D13CH3D_1, D12CH2D2_1, dD_2, d13C_2, D13CH3D_2, D12CH2D2_2):
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
    fraction = np.linspace(0,1,21) # This controls the f you want to plot with
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
def bulk_mixing(dD_1,d13C_1,dD_2, d13C_2):
    f=np.linspace(0,1,21)
    d13C=np.zeros(len(f))
    dD=np.zeros(len(f))
    for i in range(len(f)):
        d13C[i]=d13C_1*f[i]+d13C_2*(1-f[i])
        dD[i]=dD_1*f[i]+dD_2*(1-f[i])
    return d13C,dD

# # Mixing scenario 1: mixing with MeOH and low temp H2CO2 
# d13C_1, dD_1, D13CH3D_1, D12CH2D2_1 = avg_TMA_H2
# d13C_2, dD_2, D13CH3D_2, D12CH2D2_2 = avg_H2CO2_HT
# D13CH3D_msc_clump, D12CH2D2_msc_clump = clumped_mixing(dD_1, d13C_1, D13CH3D_1, D12CH2D2_1, dD_2, d13C_2, D13CH3D_2, D12CH2D2_2)
# d13C_msc, dD_msc = bulk_mixing(dD_1, d13C_1,dD_2, d13C_2)
# print(D13CH3D_msc_clump,D12CH2D2_msc_clump)
# msc_clump, ax_msc_clump = plt.subplots(figsize=(12,12))
# msc_bulk, ax_msc_bulk = plt.subplots(figsize=(12,12))
# ax_msc_clump.plot(equib['D13CH3D'],equib['D12CH2D2'],'-k', label = 'Equilibrium', linewidth = 2.5, markersize = 15)
# for i in range(len(equib)):
#     if equib['p'].iloc[i]==1:
#         ax_msc_clump.scatter(equib['D13CH3D'].iloc[i], equib['D12CH2D2'].iloc[i],color='black',s=60)
# quick_plot_mix(ax_msc_clump,'D13CH3D','D12CH2D2','13CDse','DDse')
# ax_msc_clump.plot(D13CH3D_msc_clump,D12CH2D2_msc_clump, linewidth=3.0, color='blue', 
#                   label=r"Mixing curve: CH$_3$OH + H$_2$/CO$_2$(low-temp))")
# ax_msc_clump.set_xlabel('$\Delta^{13}$CH$_3$D (\u2030)', fontdict = font_labels)
# ax_msc_clump.set_ylabel('$\Delta^{12}$CH$_2$D$_2$ (\u2030)', fontdict = font_labels)
# ax_msc_clump.xaxis.set_major_locator(MultipleLocator(5))
# ax_msc_clump.yaxis.set_minor_locator(MultipleLocator(2))
# ax_msc_clump.xaxis.set_minor_locator(MultipleLocator(1))
# ax_msc_clump.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
# ax_msc_clump.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)
# ax_msc_clump.set_ylim([-60,13.5])
# ax_msc_clump.set_xlim([-10.5,5.1])

# quick_plot_mix(ax_msc_bulk,'d13C_CH4','dD_CH4','cse','dse')
# ax_msc_bulk.plot(d13C_msc,dD_msc, linewidth=3.0, color='blue', 
#                  label=r"Mixing curve: CH$_3$OH + H$_2$/CO$_2$(low-temp)")
# ax_msc_bulk.set_xlabel('$\delta^{13}$C (\u2030)', fontdict = font_labels)
# ax_msc_bulk.set_ylabel('$\delta$D (\u2030)', fontdict = font_labels)
# #ax_clump.yaxis.set_major_locator(MultipleLocator(5))
# ax_msc_bulk.xaxis.set_major_locator(MultipleLocator(20))
# ax_msc_bulk.yaxis.set_minor_locator(MultipleLocator(10))
# ax_msc_bulk.xaxis.set_minor_locator(MultipleLocator(5))
# ax_msc_bulk.tick_params(which='major',direction='out', top=True, right=True, length=8, width=2.5, labelsize=32)
# ax_msc_bulk.tick_params(which='minor',direction='out', top=True, right=True, length=4, width=2.0, labelsize=32)