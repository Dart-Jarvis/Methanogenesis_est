%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses fluxes calculated by 'metModel_standalone.m' to
% predict Delta(13CH3D). The mass balance derivation follows largely the
% hydrogen mass balance, and the fractionation factors are expressed as
% gamma, the deviation from stochasticity. In this derivation, KFF > 1
% causes clumping (inverse KFF), and KFF < 1 causes anti-clumping (normal
% KFF).

function [D13CH3D,R13D] = ...
    isoModelClumped_13D(Tc,R2_H2O,R2_H2,R13_CO2,f,phiF,phiR,a2kinf,a13kinf)
%        
% Output
%        D13CH3D - array of methane clumped isotopologue distribution
%        R13D    - compund specific isotopologue ratios
% 
% Input
%        Tc             - Temp in degree celsius
%        R2_H2O, R2_H2, R13C_CO2 - isotopic ratios of H2O, H2 and CO2
%        f              - matrix of reversibilty values
%        phiF,phiR      - matrices of forward and backward fluxes
%        a2kinf,a13kinf - optimized H and C KFFs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve hydrogen and carbon isotope fractionation models
[~,Rh]       = isoModelHydr(Tc,phiF,phiR,a2kinf,R2_H2O,R2_H2);
[aC,a13c_dat] = isoModelCarb(Tc,f',1./a13kinf);

% Derive compound specific isotope ratios
Ru2_D   = Rh(9,:);              % F420
Rr2_D   = Rh(10,:);             % HS-CoB
Rp_D    = Rh(8,:);              % CH4
Rp_13   = R13_CO2./aC;          % CH4
Ru1_13  = Rp_13.*a13c_dat(4,:); % CH-H4MPT
Rt_13   = Rp_13.*a13c_dat(3,:); % CH2-H4MPT
Rr1_13  = Rp_13.*a13c_dat(1,:); % CH3-S-CoM
RxCin   = Rp_13.*a13c_dat(7,:); % intracellular CO2

p = length(phiF(:,1));

% Define minimal primary (p) and secondary (s) gamma range.
gammap = 0.997;
gammas = 1;

gammaf = ones(1,12);
gammaf([2 3 4 5 6 7]) = gammas(1) + (1 - gammas(1))*rand(1,6);
gammaf([1 8:12])      = gammap(1) + (1 - gammap(1))*rand(1,6);

% EQUILIBRIUM Isotope Fractionation Factors, taken from Gropp et
% al., GCA, doi.org/10.1016/j.gca.2020.10.018
b_vals = [ 5.04286 -46.55646  154.26902 -220.44125 241.20163
           0.11987  -1.46945    4.88531   12.53630   0.74283
           0.84641  -9.99803   40.49548  -84.72305  24.53775
          -1.04439  13.08131  -58.23150   75.42109 -20.49841
           0.27840  -4.97688   34.63646  -46.41678  17.36500
           0.38561  -6.30570   40.61156  -57.46741  20.95643
           0.36154  -4.50907   20.64966  -15.06848   3.60230
          -0.17407   0.60937   12.69865  -29.73599  13.05395
          -0.39285   5.37175  -26.59094   17.44264  -4.97376
          -3.35961  46.21746 -253.01593   40.37976  45.84330
          -0.07497  -0.18489   11.09870  -41.64583  17.35250
          -0.83925  10.64855  -43.85178 -149.76316  60.55202];

Tk = Tc + 273.15;

klna132eq = b_vals(:,1).*1e12./Tk.^4 + b_vals(:,2).*1e9./Tk.^3 + ...
          b_vals(:,3).*1e6./Tk.^2  + b_vals(:,4).*1e3./Tk + ...
          b_vals(:,5);
a132eq    = exp(klna132eq'./1000);

% KFFs from gamma and optimized C and H KFFs
a132kf(:,1)  = gammaf(1).*a2kinf(1).*a13kinf(1);   % Fmd, p
a132kf(:,2)  = gammaf(2).*a2kinf(2).*a13kinf(2);   % Ftr, s
a132kf(:,3)  = gammaf(3).*a2kinf(3).*a13kinf(3);   % Mch, s
a132kf(:,4)  = gammaf(4).*a2kinf(4).*a13kinf(4);   % Mtd, s (u1->tS)
a132kf(:,5)  = gammaf(5).*a2kinf(5).*a13kinf(5);   % Mer, s (tS->s)
a132kf(:,6)  = gammaf(6).*a2kinf(6).*a13kinf(5);   % Mer, s (tR->s)
a132kf(:,7)  = gammaf(7).*a2kinf(7).*a13kinf(6);   % Mtr, s
a132kf(:,8)  = gammaf(8).*a2kinf(8).*a13kinf(7);   % Mcr, s
a132kf(:,9)  = gammaf(9).*a2kinf(10).*a13kinf(4);  % Mtd, p
a132kf(:,10) = gammaf(10).*a2kinf(14).*a13kinf(4); % Hmd, p
a132kf(:,11) = gammaf(11).*a2kinf(11).*a13kinf(5); % Mer, p
a132kf(:,12) = gammaf(12).*a2kinf(13).*a13kinf(7); % Mcr, p

a132kr = a132eq.*a132kf;

% Calculate Delta13CH3D 
D13CH3D = zeros(1,p);
R13D = zeros(8,p);
for i = 1:p
    phif_i = phiF(i,:);
    phir_i = phiR(i,:);
    if phif_i(1)-phir_i(1)<0 || isnan(phif_i(1))
        D13CH3D(i) = nan;
    elseif phif_i(1)-phir_i(1)>=0
        C = [1 1 1 1 1 1 1 1].*1e-6';
        options = odeset('NonNegative',1:8,'AbsTol',1e-5,'RelTol',1e-2);
        [~,C] = ode15s(@frac_hydr_ODE,[0 1e10],C,options,R2_H2O,RxCin(i),phif_i,...
            phir_i,Ru2_D(i),Rr2_D(i),Ru1_13(i),Rt_13(i),Rr1_13(i),...
            R2_H2,a132kf,a132kr);
        R13D(1,i) = C(end,1);       % cho-mfr
        R13D(2,i) = C(end,2);       % cho-h4mpt
        R13D(3,i) = C(end,3);       % ch-h4mpt
        R13D(4,i) = C(end,4);       % ch2-h4mpt (S)
        R13D(5,i) = C(end,5);       % ch2-h4mpt (R)
        R13D(6,i) = C(end,6); % ch3-h4mpt
        R13D(7,i) = C(end,7); % ch3-s-com
        R13D(8,i) = C(end,8); % ch4
        D13CH3D(i) = (1000.*(R13D(8,i)./(Rp_13(i).*Rp_D(i)) - 1));
    end
end
end

function ddt = frac_hydr_ODE(t,C,RxH,RxCin,phif,phir,Ru2_D,...
    Rr2_D,Ru1_13,Rt_13,Rr1_13,Rxh2,akinf,akinr)
% ODE solver function.

ddt = zeros(size(C));

phinet = phif(1)-phir(1);
R_13C_w  = C(1);
R_13C_v  = C(2);
R_13C_u1 = C(3);
R_13C_tS = C(4); 
R_13C_tR = C(5);
R_13C_s  = C(6);
R_13C_r1 = C(7);
R_13C_p  = C(8);

xw    = phif(1)*RxH*RxCin*akinf(1);
wx    = phir(1)*R_13C_w*akinr(1);
wv    = phif(2)*R_13C_w*akinf(2);
vw    = phir(2)*R_13C_v*akinr(2);
vu1   = phif(3)*R_13C_v*akinf(3);
u1v   = phir(3)*R_13C_u1*akinr(3);
u1tS  = (phif(4)+phif(10))*R_13C_u1*akinf(4);
tSu1  = (phir(4)+phir(10))*R_13C_tS*akinr(4);
tSs   = phif(5)*R_13C_tS*akinf(5);
tRs   = phif(5)*R_13C_tR*akinf(6);
sr1   = 3*phif(6)*R_13C_s*akinf(7);
r1s   = 3*phir(6)*R_13C_r1*akinr(7);
r1p   = 3*phif(7)*R_13C_r1*akinf(8);
utR   = Ru1_13*(phif(4)*Ru2_D*akinf(9)+phif(10)*Rxh2*akinf(10));
tRu2  = R_13C_tR*(phir(4)*akinr(9)+phir(10)*akinr(10));
u2s   = phif(5)*Rt_13*Ru2_D*akinf(11);
r2p   = phif(7)*Rr1_13*Rr2_D*akinf(12);
pout  = 4*phinet*R_13C_p;
stS   = phir(5)*R_13C_s*akinr(5);
stR   = phir(5)*R_13C_s*akinr(6);
pr1   = 3*phir(7)*R_13C_p*akinr(8);
su2   = phir(5)*R_13C_s*akinr(11);
pr2   = phir(7)*R_13C_p*akinr(12);

ddt(1) = xw + vw - (wx + wv); % A
ddt(2) = wv + u1v - (vw + vu1); % B
ddt(3) = vu1 + tSu1 - (u1v + u1tS); % C
ddt(4) = u1tS + stS - (tSu1 + tSs); % D
ddt(5) = utR + stR - (tRu2 + tRs); % E
ddt(6) = tSs + tRs + u2s + r1s - (stS + stR + su2 + sr1); % F
ddt(7) = sr1 + pr1 - (r1s + r1p); % G
ddt(8) = r1p + r2p - (pr1 + pr2 + pout); % H

end
