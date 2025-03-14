%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses fluxes calculated by 'metModel_standalone.m' to
% predict Delta(13CH3D). The mass balance derivation follows largely the
% hydrogen mass balance, and the fractionation factors are expressed as
% gamma, the deviation from stochasticity.

function [D12CH2D2,RDD] = isoModelClumped_DD(Tc,R2_H2O,R2_H2,phiF,phiR,a2kinf)
%        
% Output
%        D12CH2D2 - array of methane clumped isotopologue distribution
%        RDD      - compund specific isotopologue ratios
% 
% Input
%        Tc             - Temp in degree celsius
%        R2_H2O, R2_H2  - isotopic ratios of H2O and H2
%        f              - matrix of reversibilty values
%        phiF,phiR      - matrices of forward and backward fluxes
%        a2kinf,a13kinf - optimized H and C KFFs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve hydrogen isotope fractionation model
[~,Rh] = isoModelHydr(Tc,phiF,phiR,a2kinf,R2_H2O,R2_H2);

% Derive compound specific isotope ratios
R_D_u  = Rh(4,:);  % CH-H4MPT
R_D_tS = Rh(5,:);  % CH2-H4MPT (S)
R_D_tR = Rh(11,:); % CH2-H4MPT (R)
R_D_r  = Rh(7,:);  % CH3-S-CoM
R_D_p  = Rh(8,:);  % CH4
R_D_i  = Rh(9,:);  % F420
R_D_j  = Rh(10,:); % HS-CoB

p = length(phiF(:,1));

% Define minimal primary (p) and secondary (s) gamma. gamma < 0 is normal. 
gammap = 0.994;
gammas = 1;

gammaf = zeros(1,8);
gammaf(1:3) = gammas(1) + (1 - gammas(1))*rand(1,3);
gammaf(4:8) = gammap(1) + (1 - gammap(1))*rand(1,5);

Tk = Tc + 273.15;

% EQUILIBRIUM Isotope Fractionation Factors, taken from Gropp et
% al., GCA, doi.org/10.1016/j.gca.2020.10.018
b = [ 9.96831E-18 -4.26734E-14  7.60881E-11 -7.28978E-08  4.01300E-05 -1.23247E-02  2.70904E+00 % Mer s
      1.55513E-18 -6.74282E-15  1.23156E-11 -1.23331E-08  7.35249E-06 -2.59468E-03  1.45398E+00 % Mtr s
      1.15066E-17 -4.88781E-14  8.61571E-11 -8.10885E-08  4.33637E-05 -1.26788E-02  2.61499E+00 % Mcr s
     -5.35819e-19  2.68382E-15 -5.95919E-12  7.66964E-09 -6.15692E-06  2.97095E-03  3.20812E-01 % Mtd p
      6.58835E-18 -2.85197E-14  5.05845E-11 -4.60941E-08  2.11152E-05 -2.90011E-03  1.06681E-01 % Hmd p
      5.13100E-18 -2.19709E-14  3.91503E-11 -3.74090E-08  2.04360E-05 -6.14045E-03  1.79016E+00 % Mer pS
      5.02110E-18 -2.15090E-14  3.83592E-11 -3.67124E-08  2.01145E-05 -6.07080E-03  1.78194E+00 % Mer pR
     -2.38222E-19  1.40757E-15 -3.65466E-12  5.44542E-09 -5.09692E-06  3.04268E-03 -2.31867E-02 % Mcr p
     ];
a22eq = permute(b(:,1).*Tk.^6 + b(:,2).*Tk.^5 + b(:,3).*Tk.^4 + b(:,4).*Tk.^3 + b(:,5).*Tk.^2 + b(:,6).*Tk + b(:,7), [2,1]); 

% Calculate KFFs from gamma and optimized C and H KFFs
a22kf(1) = gammaf(1).*a2kinf(5).*a2kinf(6);   % Mer, s (tS->s)
a22kf(2) = gammaf(2).*a2kinf(7).^2;           % Mtr, s
a22kf(3) = gammaf(3).*a2kinf(8).^2;           % Mcr, s
a22kf(4) = gammaf(4).*a2kinf(4).*a2kinf(10);  % Mtd, p
a22kf(5) = gammaf(5).*a2kinf(4).*a2kinf(14);  % Hmd, p
a22kf(6) = gammaf(6).*a2kinf(5).*a2kinf(11);  % Mer, pS
a22kf(7) = gammaf(7).*a2kinf(6).*a2kinf(11);  % Mer, pR
a22kf(8) = gammaf(8).*a2kinf(8).*a2kinf(13);  % Mcr, p

a22kr = a22eq.*a22kf;

% Calculate Delta12CH2D2 
D12CH2D2 = zeros(1,p);
RDD = zeros(4,p);
for i = 1:p
    phif_i = phiF(i,:);
    phir_i = phiR(i,:);
    if phif_i(1)-phir_i(1)<0 || isnan(phif_i(1))
        D12CH2D2(i) = nan;
    elseif phif_i(1)-phir_i(1)>=0
        C = [1 1 1 1].*1e-8';
        options = odeset('NonNegative',1:4,'AbsTol',1e-5,'RelTol',1e-2);
        [~,C] = ode15s(@frac_hydr_ODE,[0 1e10],C,options,phif_i,...
            phir_i,R2_H2,R_D_i(i),R_D_j(i),R_D_u(i),R_D_r(i),...
            R_D_tS(i),R_D_tR(i),a22kf,a22kr);
        RDD(1,i) = C(end,1);       % 12CD2-H4MPT
        RDD(2,i) = C(end,2); % 12CHD2-H4MPT
        RDD(3,i) = C(end,3); % 12CHD2-SCoM
        RDD(4,i) = C(end,4); % 12CH2D2  

        D12CH2D2(i) = (1000.*(RDD(4,i)./(R_D_p(i).^2) - 1));
    end
end
end

function ddt = frac_hydr_ODE(~,C,phif,phir,RxH2,R_D_i,R_D_j,R_D_u,R_D_r,...
                             R_D_tS,R_D_tR,akinf,akinr)
% ODE solver function.

ddt = zeros(size(C));

phinet = phif(7)-phir(7);

R_DD_t = C(1); 
R_DD_s = C(2);
R_DD_r = C(3);
R_DD_p = C(4);

u_t   = R_D_u*(phif(4)*R_D_i*akinf(4)+phif(10)*RxH2*akinf(5));
t_u   = R_DD_t*(phir(4)*akinr(4)+phir(10)*akinr(5));
t_s   = phif(5)*R_DD_t*akinf(1);
s_t   = phir(5)*R_DD_s*akinr(1);
ti_s  = phif(5)*R_D_i*(R_D_tS*akinf(6) + R_D_tR*akinf(7));
s_ti  = phir(5)*R_DD_s*(akinr(6)+akinr(7));
s_r   = 3*phif(6)*R_DD_s*akinf(2);
r_s   = 3*phir(6)*R_DD_r*akinr(2);
rj_p  = 3*phif(7)*R_D_j*R_D_r*akinf(8);
p_rj  = 3*phir(7)*R_DD_p*akinr(8);
r_p   = 3*phif(7)*R_DD_r*akinf(3);
p_r   = 3*phir(7)*R_DD_p*akinr(3);
p_out = 6*R_DD_p*phinet;

ddt(1) = u_t + s_t - (t_u + t_s);
ddt(2) = t_s + ti_s + r_s - (s_t + s_ti + s_r);
ddt(3) = s_r + p_r - (r_s + r_p);
ddt(4) = rj_p + r_p - (p_rj + p_r + p_out);
end
