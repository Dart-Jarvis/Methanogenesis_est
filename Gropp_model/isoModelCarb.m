function [a13C_net,a13Cdat,a13kr] = isoModelCarb(Tc,fji,a13kf)

% Calculate the isotopic fractionation between CO2 and CH4

% DOCUMENTATION
% Output:
%         a13C_net = net isotopic fractionation between CO2 and CH4
%         a13C_dat = net isotopic fractionation between intermediates and CH4
%         a13kr    = a vector of the reverse kinetic fractionation factors
%                    (KFFS) calculated from the forward KFFs and 
%                    equilibrium fractionation factors.
%
% Input: 
%         Tc    = temperature in degree C
%         fji   = a matrix of enzyme specific reversibility
%         a13kr = a vector of forward KFFs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EQUILIBRIUM Isotope Fractionation Factors, taken from Gropp et
% al., GCA, doi.org/10.1016/j.gca.2020.10.018
b_vals = [-0.03507  0.89189 -5.78478 16.97954 -7.27610;
           0.02543 -0.28251  0.78966  1.26127 -0.16561;
           0.02412 -0.23329  0.55401 -1.22821  0.26467;
           0.11085 -1.18789  4.54994 -0.98362  0.14324;
           0.02240 -0.57894  4.73613 -6.21676  2.17747;
           0.06685 -0.84498  4.16633 -1.80360  0.51959;
          -0.04853  0.32124  0.93815 -4.85941  1.76778];

Tk = Tc+273.15;
klna13eq = b_vals(:,1).*1e12./Tk.^4 + b_vals(:,2).*1e9./Tk.^3 + ...
        b_vals(:,3).*1e6./Tk.^2  + b_vals(:,4).*1e3./Tk + ...
        b_vals(:,5);
a13eq = exp(klna13eq'./1000);

a13eq(8) = 1; % Diffusivity EFF
a13kr    = a13eq.*(1./a13kf);

% Calculate NET CARBON ISOTOPIC FRACTIONATION
% Produce a matrix of net carbon isotopic fractionation for given conditions.
f      = fji;
f(8,:) = fji(11,:); % Diffusion
f(4,:) = fji(12,:); % Combine Mtd and Hmd

pos = [7 6 5 4 3 2 1 8];
ar_size = size(f);
for k = 1:ar_size(2)
    a13C = 1;
    for n = 1:8 
        a13C = f(pos(n),k).*(a13eq(pos(n)).*a13C - a13kf(pos(n))) + a13kf(pos(n));
        a13Cdat(n,k) = a13C;
    end
    a13C_net(k) = a13C;
end

