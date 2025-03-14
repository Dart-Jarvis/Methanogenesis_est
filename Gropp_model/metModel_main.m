function [fji,Jcell,phiF,phiR,met] = metModel_main(Tc,h2,co2,ch4,args)
                                
                                
% DOCUMENTATION
% Output:
%         fji   - matrix of reversibility factors
%         Jcell - cell-specific rates in fmol/cell/d
%         phiF  - forward flux
%         phiR  - backward flux    
%         met   - concentrations of pathway metabolites
%
% Input: 
%         Tc   - temperature in degree celcius 
%         h2   - concentration row vector in M
%         co2  - concentration in M
%         ch4  - concentration in M
%         args - factorization for metabolic parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP DEFAULT MODEL PARAMETER VALUES
R = 8.31446; % Ideal gas constant.

% TUNABLE parameter values
Tk    = Tc + 273.15; % K
p     = length(h2);
h2_v  = h2;
% co2_v = repmat(co2,1,length(h2));
% ch4_v = repmat(ch4,1,length(h2));
co2_v = co2;
ch4_v = ch4;

%% CALCULATE REVERSIBILITY FACTORS

metModel_param % load kinetic parameters from desginated script

[fji,J,phiF,phiR,met] = ...
    metModel_ODEsol(h2_v,co2_v,ch4_v,Tk,R,p,par,args);

Jcell = J.*args(1).*86400.*1e15; % Convert mol/L/s to fmol/cell/day


