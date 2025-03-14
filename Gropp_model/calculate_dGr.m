function dgr = calculate_dGr(Tc,h2,co2,ch4)

% This function takes as imput temperature, CO2, H2 and CH4 concentrations
% in M, and produces a vector of transformed Gibbs free energies (dGr')
% 
% Output:
%       dGr - transformed Gibbs free energy (kJ/mol)
% Input:
%       Tc  - temperature in degree Celsius
%       H2  - H2 concentration (M)
%       CO2 - CO2 concentration (M)
%       CH4 - CH4 concentration (M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the dG'0 of methanogenensis at specified temperature using
% Van't Hoffs equation. Values from Alberty, 2003.
T1  = 298.15;
T2  = Tc + 273.15;
Ha  = [-90.68,-286.65,-413.8,-5.02]; 
R   = 8.31446;
dG  = -193.13e3;
dH  = ((Ha(1) + 2*Ha(2)) - (Ha(3) + 4*Ha(4))).*1000;
vh  = (-dH./R).*((1./T2)-(1./(T1)));
dg0 = -R.*T2.*log(exp((-dG./R./T1)+vh));

dgr = ((dg0 + R.*(Tc+273.15).*log(ch4./(co2.*(h2.^4))))./1000);

end