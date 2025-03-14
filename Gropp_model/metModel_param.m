% Metabolic parameters for the model

% Standard transformed Gibbs free energies
par(1,:) = [(10.3e3) -3.5e3 (-4.2e3) (2.0e3) -1.7e3...
            0        0      -28.0e3  -90.8e3 (-26.0e3)];

% Calculate the dG'0 of methanogenensis at specified temperature using
% Van't Hoffs equation. Values from Alberty, 2003.
T1  = 298.15;
T2  = Tk;
Ha  = [-90.68,-286.65,-413.8,-5.02]; 
dG  = -193.13e3;
dH  = ((Ha(1) + 2*Ha(2)) - (Ha(3) + 4*Ha(4))).*1000;
vh  = (-dH./R).*((1./T2)-(1./(T1)));
dg0 = -R.*T2.*log(exp((-dG./R./T1)+vh));

% Due to uncertainty in dG'0 values of Mtr and Mcr, their dG'0 is set a the
% residual from the total dG'0 of the pathway.
totdg    = sum([par(1,1:9) par(1,8)]);
dg0tot   = dg0-totdg;
Rdg0     = 0.4;             % Allocation ratio of dG0 between Mtr and Mcr
par(1,6) = dg0tot.*(1-Rdg0); % dG0 Mtr
par(1,7) = dg0tot.*Rdg0;     % dG0 Mcr

% Maximal rate capacities of enzymes (V+)
par(2,:)   = [6.086e-3,2.853e-2,4.032e-3,1.292e-2,6.086e-3,...
              2.380e-3,1.041e-2,1.415e-2,1.902e-2,3.041e-3];
% KM values          
par(3,:)   = [30e-6,50e-6,167e-6,33e-6,  3e-6,135e-6,821e-6,12e-6, 30e-6,150e-6];
par(4,:)   = [10e-6,60e-6,   NaN,65e-6,300e-6,277e-6,204e-6,36e-6,145e-6, 50e-6];
par(5:8,9) = [10e-6,200e-6,200e-6,75e-6];
par(9,:)   = [6800e-6^3,4.7e-6^2,148.1e-6,16e-6*50e-6,24.4e-6*40e-6,...
              558.7e-6*98.3e-6,109.8e-6^2,11.7e-6,NaN,40e-6];

if args(3) == 0 % For Mcr I
    % Change Vmax of particular enzyme
    par(2,7) = 7.6e-2;     % Total Mcr activity from Pennings 2000
    par(2,8) = par(2,8)*2; % Pennings 2000
    par(2,9) = par(2,9)*3; % Pennings 2000
    % Change Km of particular enzyme - Mcr I
    par(3,7) = 280e-6;
    par(4,7) = 75e-6;
end

% Application of the Uviv/vit parameter
par(2,:) = 10.6.*par(2,:); % Total Vmax, scaling net flux along y-axis

% Temperature dependence of enzyme activities
vplus_scaling = args(4);
par(2,:) = par(2,:).*vplus_scaling.^((Tc-60)./10);




