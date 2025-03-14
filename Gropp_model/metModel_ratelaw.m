% Define kinetic parameters
%% ATPase
yatp = 0.2;

%% Calculate standard transformed dG for ATP formation for a given Tk
dH    = 20.5e3;
VH    = exp((-dH./R).*(1./Tk-1./298.15));
K1    = exp(-43.5e3/R./298.15);
dGatp = yatp*-R.*Tk.*log(VH.*K1);

%% Load metabolic parameters from 'metModel_param.m'. 
dG  = par(1,:);
v   = par(2,:);
k1  = par(3,:);
k2  = par(4,:);
k3  = par(5,:);
k4  = par(6,:);
k5  = par(7,:);
k6  = par(8,:);
k1a = par(9,:);

%% Full rate laws
q     = length(h2i);
Fk    = zeros(1,q,10);
Ft    = zeros(1,q,10);
J_net = zeros(1,q,10);

% Fmd (calculated for reverse direction)
Fk(:,:,1) = (fdox.*formmfr./(k1(1).*k2(1)))./(1+(fdred.*mfr.*co2in./(k1a(1)))+(fdox.*formmfr./(k1(1).*k2(1))));
Ft(:,:,1) = 1 - (fdred.*mfr.*co2in./(fdox.*formmfr)).*exp(-dG(1)./(R.*Tk));
% Ftr
Fk(:,:,2) = ((formmfr.*h4mpt./(k1(2).*k2(2)))./(1 + formmfr.*h4mpt./(k1(2).*k2(2)) + formh4mpt.*mfr./(k1a(2))));
Ft(:,:,2) = 1 - (formh4mpt.*mfr./(formmfr.*h4mpt)).*exp(dG(2)./(R.*Tk));
% Mch (calculated for reverse direction)
Fk(:,:,3) = ((menylh4mpt./k1(3))./(1 + menylh4mpt./k1(3) + formh4mpt./k1a(3)));
Ft(:,:,3) = 1 - (formh4mpt./menylh4mpt).*exp(-dG(3)./(R.*Tk));
% Mtd (calculated for reverse direction)
Fk(:,:,4) = ((mleneh4mpt.*f420./(k1(4).*k2(4)))./(1+mleneh4mpt.*f420./(k1(4).*k2(4))+menylh4mpt.*f420h2./(k1a(4))));
Ft(:,:,4) = 1 - (menylh4mpt.*f420h2./(mleneh4mpt.*f420)).*exp(-dG(4)./(R.*Tk));
% Hmd (calculated for reverse direction)
Fk(:,:,10) = (mleneh4mpt./k1a(10))./...
    (1 + (mleneh4mpt./k1a(10)) + (menylh4mpt.*h2i./(k1(10).*k2(10))));
Ft(:,:,10) = 1 - (menylh4mpt.*h2i./mleneh4mpt).*exp(-dG(10)./(R.*Tk));
% Mer
Fk(:,:,5) = ((mleneh4mpt.*f420h2./(k1(5).*k2(5)))./(1+mleneh4mpt.*f420h2./(k1(5).*k2(5)) + mh4mpt.*f420/(k1a(5))));
Ft(:,:,5) = 1 - (mh4mpt.*f420./(mleneh4mpt.*f420h2)).*exp(dG(5)./(R.*Tk));
% Mtr
Fk(:,:,6) = ((mh4mpt.*com./(k1(6).*k2(6)))./(1 + mh4mpt.*com./(k1(6).*k2(6)) + mcom.*h4mpt./(k1a(6))));
Ft(:,:,6) = 1 - (mcom.*h4mpt./mh4mpt./com).*exp((dG(6)+dGatp)./(R.*Tk));
% Mcr
Fk(:,:,7) = ((mcom.*cob./(k1(7).*k2(7)))./(1+mcom.*cob./(k1(7).*k2(7))+ch4i.*hsfd./(k1a(7))));
Ft(:,:,7) = 1 - (ch4i.*hsfd./(mcom.*cob)).*exp(dG(7)./(R.*Tk));
% Frh
Fk(:,:,8) = ((f420.*h2i./(k1(8).*k2(8)))./(1 + f420.*h2i./(k1(8).*k2(8)) + f420h2./k1a(8)));
Ft(:,:,8) = 1 - (f420h2./(f420.*h2i)).*exp(dG(8)./(R.*Tk));
% Hdr
Fk(:,:,9) = ((h2i.^2.*hsfd.*fdox./(k1(9).^2.*k2(9).*k3(9)))./...
    (1 + h2i.^2.*hsfd.*fdox./(k1(9).^2.*k2(9).*k3(9)) + com.*cob.*fdred./(k4(9).*k5(9).*k6(9))));
Ft(:,:,9) = 1 - (com.*cob.*fdred./(h2i.^2.*hsfd.*fdox)).*exp(dG(9)./(R.*Tk));

% CO2 diffusion
Vcell     = args(1);
dX        = 0.5e-9;                       % Membrane thickness, m. (~0.5 nm)
D         = 2.9e-9;                       % Diffusivity constant, Missner 2008, m^2 s^-1
rad       = (Vcell*1e-3*3/4/pi)^(1/3);    % Cell radius, m
A         = 4*pi*rad^2;                   % Cell surface area (m^2)  
Jco2_cell = D.*A.*(co2i-co2in)./1000./dX; % mol/cell/s

J_net(:,:,1)  = -v(1) .*Fk(:,:,1) .*Ft(:,:,1);
J_net(:,:,2)  =  v(2) .*Fk(:,:,2) .*Ft(:,:,2);
J_net(:,:,3)  = -v(3) .*Fk(:,:,3) .*Ft(:,:,3);
J_net(:,:,4)  = -v(4) .*Fk(:,:,4) .*Ft(:,:,4);
J_net(:,:,5)  =  v(5) .*Fk(:,:,5) .*Ft(:,:,5);
J_net(:,:,6)  =  v(6) .*Fk(:,:,6) .*Ft(:,:,6);
J_net(:,:,7)  =  v(7) .*Fk(:,:,7) .*Ft(:,:,7);
J_net(:,:,8)  =  v(8) .*Fk(:,:,8) .*Ft(:,:,8);
J_net(:,:,9)  =  v(9) .*Fk(:,:,9) .*Ft(:,:,9);
J_net(:,:,10) = -v(10).*Fk(:,:,10).*Ft(:,:,10);
J_net(:,:,11) = Jco2_cell./Vcell;

% Simulate model without Hmd activity
if args(2) == 0
    J_net(:,:,10) = 0;
end

