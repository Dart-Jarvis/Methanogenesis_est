% A model to run the differential equations
function [fji,J,phi_for,phi_rev,met] = ...
                    metModel_ODEsol(h2mat,co2mat,ch4mat,Tk,R,p,par,args)

met = zeros(length(co2mat(:,1)),p,16);
for i = 1:p
    for j = 1:length(co2mat(:,1))
        h2i  = h2mat(j,i);
        co2i = co2mat(j,i);
        ch4i = ch4mat(j,i);
        % Initial concentrations of metabolites for ODE solver
        C = [1.0e-7 1.8e-3 5.0e-3 1e-7   1.2e-3 1.5e-5 7.5e-5 2.5e-2...
             5.8e-4 1e-7   1.5e-2 6.0e-3 1.0e-7 6.0e-3 1.0e-7 1e-2]'; % [M]
        n = length(C);
        
        TIMEVECTOR = [0 1e8];

        options = odeset('NonNegative',1:n,'AbsTol',1e-20,'RelTol',1e-4);
        [t,C] = ode15s(@metModel_difeqns,TIMEVECTOR,C,options,Tk,R,par,h2i,co2i,ch4i,args);
%         keyboard
        met(j,i,:) = C(end,1:n);
    end
end

[J,fji,phi_for,phi_rev] = metModel_FluxRev(met,h2mat,co2mat,ch4mat,Tk,R,par,args); 

% Remove simulations with no net methanogenesis flux (due to positive dGr)
neg = isnan(J(:,:,1));
[negR,negC] = find(neg == 1);
for ineg = 1:length(negR)
    met(negR(ineg),negC(ineg),:) = NaN;
end

end

function ddt = metModel_difeqns(t,C,Tk,R,par,h2i,co2i,ch4i,args)

% This function generates differential equations for the change of metabolite
% concentrations through time to solve in the main function 
% metabol_netFlux.
%
% Input:           
%                  C - a matrix of initial concentrations of metabolites, 
%                      listed in metabol_param.
%
% Output:          ddt - a matrix containing the solutions of the
%                        differential equations for 15 metabolites.
%                        This list goes to solver

ddt = zeros(size(C));

formmfr    = C(1);
mfr        = C(2);
fdox       = C(3);
fdred      = C(4);
h4mpt      = C(5);
formh4mpt  = C(6);
menylh4mpt = C(7);
mleneh4mpt = C(8);
f420       = C(9);
f420h2     = C(10);
mh4mpt     = C(11);
com        = C(12);
mcom       = C(13);
cob        = C(14);
hsfd       = C(15);
co2in      = C(16);

metModel_ratelaw

ddt(1)  = J_net(:,:,1) - J_net(:,:,2);
ddt(2)  = -ddt(1);
ddt(3)  = J_net(:,:,1) - J_net(:,:,9);
ddt(4)  = -ddt(3);
ddt(5)  = J_net(:,:,6) - J_net(:,:,2);
ddt(6)  = J_net(:,:,2) - J_net(:,:,3);
ddt(7)  = J_net(:,:,3) - J_net(:,:,4) - J_net(:,:,10);
ddt(8)  = J_net(:,:,4) + J_net(:,:,10) - J_net(:,:,5);
ddt(9)  = J_net(:,:,4) + J_net(:,:,5) - J_net(:,:,8);
ddt(10) = -ddt(9);
ddt(11) = J_net(:,:,5) - J_net(:,:,6);
ddt(12) = J_net(:,:,9) - J_net(:,:,6);
ddt(13) = J_net(:,:,6) - J_net(:,:,7);
ddt(14) = J_net(:,:,9) - J_net(:,:,7);
ddt(15) = -ddt(14);
ddt(16) = J_net(:,:,11) - J_net(:,:,1);
end

function [J,f,phi_for,phi_rev] = ...
                metModel_FluxRev(met,h2mat,co2mat,ch4mat,Tk,R,par,args)

% Description: this function calculates the net flux from the metabolites
% concentrations provided by the ODE solver, and from that the
% reversibility of each reaction.

h2i  = h2mat;
co2i = co2mat;
ch4i = ch4mat;

% Find the net fluxes according to the solution:
formmfr    = met(:,:,1);
mfr        = met(:,:,2);
fdox       = met(:,:,3);
fdred      = met(:,:,4);
h4mpt      = met(:,:,5);
formh4mpt  = met(:,:,6);
menylh4mpt = met(:,:,7);
mleneh4mpt = met(:,:,8);
f420       = met(:,:,9);
f420h2     = met(:,:,10);
mh4mpt     = met(:,:,11);
com        = met(:,:,12);
mcom       = met(:,:,13);
cob        = met(:,:,14);
hsfd       = met(:,:,15);
co2in      = met(:,:,16);

metModel_ratelaw

% REVERSIBILITY
f = 1 - Ft;
% f for reactions that go in reverse is the inverse of Ft
f(:,:,1)  = 1./f(:,:,1);
f(:,:,3)  = 1./f(:,:,3);
f(:,:,4)  = 1./f(:,:,4);
f(:,:,10) = 1./f(:,:,10);
f(:,:,11) = co2in./co2i;

% Net Flux
J = J_net;
% Remove simulations with no net methanogenesis flux (due to positive dGr)
neg = find(J(:,:,2) <= 0);
J(:,neg,:) = NaN;
f(:,neg,:) = NaN;

% Forward and reverse fluxes
phi_rev = f.*J./(1-f);
phi_for = phi_rev./f;
f(:,:,12) = (phi_rev(:,:,4)+phi_rev(:,:,10))./(phi_for(:,:,4)+phi_for(:,:,10));
end

