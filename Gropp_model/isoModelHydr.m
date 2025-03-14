function [a2H_net,R2H,a2kr] = isoModelHydr(Tc,phiF,phiR,a2kf,R2_H2O,R2_H2)

% Calculate the isotopic fractionation between CH4 and H2O

%
% Output:
%         a2H_net = net isotopic fractionation between CH4 and H2O
%         R2H     = hydrogen isotopic ratios of metabolites
%         a13kr   = a vector of the reverse KFFs
%
% Input: 
%         Tc        = yemperature in degree C
%         phiF,phiR = gross forward and backward fluxes
%         fji       = a matrix of enzyme specific reversibility
%         a2kf      = a vector of forward KFFs
%         R2_H2O    = hydrogen isotope ratio of H2O
%         R2_H2     = hydrogen isotope ratio of H2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters and loop for ODE solver

R2H     = zeros(11,size(phiF,1));
a2H_net = zeros(1,size(phiF,1));

% EQUILIBRIUM Isotope Fractionation Factors, taken from Gropp et al., GCA,
% doi.org/10.1016/j.gca.2020.10.018. Reactions 1-7 are for the consecutive
% steps in the pathway, 8-10 are for the branching reactions and 11-12 are
% for EC recycling.

b_vals = [ 5.2327  -49.0324  166.1621 -245.9147  252.9625
           0.0797   -1.0541    3.6232   11.9088    0.6043
           0.8225   -9.7259   39.5466  -82.2848   23.5085
          -1.1204   13.8748  -61.0661   73.3352  -18.9078
           0.2246   -4.0886   28.8637  -38.4722   14.1934
           0.3687   -5.7787   36.1315  -51.4672   18.8685
           0.3490   -4.2149   18.4743  -16.3863    4.8821
          -0.1854    0.8537    9.8527  -22.0439    9.7272
           5.3661  -50.6782  171.3311 -248.3334  259.0032
          -0.4436    5.9171  -28.5163   15.6731   -4.1018
          -0.1139    0.5054    6.3783  -34.0472   13.9050
           6.1915  -63.5921  249.1940 -174.4374  228.2042
          -0.7896   10.2625  -44.1298 -144.5034   58.0480
          -3.4212   46.8356 -255.0463   38.4767   46.9858];

Tk = Tc + 273.15;
klna2eq = b_vals(:,1).*1e12./Tk.^4 + b_vals(:,2).*1e9./Tk.^3 + ...
          b_vals(:,3).*1e6./Tk.^2  + b_vals(:,4).*1e3./Tk + ...
          b_vals(:,5);
a2eq    = exp(klna2eq'./1000);
a2kr = a2eq.*a2kf;

for j = 1:size(phiF,1)
    phi_Rev = squeeze(phiR(j,:));
    phi_For = squeeze(phiF(j,:));
    C = [1 1 1 1 1 1 1 1 1 1].*1e-4';
    options = odeset('NonNegative',1:10,'AbsTol',1e-5,'RelTol',1e-2);
    %         options = odeset('NonNegative',1:10,'AbsTol',1e-15,'RelTol',1e-4);
    if isnan(phi_For(1))
        a2H_net(j) = NaN; 
    else
        [~,C] = ode15s(@frac_hydr_ODE,[0 1e10],C,options,phi_For,phi_Rev,a2kf,a2kr,R2_H2O,R2_H2);
        R2H(2,j)  = C(end,1);  % CHO-MFR
        R2H(3,j)  = C(end,2);  % CHO-H4MPT
        R2H(4,j)  = C(end,3);  % CH-H4MPT
        R2H(5,j)  = C(end,4);  % CH2-H4MPT (S)
        R2H(6,j)  = C(end,5);  % CH3-H4MPT
        R2H(7,j)  = C(end,6);  % CH3-S-CoM
        R2H(8,j)  = C(end,7);  % CH4
        R2H(9,j)  = C(end,8);  % F420H2
        R2H(10,j) = C(end,9);  % HS-CoB
        R2H(11,j) = C(end,10); % CH2-H4MPT (R)
        
        a2H_net(j) = R2H(8,j)./R2_H2O; % CH4 to H2O
    end
end
R2H(1,:) = R2_H2O;
ind = (phiF(:,1)-phiR(:,1) < 0);
a2H_net(ind) = NaN;
end

function ddt = frac_hydr_ODE(~,C,phi_For,phi_Rev,akf,akr,Rx,Rxh2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential equations describing the change of isotope ratio (R) of the
% metabolites in the pathway. This version is solvong assuming that dC/dt =
% 0 (where C is the concentration of a metabolite) bsed on the steady state
% solution of the metabolic model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddt = zeros(size(C));

% Assignment of the fluxes
jf = phi_For;
jr = phi_Rev;
jnet = jf(1) - jr(1);

Rw = C(1); Rv = C(2); Ru1 = C(3); RtS = C(4); 
RtR = C(10); Rs = C(5); Rr1 = C(6); Rp = C(7);
Ru2 = C(8); Rr2 = C(9);

xw    = jf(1).*Rx.*akf(1);      % h2o - cho-mfr
wx    = jr(1).*Rw.*akr(1);    % cho-mfr - h2o 
wv    = jf(2).*Rw.*akf(2);    % cho-mfr - cho-h4mpt
vw    = jr(2).*Rv.*akr(2);    % cho-h4mpt - cho-mfr
vu1   = jf(3).*Rv.*akf(3);    % cho-h4mpt - ch-h4mpt
u1v   = jr(3).*Ru1.*akr(3);    % ch-h4mpt - cho-h4mpt
u1tS  = (jf(4)+jf(10)).*Ru1.*akf(4);    % ch-h4mpt - ch2-h4mpt(S)
tSu1  = (jr(4)+jr(10)).*RtS.*akr(4);    % ch2-h4mpt(S) - ch-h4mpt
tSs   = jf(5).*RtS.*akf(5);    % ch2-h4mpt(S) - ch3-h4mpt
tRs   = jf(5).*RtR.*akf(6);   % ch2-h4mpt(R) - ch3-h4mpt
stS   = jr(5).*Rs.*akr(5);    % ch3-h4mpt - ch2-h4mpt(S)
stR   = jr(5).*Rs.*akr(6);    % ch3-h4mpt - ch2-h4mpt(R)
sr1   = 3.*jf(6).*Rs.*akf(7); % ch3-h4mpt - ch3-s-com
r1s   = 3.*jr(6).*Rr1.*akr(7); % ch3-s-com - ch3-h4mpt
r1p   = 3.*jf(7).*Rr1.*akf(8); % ch3-s-com - ch4
pr1   = 3.*jr(7).*Rp.*akr(8); % ch4 - ch3-s-com

xu2   = jf(8).*Rx.*akf(9);      % h2o - f420h2
u2x   = jr(8).*Ru2.*akr(9);    % f420h2 - h2o
u2s   = jf(5).*Ru2.*akf(11);   % f420h2 - ch3-h4mpt
su2   = jr(5).*Rs.*akr(11);   % ch3-h4mpt - f420h2
r2p   = jf(7).*Rr2.*akf(13);   % hscob - ch4
pr2   = jr(7).*Rp.*akr(13);   % ch4 - hscob
xr2   = jf(9).*Rx.*akf(12);     % h2o - hscob
r2x   = jr(9).*Rr2.*akr(12);   % hscob - h2o
u2tR  = jf(4).*Ru2.*akf(10);    % f420h2 - ch2-h4mpt(R)
xh2tR = jf(10).*Rxh2.*akf(14);  % h2 - ch2-h4mpt(R)
tRu2  = jr(4).*RtR.*akr(10);   % ch2-h4mpt(R) - f420h2
tRxh2 = jr(10).*RtR.*akr(14); % ch2-h4mpt(R) - h2
pout  = 4.*jnet.*Rp;

%% 
ddt(1)  = xw + vw - (wx + wv); % CHO-MFR
ddt(2)  = wv + u1v - (vw + vu1); % CHO-H4MPT
ddt(3)  = vu1 + tSu1 - (u1v + u1tS); % CH-H4MPT
ddt(4)  = u1tS + stS - (tSu1 + tSs); % CH2-H4MPT - S
ddt(10) = u2tR + xh2tR + stR - (tRu2 + tRxh2 + tRs); % CH2-H4MPT - R
ddt(5)  = tSs + tRs + u2s + r1s - (stS + stR + su2 + sr1); % CH3-H4MPT
ddt(6)  = sr1 + pr1 - (r1s + r1p); % CH3-S-CoM
ddt(7)  = r1p + r2p - (pr1 + pr2 + pout); % CH4
ddt(8)  = xu2 + tRu2 + su2 - (u2x + u2tR + u2s); % F420H2 with chiral center
ddt(9)  = xr2 + pr2 - (r2x + r2p); % HS-CoB
end





