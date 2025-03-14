%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A code to run a bio-isotopic model of hydrogenotrophic methanogenesis.
% Last edit: Dec 7th, 2021 (Jonathan Gropp)

clear

% DEFINE ENVIRONMENTAL CONDITIONS
sims = 30;
H2   = logspace(-8,-2,sims);           
CO2  = repmat(1e-2,1,sims);
CH4  = repmat(1e-5,1,sims);
Tc   = 35;
Tk   = Tc + 273.15;
dGr  = calculate_dGr(Tc,H2,CO2,CH4);

% DEFINE TUNABLE PARAMETERS
cell_vol   = 2e-15; % L
Hmd_act    = 1;     % 1 to include Hmd in metabolic simulation, 0 to ignore.
Mcr_isoenz = 1;     % 1 for McrII (growth with replete H2), 0 for McrI (growth under H2 limitation)
Q10Scale   = 1;     % An Arrhenius type scaling parameter for enzyme activities, we recommend using values between 1-2.

% DEFINE ISOTOPIC COMPOSITIONS FOR CO2 AND H2O
d13CCO2   = -36.01; % permil
dDH2O     = -50; % permil

args = [cell_vol,Hmd_act,Mcr_isoenz,Q10Scale];

% SOLVE METABOLIC MODEL
[Rev,J_net,J_F,J_R] = metModel_main(Tc,H2,CO2,CH4,args);
% Order of reactions in marrices:
% (1) Fmd, (2) Ftr, (3) Mch, (4) Mtd, (5) Mer, (6) Mtr, (7) Mcr, (8) Frh,
% (9) Mvh/Hdr, (10) Hmd, (11) CO2 diffusion, (12) Mtd+Hmd
Rev   = squeeze(Rev);          % Reversibility of reactions
csMR  = squeeze(J_net(1,:,1)); % Cell-specific methanogenesis rates in fmol/cell/day
J_net = squeeze(J_net);
J_F   = squeeze(J_F);          % Gross forward fluxes of reactions
J_R   = squeeze(J_R);          % Gross reverse fluxes of reactions

% SOLVE ISOTOPIC MODEL
RVPDB    = 0.011202; % Standard carbon isotope ratio (VPDB)
RVSMOW   = 1.5576e-4; % Standard hydrogen isotope ratio (VSMOW)
R_13CCO2 = (d13CCO2./1000 + 1).*RVPDB; % Carbon isotope ratio of CO2
R_H2O    = (dDH2O./1000 + 1).*RVSMOW; % Hydrogen isotope ratio of H2O

% Find temperature dependent equilibrium isotope fractionation between H2O
% and H2
aH2Ol_H2_eq = 0.0334.*1e12./Tk.^4 - 0.2513.*1e9./Tk.^3 + ...
          1.0267.*1e6./Tk.^2 - 1.2166.*1e3./Tk + 1.7321;
R_H2 = R_H2O/aH2Ol_H2_eq; % Hydrogen isotope ratio of H2, assuming rapid equilibration between H2 and H2O. 

% Load distributions of kinetic fractionation factors (KIEs)
load('KFF_distributions_new.mat','KFF13C_FOR','KFF2H_FOR')

% load("a2H_post_ALL_DATA.mat")

% Calculate median KIE values
% randomly choose from the KIE distribution
nsims=200;
KFF_2H=zeros(1,14);
KFF_13C=zeros(1,8);
a13C=zeros(nsims,sims);
aD=zeros(nsims,sims);
e13C=zeros(nsims,sims);
eD=zeros(nsims,sims);
D13CH3D=zeros(nsims,sims);
D12CH2D2=zeros(nsims,sims);
for i=1:nsims
    idx_H=randi(length(KFF2H_FOR),1,14);
    idx_C=randi(length(KFF13C_FOR),1,8);
    for j=1:14
        KFF_2H(1,j)  = KFF2H_FOR(idx_H(j),j);
    end
    for k=1:8
        KFF_13C(1,k) = KFF13C_FOR(idx_C(k),k);
    end

    % CALCULATE ISOTOPE EFFECTS
    
    % CARBON ISOTOPES
    a13C(i,:) = isoModelCarb(Tc,Rev',KFF_13C);
    % HYDROGEN ISOTOPES
    aD(i,:) = isoModelHydr(Tc,J_F,J_R,KFF_2H,R_H2O,R_H2);
    %CLUMPED ISOTOPES
    D13CH3D(i,:) = isoModelClumped_13D(Tc,R_H2O,R_H2,R_13CCO2,Rev,J_F,J_R,KFF_2H,1./KFF_13C);
    D12CH2D2(i,:) = isoModelClumped_DD(Tc,R_H2O,R_H2,J_F,J_R,KFF_2H);
    
    % CONVERT ALPHA VALUES TO EPSILON
    e13C(i,:) = 1000.*(a13C(i,:)-1);
    eD(i,:) = 1000.*(aD(i,:)-1);
end

% Calculate average and 1-sigma uncertainties
y(1,:)=mean(e13C);
y(2,:)=mean(eD);
y(3,:)=mean(D13CH3D);
y(4,:)=mean(D12CH2D2);
y(5,:)=std(e13C);
y(6,:)=std(eD);
y(7,:)=std(D13CH3D);
y(8,:)=std(D12CH2D2);
y(9,:)=dGr;
output=y.'; % Transpose y
% PLOT FIGURE
plt_ylabels = {['' char(949) '_{CO_2-CH_4} (' char(8240) ')'],...
               ['' char(949) '_{CH_4-H_2O} (' char(8240) ')'],...
               ['\Delta{}^{13}CH_3D (' char(8240) ')'],...
               ['\Delta{}^{12}CH_2D_2 (' char(8240) ')']};
clf
x = -dGr;
EFF_lines      = calc_EFFs(Tc);
EFF_lines(1:2) = 1000.*(exp(EFF_lines(1:2)./1000)-1);

for i = 1:4
    ax = subplot(2,2,i);
    plot(x,y(i,:),'k')
    ylabel(plt_ylabels{i})
    if i > 2
        xlabel(['' char(8722) '\DeltaG_r (kJ mol^{' char(8722) '1})'])
    end
    if i == 1
        yline(EFF_lines(i),'--k','Equilibrium')
    elseif i > 1
        yline(EFF_lines(i),'--k')
    end
    ax.FontSize = 13;
    box off
end
