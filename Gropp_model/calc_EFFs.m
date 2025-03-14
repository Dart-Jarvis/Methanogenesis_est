function EFFs = calc_EFFs(Tc)

% A function to find EFFs based on DFT calculations in the Mtaken from
% Gropp et al., GCA (doi.org/10.1016/j.gca.2020.10.018), in the
% M06-L/def2-TZVP level of theory with the SMD implicit solvation model

% Input:
%         Tc - temperature in degree celcuis
% Output:
%         EFFs - a matrix of 4 X Tc. Ordered by 1000ln13a, 1000ln2a, D13CH3D and D12CH2D2.          

% Coefficients (updated as of 1.1.2020)
b_vals = [ 0.17825   -2.03421    10.3635    2.55042    -2.26002;
          -5.39537   53.2932   -205.032   319.039   -286.391;
           0.0307659 -0.370409    1.80240  -1.42608    0.321113;
          -0.1103590  0.734810    1.15068  -2.83633    1.319360];
Tk = Tc+273.15;

% Calculate EFFs
values = b_vals(:,1).*1e12./Tk.^4 + b_vals(:,2).*1e9./Tk.^3 + ...
         b_vals(:,3).*1e6./Tk.^2  + b_vals(:,4).*1e3./Tk + ...
         b_vals(:,5);
     
kln13a   = values(1,:);
kln2a    = values(2,:);
D13CH3D  = values(3,:);
D12CH2D2 = values(4,:);
EFFs    = [kln13a;kln2a;D13CH3D;D12CH2D2];    