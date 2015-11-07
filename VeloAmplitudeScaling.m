function [V, Fc] = VeloAmplitudeScaling(Gb,M)

Gb = abs(Gb);

rhoi = M.G_ice / M.c_ice^2 * 1e6;
zi = M.G_ice/sqrt(M.G_ice*1e6/rhoi);
zb = Gb./sqrt(Gb*1e6/M.rho_till);
Z = zb./(zi+zb);

Gs = bimat(M.G_ice,Gb,M.nui,M.nub,rhoi,M.rho_till);
eta = 1/(1/zi + 1/zb);
% Fc = sqrt(2) *Gs./eta;
Fc = Gs./eta;
V = M.D / M.H / M.c_ice * Z' .* Fc.^2 * 1e9;

% disp (num2str(Gb));
% disp (num2str(abs(V-100)));
% disp(' ');