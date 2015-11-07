function M = DerivedParams(M)

% Derived quantities
M.c_till = sqrt(M.G_till*10^6/M.rho_till);

M.z_till= M.G_till/M.c_till;
M.z_ice = M.G_ice/M.c_ice;

M.eta = 1 / ( 1/M.z_till + 1/M.z_ice );

M.fs=1/M.dt;



% Use the bimaterial solution for the spring stiffness
M.k = bimat(M.G_ice,M.G_till,M.nui,M.nub,1e6*M.G_ice/M.c_ice^2,M.rho_till)./M.R; 

% Wrong bimaterial solution
% M.k = 7*pi*M.G./(16*M.R);
% M.G = 2 * (   1/M.G_ice + 1 / M.G_till   )^-1;
% M.c = M.G/M.eta;