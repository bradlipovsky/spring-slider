function M = LoadParams(varargin)
%
% This file sets the parameters for various cases.  The default parameters
% are listed first.  These are then modified for each particular case
% below.
%

verbose = 1;
if numel (varargin) == 0
    situ = 'CompliantBed';
else
    situ = varargin{1};
    if numel (varargin) > 1
        if strcmpi(varargin{2},'verbose')
            verbose = varargin{3};
        end
    end
end

% Default parameters

M.f0 = 0.4;     % static coefficient of friction
M.FarfieldVelocity = 1; % Interpret output velocities as far field velocities (set to zero to interpret as sliding velocities)
M.Q_ice = 500;  % quality factor
M.c_ice = 2000; % wavespeed
M.G_ice = 3664; % shear modulus (MPa)
M.H = 800;      % source-reciever distance
M.V0 = 10^-5;
M.dt = 0.5e-3;  % default time step
M.nui = 0.33;   % Poisson ratio in upper material
M.nub = 0.49;   % Poisson ratio in lower material
M.WindowDuration = 10;
M.a = 0.005;
M.b = 0.015;
M.G_till = 100;
M.rho_till = 1700;
M.L = 1e-6;

M.nstate=1;
M.b2 = 0;
M.L2 = M.L;

switch situ
    % Particular cases
    
        case 'SanAndreas'
   		M.FarfieldVelocity = 0; % Interpret velocities as sliding velocities (much faster...)
        M.Vs = 0.4/3.14e7;      % Far-field loading rate 40cm/yr
        M.c_ice = 4000;         % wavespeed
        M.G_ice = 50e4;         % Shear modulus (MPa)
        M.G_till = M.G_ice;     % Shear modulus (MPa) is identical on both sides
        M.R = 10e3;             % Fault Size (m)
        M.N = 200;              % Effective normal stress (MPa)
        M.L = 1e-3;             % Frictional state evolution distance (aka D sub C)
        
        case 'DavidGlacier'
            
        M.L = 1e-4;
        M.G_till = 3664; 
        M.R = 1500;
        M.N = 1;
        
        M.H = 180e3; % Source-to-receiver distance
        M.Q_ice = 100;  % quality factor
        M.Vs = 1.6e-5;
        M.G_ice  = M.G_till;
        M.nui = 0.25;   M.nub = M.nui;
        M.rho_till = 2000; M.c_ice = 2000;
      
        case 'CompliantBed'
        M.Vs = 4.86e-4;
        M.G_till = 10;
        M.R = 3.3;
        M.N = 0.04;
        
        
        case 'RigidBed'
        M.Vs = 4.86e-4;
        M.R = 3.3;
        M.N = 0.1;
        M.D = M.Vs/M.f1;
        
        
        case 'SubductionZone'
        
        % Both sides of the fault have the same material
        M.G_ice = 1e4; M.G_till = M.G_ice;
        M.nui = 0.25; M.nub = 0.25;
        
        % 10m fault radius
        M.R = 10;
        
        % Normal stress is 1 kPa
        M.N = 1;
        
        % source-reciever distance
        M.H = 15e3;
        


end




% Derived quantities
M = DerivedParams(M);



end

