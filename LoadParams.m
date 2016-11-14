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
            
        M.Vs = 0.4/3.14e7;       % Far-field loading rate
        M.Q_ice = 100;          % quality factor
        M.c_ice = 4000;         % wavespeed
        M.G_ice = 50e4;         % Shear modulus (MPa)
        M.G_till = M.G_ice;     % Shear modulus (MPa) is identical on both sides
        M.R = 200;              % Fault Size (m)
        M.N = 1;                % Effective normal stress (MPa)
        M.H = 10e3;             % Source-to-receiver distance
        M.L = 1e-4;             % Frictional state evolution distance
        
        case 'DavidGlacier'
        
        M.Vs = 4.86e-4;
        M.G_till = 100;
        M.R = 200;
        M.N = 0.04;
        M.H = 10e3; % Source-to-receiver distance
        M.Q_ice = 100;  % quality factor

      
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

