function M = LoadParams(varargin)

% M.Vs = 7e-4; % 100s GPS lowpass (Event 49)
% M.Vs = 6.02e-4; % 1000s GPS lowpass (Event 49)
% M.Vs = 4.5e-4; % 100s GPS lowpass (Event 29)

verbose = 1;
if numel (varargin) == 0
    situ = 'CompliantBed2';
else
    situ = varargin{1};
    if numel (varargin) > 1
        if strcmpi(varargin{2},'verbose')
            verbose = varargin{3};
        end
    end
end

% These parameters are usually not changed
M.f0 = 0.4;
M.Q_ice = 500;
M.c_ice = 2000;
M.G_ice = 3664;
M.H = 800;
M.V0 = 10^-5;
M.dt = 0.5e-3;
M.nui = 0.33;
M.nub = 0.49;
M.WindowDuration = 10;

% Default Tunable parameters
M.a = 0.005;
M.b = 0.015;
M.G_till = 100;
M.rho_till = 1700;
M.L = 1e-6;

M.nstate=1;
M.b2 = 0;
M.L2 = M.L;

FullSim=0;
switch situ
    % Tunable parameters

      
        case 'CompliantBed'
        
        % f1, Vs, and D are measured from the data and are 
        %  not used in the simulation, 
        M.f1 = 14.25;
        M.Vs = 4.86e-4; % f0=13.7
        M.D = M.Vs/M.f1;
        
        M.G_till = 10;
        M.R = 3.3;
        M.N = 0.01;
        
        
        
        
        case 'RigidBed'
        M.f1 = 14.25;
        M.Vs = 4.86e-4; % f0=13.7
        M.G_till = 100;
        
        M.R = 3.3;
        M.N = 0.1;
        M.D = M.Vs/M.f1;


end




% Derived quantities
M = DerivedParams(M);

if ~verbose
    return;
end

%
% Check for stable sliding
%
disp(' ');
disp(' ');
disp('----------------------------------------------');
if ~FullSim
StableSliding = (M.b-M.a)*M.N <= M.eta * M.Vs + M.k*M.L;
if StableSliding
    tr = 1/eps; dt= 0;
    disp('Sliding is stable because...');
    w=0;
    if (M.b-M.a)*M.N < M.k*M.L
        Rc = M.k*M.R*M.L/(M.b-M.a)/M.N;
        disp(['   The patch size ' num2str(M.R,3) ...
            ' is less than the critical patch size ' num2str(Rc,3) '.']);
        w=1;
    end
    if (M.b-M.a)*M.N < M.eta * M.Vs
        Vc = (M.b-M.a)*M.N / M.eta;
        disp(['   V = ' num2str(M.Vs,3) ...
            ' > Vc =  ' num2str(Vc,3)]);
        w=1;
    end
    if w==0
        disp('   The stability is in the intermediate regime.');
    end
else
    disp('----------------------------------------------');
disp(' ');
disp(' Sliding is in the unsteady regime.');
disp(' ');
end


% Useful params
Fc = M.k/M.eta;
Vc = (M.b-M.a)*M.N / M.eta - Fc*M.L;
D = M.N*(M.b-M.a)/M.k;
Z = M.z_till/(M.z_ice+M.z_till);

% vel amp
V = VeloAmplitudeScaling(M.G_till,M);

disp('----------------------------------------------');
disp(' ');
disp(['   Vc                ' num2str(Vc,3)]);
disp(['   D (stress drop)   ' num2str(D*1e6,3)]);
disp(['   D (data)          ' num2str(M.Vs/M.f1*1e6,3)]);
disp(['   Fc                ' num2str(Fc,3)]);
disp(' ');
disp(['   Z                 ' num2str(Z,3)]);
disp(['   Eta               ' num2str(M.eta,3)]);
disp(['   k                 ' num2str(M.k,3)]);
disp(' ');
disp('----------------------------------------------');
disp(' ');
disp(['   Inertial          ' num2str(M.eta * M.Vs*1e6,3)]);
disp(['   Static            ' num2str(M.k*M.L*1e6,3)]);
disp(['   Stress Drop       ' num2str((M.b-M.a)*M.N*1e6,3)]);
disp(' ');
disp('----------------------------------------------');
disp(' ');
disp([' Far fie par vel amp ' num2str(V,3)]);
disp(' ');
disp('----------------------------------------------');


end

