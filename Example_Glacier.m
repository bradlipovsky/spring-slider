%
% A simple rate- and state-friction spring slider system.  
%

%
% Please cite Lipovsky and Dunham (2016) if using this code 
% (doi:10.5194/tcd-9-1-2015).  
%
% The code is unique in that implements quasi-dynamic elasticity 
% for a bi-material  interface as described in the previously 
% mentioned reference.  
%
% This code is based on an earlier code by Dmitrieva et al. 
% (2015, doi:10.1038/NGEO1879).
%


clear;

% Load the parameters
% M = LoadParams('RigidBed');
M = LoadParams('DavidGlacier');

% Run the slider
M.SimulationDuration = 60*60*2;
tic; [vtr,dt,D] = RunSlider(M); tt=toc;
disp(['Simulation done in ' num2str(tt) ' s.']);

% Make a time vector
t = dt*(0:numel(vtr)-1);

% Take the fourier transform
[ft,f] = bft( vtr,dt );

% Make plots
figure(1);
plot(t/60,vtr,'linewidth',2); 
xlabel('Time (min)'); ylabel('Seismic Particle Velocity (nm/s)');
set(gca,'fontsize',18); axis tight; hold on;
yl=ylim; ylim( [min(yl(1)/2,1.25*yl(1)) 1.25*yl(2)]);
