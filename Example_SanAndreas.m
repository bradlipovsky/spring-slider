%
% A simple rate- and state-friction spring slider system.  
%
% The code is unique in that implements quasi-dynamic elasticity 
% for a bi-material  interface as described in the previously 
% mentioned reference. doi:10.5194/tcd-9-1-2015.  
%

clear;
% Load the parameters
M = LoadParams('SanAndreas');

% How long to run the simulation (in seconds)
M.WindowDuration = 86400 * 365;

% Run the slider
tic; [vtr,dt,D,t] = RunSlider(M); tt=toc;
disp(['Simulation done in ' num2str(tt) ' s.']);

% Make plots
figure(1);
plot(t,detrend(vtr),'linewidth',2); 
xlabel('Time (s)'); ylabel('Fault Sliding Velocity (nm/s)');
set(gca,'fontsize',18); axis tight; 
yl=ylim; ylim( [min(yl(1)/2,1.25*yl(1)) 1.25*yl(2)]);