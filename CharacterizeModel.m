clear;


% Load the parameters
% M = LoadParams('RigidBed');
% M = LoadParams('CompliantBed','verbose',1);
% M = LoadParams('CompliantBed');

% Run the slider
M.WindowDuration = 2;
tic; [vtr,dt,D] = RunSlider(M); tt=toc;
disp(['Simulation done in ' num2str(tt) ' s.']);

% Make a time vector
t = dt*(0:numel(vtr)-1);

% Take the fourier transform
[ft,f] = bft( vtr,dt );

% Make plots
figure(1); clf;
subplot(1,2,1);
plot(t,detrend(vtr),'linewidth',2); 
xlabel('Time (s)'); ylabel('Seismic Particle Velocity (nm/s)');
set(gca,'fontsize',18); axis tight; 
yl=ylim; ylim( [min(yl(1)/2,1.25*yl(1)) 1.25*yl(2)]);

subplot(1,2,2);
plot(f,abs(ft),'-','linewidth',2);
xlim([5 90]); ylim([0 1.25*max(abs(ft(f<90)))]); % Event 49
xlabel('Frequency (Hz)'); ylabel('Spectral Power (m/s)^2 / Hz');
set(gca,'fontsize',18);