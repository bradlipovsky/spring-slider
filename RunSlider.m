function [tr,dt,De] = RunSlider(M)

% specify loading
M.load = 'constant';
M.alpha = M.k.*M.Vs;

% Use the bimaterial solution for the spring stiffness
% M.k = 7*pi*M.G./(16*M.R);
M.k = bimat(M.G_ice,M.G_till,M.nui,M.nub,1e6*M.G_ice/M.c_ice^2,M.rho_till)./M.R; 

% Solve problem, determine slip and recurrence interval
tmin = 0;                       % start time (s)
tmax0 = M.WindowDuration + 1/M.fs;

tmax=tmax0;

% Set initial conditions
Vini = M.Vs;
M.Q0 =  -M.b*log(Vini/M.V0); % initial state
T0 = (M.f0+(M.a-M.b)*log(Vini/M.V0))*M.N; % initial shear stress (MPa)
M.D0 = -T0/M.k; % initial slip (m), set from initial stress

% Run the simulation
[t,D,Q,V,dQdt,~] = slider(M,tmin,tmax);

T = M.a.*log(V/M.V0) + Q;

% Check to make sure that stick-slip events are happening.
% T = calcT(t,V);
% if T == 0
%     tr = 1/eps; dt= 0;
%     disp('Solutions are stable; no earthquakes are going to occur!');
%     return;
% end
    
% calculate the shear stress
% tau = -M.k*D - M.eta*V;

% High-frequency amplitude calculation
fs2 = 1/min(diff(t));  dt = 1/fs2; 
t2 = min(t):dt:max(t);
N2 = numel(t2); if mod(N2,2), t2 = min(t2):dt:max(t2)-dt; end
A = (M.alpha-M.k*V-M.N*dQdt)./(M.a*M.N./V+M.eta);   % slip acceleration A=dV/dt
A_interp = interp1(t,A,t2);


% Realistic Q values in Whillans suggest that Q dosen't matter.  So it's
% easier to calculate the farfield velocity this way:
v2 = M.z_till/(M.z_ice + M.z_till) * M.R^2/(M.c_ice*M.H) * farfield(A_interp,dt,M.R,M.H,M.c_ice,M.Q_ice);
% correctly account for the bimaterial interface (inf Q)
% v2 = M.z_till/(M.z_ice + M.z_till) * M.R^2/(M.c_ice*M.H) * A_interp;

% Don't decimate.
tr = v2;
% time = t2;



% Decimate to the seismometer frequency
% tr = decimate(v2,round(1/(dt*M.fs)),12);
% tr = resample(v2,M.fs,round(fs2));
% dt = 1/M.fs;

% Convert to nm/s
tr = tr * 10^9;


% Return D/event
[~,val] = findpeaks(V);
De = mean(diff(D(val)));


%
% Plots
%

% % Predicted far field velocity amplitude
% [ft,f] = bft(v3,dt);
% [~,fmax] = max(ft);
% Fc = f(fmax);
% vff = M.R^2/(4*M.c_ice*M.H) * M.Vs * Fc;
% disp(['Predicted max amplitude:  ' num2str(vff)]);

% [ft,f] = bft(tr,dt);
% figure(1); plot(f, abs(ft)/max(abs(ft)) ,'-k');
% time = min(t2):dt:max(t2);
% figure; plot(time,tr,'-ob');
% figure; plot(t,D);


% Q-V space phase portrait
% figure(3); clf;
% tau = M.N*(M.f0 + M.a*log(V/M.V0)+Q);
% Vstar = V/M.Vs;
% % Tss = M.N * (M.f0 - (M.b-M.a)*log(M.Vs/M.V0));
% Tss = M.N * (M.f0 - (M.b-M.a)*log(V/M.V0));
% TauStar = (tau - Tss)/(M.a/M.N);
% % semilogx(Vstar,TauStar);hold on; 
% % vl = logspace(log10(min(abs(Vstar))),log10(max(abs(Vstar))),10);
% 
% nQ = numel(Q);
% col = parula(nQ);
% for i = 1:nQ
%     semilogx(Vstar(i),TauStar(i),'o','markerfacecolor',col(i,:),...
%         'markeredgecolor',col(i,:));hold on; 
% end

% figure;
% for i = 1:nQ
%     plot(A(i) * M.L/M.Vs^2,TauStar(i),'o','markerfacecolor',col(i,:),...
%         'markeredgecolor',col(i,:));hold on; 
% end

% The steady state line is now x=0.
% tau_ss = M.N*(M.f0 - (M.b - M.a)*log(vl/M.V0));
% tau_ss = tau_ss - Tss)/A;
% semilogx(vl,tau_ss); axis tight;

% xlabel('V/V_L');
% ylabel('[ \tau - \tau_{SS} ] / a\sigma');