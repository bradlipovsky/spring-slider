function [tr,dt,De,t] = RunSlider(M)

% Integration tolerances
rtol = 1e-7; % relative tolerance
atol = 1e-7; % absolute tolerance

% specify loading
M.load = 'constant';
M.alpha = M.k.*M.Vs;

% Use the bimaterial solution for the spring stiffness
% M.k = 7*pi*M.G./(16*M.R);             % "Eshelby"
M.k = bimat(M.G_ice,M.G_till,M.nui,M.nub,1e6*M.G_ice/M.c_ice^2,M.rho_till)./M.R; 

% Set initial conditions
M.Q0    = -M.b*log(M.Vs/M.V0);          % initial state
T0      = (M.f0+(M.a-M.b)*log(M.Vs/M.V0))*M.N; % initial shear stress (MPa)
M.D0    = -T0/M.k;                      % initial slip (m), set from initial stress
Y0      = [M.D0; M.Q0];                 % intial conditions

% Run the slider!
opt = odeset('InitialStep',M.dt,'Vectorized','on',...
    'RelTol',rtol,'AbsTol',atol,'OutputFcn',@outfunc);
sol = ode45(@VQode,[0 M.SimulationDuration],Y0,opt,M);

% Evaluate fields.
t = sol.x; 
[Y,dY] = deval(t,sol);
D = Y(1,:); Q = Y(2,:); V = dY(1,:); dQdt = dY(2,:); 




%
% Calculate Velocities
%
% The velocity returned in "tr" can either be interpreted as a far-field velocity (e.g., as recorded by a seismometer), or it can be interpreted as a fault sliding velocity.  
% If interpreted as a far-field velocity (the default), then wave propagation is calculated with attenuation.  If a sliding velocity, then we just return the output of the ode solver.

if M.FarfieldVelocity
	% High-frequency amplitude calculation
	fs2 = 1/min(diff(t));  dt = 1/fs2; 
	t2 = min(t):dt:max(t);
	N2 = numel(t2); if mod(N2,2), t2 = min(t2):dt:max(t2)-dt; end
	A = (M.alpha-M.k*V-M.N*dQdt)./(M.a*M.N./V+M.eta);   % slip acceleration A=dV/dt
	A_interp = interp1(t,A,t2);


	% Valid for large Q
	% v2 = M.z_till/(M.z_ice + M.z_till) * M.R^2/(M.c_ice*M.H) * farfield(A_interp,dt,M.R,M.H,M.c_ice,M.Q_ice);

	% Account for the bimaterial interface (inf Q)
 	v2 = M.z_till/(M.z_ice + M.z_till) * M.R^2/(M.c_ice*M.H) * A_interp;

	% Don't decimate.
	tr = v2;
	% time = t2;
	
	% Decimate to the seismometer frequency
	% tr = decimate(v2,round(1/(dt*M.fs)),12);
	% tr = resample(v2,M.fs,round(fs2));
	% dt = 1/M.fs;
else
	tr = V;
    dt = -1; % return the full time vector instead.
end



% Convert to nm/s
tr = tr * 10^9;


% Return slip (D) per event as De
[~,val] = findpeaks(V);
De = mean(diff(D(val)));






%
% Some nice plots to make
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
