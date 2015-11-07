function [t,D,Q,V,dQdt,sol,Q2] = slider(M,tmin,tmax)
% function [t,D,Q,V,dQdt,sol] = slider(M,tmin,tmax,Vini)

% initial conditions
% Q0 = M.f0-M.b*log(Vini/M.V0); % initial state
% T0 = (M.f0+(M.a-M.b)*log(Vini/M.V0))*M.N; % initial shear stress (MPa)
% D0 = -T0/M.k; % initial slip (m), set from initial stress

% integrate forward in time from tmin to tmax

% Precise
% rtol = 1e-9; % relative tolerance
% atol = eps; % absolute tolerance

% Medium
% rtol = 1e-9; % relative tolerance
% atol = 1e-9; % absolute tolerance

% Less precise
rtol = 1e-7; % relative tolerance
atol = 1e-7; % absolute tolerance

opt = odeset('InitialStep',M.dt,'Vectorized','on',...
    'RelTol',rtol,'AbsTol',atol);%,...
%     'OutputFcn',@odephas2);

% Two state version
if M.nstate == 2
    Y0 = [M.D0; M.Q0; M.Q20]; % intial conditions
elseif M.nstate == 1
    % One state version
    Y0 = [M.D0; M.Q0]; % intial conditions
end

% Solve it.
sol = ode45(@VQode,[tmin tmax],Y0,opt,M);

% evaluate fields
t = sol.x; % time vector (can replace with any t within [tmin tmax])
[Y,dY] = deval(t,sol);
D = Y(1,:); Q = Y(2,:); V = dY(1,:); dQdt = dY(2,:); 
if M.nstate == 2
    Q2 = Y(3,:);
else
    Q2 = nan;
end

function status = func(t,y,stat,varargin)

status=0;

if strcmpi(stat,'done')
elseif strcmpi(stat,'init')
else

    V = gradient(y(1,:),t);
    F = y(2,:);
    
    plot(V,F,'or');
    hold on; drawnow;
end