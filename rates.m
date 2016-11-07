function [V,dQ,n] = rates(D,Q,t,M,Vo)
%
% RATES sets slider velocity V and state Q
% by solving F(V)=stress-strength=0
%

% disp(num2str(t))

% if Q is incorrect (due to too large a time step),
% set rates to absurd values so error estimator will reject step
if isnan(Q),  V=1e9; dQ=1e9; n=0;  return, end ;

% Set tauL as a function of time
switch M.load
    case 'constant'
        tauL = M.alpha*t;
    case 'fit'
        tauL=ppval(M.AlphaFit,t) * k;
    otherwise
        disp('invalid M.load'); return;
end


%
% Solve for V in stress=strength with Newton's method,
% write as F(V)=stress(V)-strength(V) and solve F(V)=0
%

atol = 0;       % absolute tolerance for convergence
rtol = 1e-6;    % relative tolerance for convergence
nmax = 1000;

if nargin<6, Vo=M.V0; end % guess at V
V = Vo; % initial guess

for n=1:nmax % evaluate F = stress-strength
    stress = tauL-M.k*D-M.eta*V;

    % Standard log form of strength
%     strength = M.N*( M.a*log(V/M.V0) + Q );
    strength = M.N*( M.f0 + M.a*log(V/M.V0) + Q );

    % Regularized arcsinh form of strength
    % Omega = exp(Q/M.a)/(2*M.V0); 
    % strength = M.N*M.a*asinh(Omega*V);
    F = stress-strength;
    
% N* a* asinh( exp(Q/a) * V/(2*M.V0))
    % Test for convergence on F
    if abs(F)<atol
        dQ = evolveQ(V,Q,M);
        return
    end

    dstressdV = -M.eta;
    dstrengthdV = M.N*M.a/V; % log form
    %dstrengthdV = M.N*M.a*Omega/sqrt(1+(Omega*V)^2); % arcsinh form

    dFdV = dstressdV-dstrengthdV;
    dV = -F/dFdV; % correction    

    % check that update doesn't make V<=0, otherwise update to V/2
    if V+dV<=0, dV=-V/2; end
    
    % %
    % Variation of NM (VNM) proposed by Weerakoon and Fernando.
    %
    Vstar = V+dV;
    dstrengthdVstar = M.N*M.a/Vstar; % log form
    dFdVstar = dstressdV-dstrengthdVstar;
    dV2 = -2*F/(dFdV + dFdVstar);
    
    % check that update doesn't make V<=0, otherwise update to V/2
    if V+dV2<=0, dV2=-V/2; end
    V = V+dV2;  
    % V = V+dV; % NM
    %
    % End VNM
    % %
    

    % test for convergence on V

    if abs(dV)<atol+abs(V)*rtol
        dQ = evolveQ(V,Q,M);

        return
    end

end

dQ =  evolveQ (V,Q ,M);


end


function dQ = evolveQ(V,Q,M)

% slip law
% Remove f0 from the equation and put it in the definition of strength in
% the main function.  This allows the use of two state variables.
logV = log(V/M.V0);
dQ = - V/M.L*( Q + M.b*logV ); 

% ageing law
% dQ = M.b*M.V0/M.L*(exp((M.f0-Q)/M.b)-V/M.V0);

end

