function dY = VQode(t,Y,M)
  
% VQODE returns V, dQ, dV as a column vector
      
D = Y(1);
Q = Y(2);

% Two-state version
if M.nstate == 2
    Q2= Y(3);
    [V,dQ,dQ2,~] = rates(D,Q,Q2,t,M);
    dY = [V; dQ; dQ2];
else
    % One-state version
    [V,dQ,~] = rates(D,Q,t,M);
    dY = [V; dQ];
end
  
t
