function [k,eta] = bimat(G1,G2,nu1,nu2,rho1,rho2)
% Returns the "patch modulus," e.g., the spring constant times the patch 
%  size.  For homogeneous materials this gives the Eshelby solution,
%
%       16  (1 - nu1)  G1
%     --------------------
%         3 pi (2 - nu1)
%
% Also returns the radiation damping parameter which depends on the shear
% impedence of each material.  This is the only dependence on rho1 and
% rho2.
%
% Currently, G2 can be a vector, G1 cannot.

nG = numel(G2);
C = zeros(nG,1);

for i = 1:nG;
    G2a = G2(i);

    b = (1 - nu1)/2/pi/G1 + (1 - nu2)/2/pi/G2a;
    c = nu1/2/pi/G1 + nu2/2/pi/G2a;
    d = (1 - 2*nu1)/4/pi/G1 - (1 - 2*nu2)/4/pi/G2a;
    ee = (b*c + d^2) / (b+c);

    k = 1/2/pi * log ( (b+d) / (b-d));

    if d < 1e-6
        % use correct taylor series
        SpecialTerm = pi*b;
    else
        SpecialTerm = d/k/(1+k^2);
    end

    % C is the compliance
    C(i) = 16*pi/3* (b^2-d^2) / (  pi*(b-ee)  +  SpecialTerm  );

end

k = 1./C;
% esh = 16/3/pi * (1-nu1)/(2-nu1) ./ G1;
% k = esh./k;
% k = 1./k;

c1 = sqrt(G1/rho1);
c2 = sqrt(G2/rho2);
z1 = G1/c1;
z2 = G2/c2;
eta = 1./(1./z1 + 1./z2);
eta = eta';




% Plot I gave to Eric on 7/22/15. 
% val = C;
% figure(1);
% semilogx(G2,c./esh,'-k'); hold on;
% xlim([min(G2) max(G2)]);
% xlabel('Shear modulus ratio G_+/G_-');
% ylabel('Normalized slip per unit applied stress');
% 
% ll=log(3-4*nu1);
% gam=(ll * (1 + ll^2/4/pi^2))^-1;
% gam2 = 1 / (1 + 4*gam*(1-2*nu1)/(3-4*nu1)) ;
% val = 8/3/pi  * gam2;
% 
% line(xlim,val/esh*[1 1],'linesty','--','color','k');