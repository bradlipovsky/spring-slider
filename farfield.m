function [v,omega,vhat] = farfield(A,dt,R,r,c,Q)

% given slip acceleration A on circular fault of radius R,
% stored in vector A at equally spaced time steps dt,
% return ground velocity v=v(t) and spectrum vhat=vhat(omega),
% where omega is angular frequency (omega/(2*pi) is frequency)
%
% source is located distance r from receiver in medium with
% wave speed c and quality factor Q

    N = length(A);
    omegaN = pi/dt; % Nyquist angular frequency
    omega = zeros(1,N);
    omega(1:N/2+1) = 2*[0:N/2]/N;
    omega(N:-1:N/2+2) = -omega(2:N/2);
    omega = omega*omegaN;
    Ahat=dt*fft(A);

%     MysteryNumber = 0.001;
    MysteryNumber = dt;
    tstar=r/(c*Q); 
    omega0=2*pi/MysteryNumber; % omega0 is frequency at which phase velocity=c
    Ahat = Ahat.*exp(-abs(omega)*tstar/2);
    I = 2:N; Ahat(I) = Ahat(I).*exp(-1i*omega(I)*tstar/pi.*log(omega0./abs(omega(I))));
    Ahat(N/2+1)=0;

    % Leave the geometrical spreading factor out for now
%     vhat=R^2/(4*c*r)*Ahat; 
    vhat = Ahat;
    vhat(N/2+1)=0;

    v=ifft(vhat)/dt;
