function [ft,f] = bft(x,dt)
%better fourier transform function

Fs = 1/dt;
L = numel(x);
NFFT = 2^nextpow2(L);
if NFFT>1
    Y = fft(detrend(x),NFFT)/L;
%     Y = fft( (x-mean(x)) ,NFFT)/L;
%     ft = 2*Y(1:NFFT/2+1);
    ft = Y(1:NFFT/2+1);
    f = Fs/2*linspace(0,1,NFFT/2+1) ;
else
    ft = Inf;
    f = Inf;
end