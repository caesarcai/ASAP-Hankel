function y=fhmvmultiply_1D(h,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fast Hankel matrix vector mulitplication.
%
%Inputs
% h: column vector that generates the Hankel matrix H.
% x: column vector to be left multiplied by H.
%
%Output
% y: H*x.  
%
%Example
% h=rand(128,1);x=rand(65,1);y=fhmvmultiply_1D(h,x), where y is the 
% multiplication of the Hankel matrix formed by h of size 64*65 with the 
% vector x.
%
%Reference: Lu L, Xu W, Qiao S. A fast SVD for multilevel block Hankel
%matrices with minimal memory storage[J]. Numerical Algorithms, 2015, 
%69(4): 875-891.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('.\PROPACK', 'dir')==7
    addpath PROPACK;
else
    fprintf('No PROPACK installed.\n');
    error('Break');
end

lh=length(h);

lx=length(x);
xrev=x(lx:-1:1);

L=2^nextpow2(lh);

fft_h=fft(h,L);
fft_xx=fft(xrev,L);

yy=ifft(fft_h.*fft_xx);
y=yy(lx:lh);