function y = fhmvmultiply_2D(h,x,q1,q2,ind1,ind2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fast two-level Hankel matrix vector mulitplication.
%
%Inputs
% h: 2D matrix that generates the two-level Hankel matrix H.
% x: column vector to be left multiplied by H.
% q1: column dimension of level 1.
% q2: column dimension of level 2.
% ind1 (optional): index for x after padding and reversing order. 
% ind2 (optional): index for y, to be extracted after ifft.
%
%Output
% y: H*x.  
%
%Example
% h = rand(64,32);x = rand(33*17,1);y = fhmvmultiply_2D(h,x,33,17), where y is 
% the multiplication of the two-level Hankel matrix of size (32*16)*(33*17) 
% formed by h with the vector x.
%
%Reference: Lu L, Xu W, Qiao S. A fast SVD for multilevel block Hankel
%matrices with minimal memory storage[J]. Numerical Algorithms, 2015, 
%69(4): 875-891.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,n2] = size(h);

lh = n1*n2;
h = reshape(h,lh,1);

lx = length(x);
xrev = x(lx:-1:1);

p1 = n1+1-q1;
p2 = n2+1-q2;

nargs = nargin;
if nargs == 4
    ind1 = zeros(q1*q2,1);
    for i = 1:q2
        ind1((i-1)*q1+1:i*q1) = (i-1)*n1+1:(i-1)*n1+q1;
    end
    ind2 = zeros(p1*p2,1);
    for i = 1:p2
        ind2((i-1)*p1+1:i*p1) = (q2+i-2)*n1+q1:(q2+i-1)*n1;
    end
end

xx = zeros(lh,1);
xx(ind1) = xrev;

L = 2^nextpow2(lh);

fft_h = fft(h,L);
fft_xx = fft(xx,L);

yy = ifft(fft_h.*fft_xx);
y = yy(ind2);