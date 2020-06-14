function U = conv_fft(U,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FFT based (multi-dimensional) convolution.
%
%Inputs
% U,V: (multi-dimensional) arrays of the same dimension.
%
%Output
% U: (multi-dimensional) convolution of U and V, saved in U to save memory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DIM1 = ndims(U);
DIM2 = ndims(V);

if DIM1 ~= DIM2
    error('Dimension Mismatch!');
else
    DIM = DIM1;
end

for dim = 1:DIM
    
    p = size(U,dim);
    q = size(V,dim);
    
    n = p+q-1;

    l = 2^nextpow2(n);
    
    U = fft(U,l,dim);
    V = fft(V,l,dim);
    
    Block{dim} = 1:n;
    
end
 
U = U.*V;
clear V;

for dim = 1:DIM
    U = ifft(U,[],dim);
end
  
U = U(Block{:});