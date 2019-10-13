clear all;
close all;

n     = 10000;   % signal size
r     = 10;      % rank
alpha = 0.1;     % amount of corruptions
c = 1;           % magnitude of corruptions

% Generate an 1D spectrally sparse signal
[~,ox,f] = generate_signal(n,r,n);

% Add corruptions
temp = rand(1,n);
IND = find(temp<alpha);
m = length(IND);
os = zeros(m,1);
a = c*mean(abs(real(ox)));
b = c*mean(abs(imag(ox)));
for i = 1:m
    v1 = a-2*a*(1-rand(1));
    v2 = b-2*b*(1-rand(1));
    os(i) = v1+1i*v2;
end

z  = ox;
z(IND) = ox(IND)+os;

% Solving the 2D problem
[x,s] = AAP_Hankel_1D(z,r,0.5);
recovery_err = norm(x-ox,'fro')/norm(ox,'fro')
