function xx = cadzow_2D(x,r)
% cadzow denoising for 2D signal

[n1,n2] = size(x);

if mod(n1,2)
    p1 = (n1+1)/2;
    d1 = [1:p1 p1-1:-1:1]';
else
    p1 = n1/2;
    d1 = [1:p1 p1 p1-1:-1:1]';
end
if mod(n2,2)
    p2 = (n2+1)/2;
    d2 = [1:p2 p2-1:-1:1]';
else
    p2 = n2/2;
    d2 = [1:p2 p2 p2-1:-1:1]';
end

DD = kron(d2,d1);
DD = reshape(DD,n1,n2);

q1 = n1+1-p1;
q2 = n2+1-p2;

l1 = p1*p2;
l2 = q1*q2;

% indecies pre-computed for fhmvmultiply_2D to use
ind1=zeros(l2,1);
for i=1:q2
    ind1((i-1)*q1+1:i*q1) = (i-1)*n1+1:(i-1)*n1+q1;
end
ind2 = zeros(l1,1);
for i=1:p2
    ind2((i-1)*p1+1:i*p1) = (q2+i-2)*n1+q1:(q2+i-1)*n1;
end
ind3=zeros(l1,1);
for i=1:p2
    ind3((i-1)*p1+1:i*p1) = (i-1)*n1+1:(i-1)*n1+p1;
end
ind4 = zeros(l2,1);
for i=1:q2
    ind4((i-1)*q1+1:i*q1) = (p2+i-2)*n1+p1:(p2+i-1)*n1;
end

opts = [];
opts.eta = 1e-16; % set tolerance for lansvd

xx = x;

for iter = 1:10
    
    Yforward = @(z) fhmvmultiply_2D(xx,z,q1,q2,ind1,ind2);
    Ytranspose = @(z) fhmvmultiply_2D(conj(xx),z,p1,p2,ind3,ind4);
    [U,S,V] = lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);
    s = diag(S);
    
    xx = zeros(n1,n2);
    for i = 1:r
        ui = reshape(U(:,i),p1,p2);
        vi = reshape(V(:,i),q1,q2);
        xx = xx+s(i)*conv_fft(ui,conj(vi));
    end
    xx = xx./DD;
    
    fprintf('%2d \n',iter)
end