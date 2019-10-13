function [xx,U,V] = cadzow_1D(x,r)
% cadzow denoising for 1D signal

n = size(x,1);

if mod(n,2)
    p=(n+1)/2;
    DD=[1:p p-1:-1:1]';
else
    p=n/2;
    DD=[1:p p p-1:-1:1]';
end

q=n+1-p;

opts=[];opts.eta=1e-16; % set tolerance for lansvd

xx=x;

for iter = 1:50
    
    Yforward=@(z) fhmvmultiply_1D(xx,z);
    Ytranspose=@(z) fhmvmultiply_1D(conj(xx),z);
    [U,S,V]=lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
    s=diag(S);
    
    xx=zeros(n,1);
    for i=1:r
        ui=U(:,i);
        vi=V(:,i);
        xx=xx+s(i)*conv_fft(ui,conj(vi));
    end
    xx=xx./DD;
    
    fprintf('%2d \n',iter)
    
end