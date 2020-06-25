function [x,s] = ASAP_Hankel_2D(z,r,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x,s] = AAP_Hankel_2D(z,r,gamma)
% 
% Inputs:
% z     : Observed 2D signal. 
% r     : Target rank of underlying low rank 2-level Hankel matrix.
% gamma : Parameter for desired convergence rate. Value should between 0
%            and 1. Turn this parameter bigger will slow the convergence
%            speed but tolerate harder problem, such as higher p, r or mu.    
%
% Outputs:
% x : Estimated 2D spectrally r-sparse signal.
% s : Estimated sparse corruptions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('.\PROPACK', 'dir')==7
    addpath PROPACK;
    propack_exist = true;
else
    propack_exist = false;
    fprintf('No PROPACK installed.\n');
    return;
end

max_iter = 100;
err    = -1*ones(max_iter,1);
timer  = -1*ones(max_iter,1);
tol    = 1e-10; 
norm_of_z = norm(z,'fro');

[n1,n2] = size(z);

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

N = n1*n2;
c_s = max(N/l1,N/l2);

% indecies pre-computed for fhmvmultiply_2D to use
ind1 = zeros(l2,1);
for i = 1:q2
    ind1((i-1)*q1+1:i*q1) = (i-1)*n1+1:(i-1)*n1+q1;
end
ind2 = zeros(l1,1);
for i = 1:p2
    ind2((i-1)*p1+1:i*p1) = (q2+i-2)*n1+q1:(q2+i-1)*n1;
end
ind3 = zeros(l1,1);
for i = 1:p2
    ind3((i-1)*p1+1:i*p1) = (i-1)*n1+1:(i-1)*n1+p1;
end
ind4 = zeros(l2,1);
for i = 1:q2
    ind4((i-1)*q1+1:i*q1) = (p2+i-2)*n1+p1:(p2+i-1)*n1;
end

opts = []; opts.eta = 1e-16;

% use Cadzow to estimate incoherence and sigma_L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(n1,n2);
Yforward = @(y) fhmvmultiply_2D(z,y,q1,q2,ind1,ind2);
Ytranspose = @(y) fhmvmultiply_2D(conj(z),y,p1,p2,ind3,ind4);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r)); sigma_D = ss(1);
U = U(:,1:r);
V = V(:,1:r);

for i = 1:r
    ui = reshape(U(:,i),p1,p2);
    vi = reshape(V(:,i),q1,q2);
    x = x+ss(i)*conv_fft(ui,conj(vi));
end
x = x./DD;

Yforward = @(y) fhmvmultiply_2D(x,y,q1,q2,ind1,ind2);
Ytranspose = @(y) fhmvmultiply_2D(conj(x),y,p1,p2,ind3,ind4);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r)); sigma_L = ss(1);
U = U(:,1:r);
V = V(:,1:r);

row_norm = zeros(l1,1);
for i = 1:l1
    row_norm(i) = norm(U(i,:))^2;
end

col_norm = zeros(l2,1);
for j = 1:l2
    col_norm(j) = norm(V(j,:))^2;
end

mu = max(max(row_norm),max(col_norm))*N/(c_s*r);
 
beta = mu*c_s*r/(2*N);

c_init = 2; % between 1 and 3
beta_init = c_init*mu*c_s*r*sigma_L/(N*sigma_D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% two-step initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
eta = beta_init*sigma_D;
ind = find(abs(z) > eta);
s = zeros(n1,n2);
s(ind) = z(ind);

x = z-s;
Yforward = @(y) fhmvmultiply_2D(x,y,q1,q2,ind1,ind2);
Ytranspose = @(y) fhmvmultiply_2D(conj(x),y,p1,p2,ind3,ind4);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r));
U = U(:,1:r);
V = V(:,1:r);

x = zeros(n1,n2);
for i = 1:r
    ui = reshape(U(:,i),p1,p2);
    vi = reshape(V(:,i),q1,q2);
    x = x+ss(i)*conv_fft(ui,conj(vi));
end
x = x./DD;

eta = beta*ss(1);
temp = z-x;
ind = find(abs(temp)>eta);
s = zeros(n1,n2);
s(ind) = temp(ind);
init_timer = toc;
init_err = norm(z-x-s,'fro')/norm_of_z;
fprintf('Init: error: %e, timer: %f \n', init_err, init_timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


threshold = mu*c_s*r/N;

% iterations
for iter = 1:max_iter
    tic;
    % trim
    [U,V] = trim(U,diag(ss),V,threshold,threshold);
    
    % estimate the signal
    x = z-s;
    
    UtHx = zeros(r,l2);
    HxV = zeros(l1,r);
    Z = zeros(r,r); 
    
    for i = 1:r
        ui = U(:,i);
        UtHx(i,:) = (fhmvmultiply_2D(conj(x),ui,p1,p2,ind3,ind4))';
        vi = V(:,i);
        HxV(:,i) = fhmvmultiply_2D(x,vi,q1,q2,ind1,ind2);
    end
    
    C = UtHx*V;
    Xt = UtHx-C*(V');
    X = Xt';
    Y = HxV-U*C;
    [Q1,R1] = qr(X,0);
    [Q2,R2] = qr(Y,0);
    M = [C R1';R2 Z];
    
    [Uc,SS,Vc] = svd(M);
    ss = diag(SS(1:r,1:r));
    Uc = Uc(:,1:r);
    Vc = Vc(:,1:r);
    
    U = [U Q2]*Uc;
    V = [V Q1]*Vc;
    
    x = zeros(n1,n2);
    for i = 1:r
        ui = reshape(U(:,i),p1,p2);
        vi = reshape(V(:,i),q1,q2);
        x = x+ss(i)*conv_fft(ui,conj(vi));
    end
    x = x./DD;
    
    
    % estimate the outlier
    eta = beta*(gamma^iter)*SS(1,1);
    s = wthresh( z - x ,'h', eta);
    
    err(iter) = norm (z-x-s,'fro')/norm_of_z;
    timer(iter) = toc;
    
    if err(iter) < tol
        fprintf('Total %d iteration, final error: %e, total time without init: %f , with init: %f\n======================================\n', iter, err(iter), sum(timer(timer>0)),sum(timer(timer>0))+init_timer);       
        return;
    elseif err(iter) > 1.5 % blow up
        return;
    else
        fprintf('Iteration %d: error: %e, timer: %f \n', iter, err(iter), timer(iter));
    end
end

fprintf('Maximum iterations reached, final error: %e.\n======================================\n', err(iter));
end
