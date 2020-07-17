function [x,s] = ASAP_Hankel_1D(z,r,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x,s] = AAP_Hankel_1D(z,r,gamma)
% 
% Inputs:
% z     : Observed 1D signal. 
% r     : Target rank of underlying low rank Hankel matrix.
% gamma : Parameter for desired convergence rate. Value should between 0
%            and 1. Turn this parameter bigger will slow the convergence
%            speed but tolerate harder problem, such as higher p, r or mu.    
%
% Outputs:
% x : Estimated 1D spectrally r-sparse signal.
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
norm_of_z = norm(z);

n = size(z,1);
if mod(n,2)
    p = (n+1)/2;
    DD = [1:p p-1:-1:1]';
else
    p = n/2;
    DD = [1:p p p-1:-1:1]';
end
q = n+1-p;

c_s = max(n/p,n/q);

opts = []; opts.eta = 1e-16;

% use Cadzow to estimate incoherence and sigma_L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(n,1);
Yforward = @(y) fhmvmultiply_1D(z,y);
Ytranspose = @(y) fhmvmultiply_1D(conj(z),y);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r)); sigma_D = ss(1);
U = U(:,1:r);
V = V(:,1:r);

for i = 1:r
    ui = U(:,i);
    vi = V(:,i);
    x = x+ss(i)*conv_fft(ui,conj(vi));
end
x = x./DD;

Yforward = @(y) fhmvmultiply_1D(x,y);
Ytranspose = @(y) fhmvmultiply_1D(conj(x),y);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r)); 

sigma_L = ss(1);

U = U(:,1:r);
V = V(:,1:r);

row_norm = zeros(p,1);
for i = 1:p
    row_norm(i) = norm(U(i,:))^2;
end

col_norm = zeros(q,1);
for j = 1:q
    col_norm(j) = norm(V(j,:))^2;
end

mu = max(max(row_norm),max(col_norm))*n/(c_s*r);

beta = mu*c_s*r/(2*n);

c_init = 2; % between 1 and 3
beta_init = c_init*mu*c_s*r*sigma_L/(n*sigma_D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% two-step initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
eta = beta_init*sigma_D;
ind = find(abs(z)>eta);
s = zeros(n,1);
s(ind) = z(ind);

x = z-s;
Yforward = @(y) fhmvmultiply_1D(x,y);
Ytranspose = @(y) fhmvmultiply_1D(conj(x),y);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r));
U = U(:,1:r);
V = V(:,1:r);

x = zeros(n,1);
for i = 1:r
    ui = U(:,i);
    vi = V(:,i);
    x = x+ss(i)*conv_fft(ui,conj(vi));
end
x = x./DD;

eta = beta*ss(1);
temp = z-x;
ind = find(abs(temp)>eta);
s = zeros(n,1);
s(ind) = temp(ind);
init_timer = toc;
init_err = norm(z-x-s)/norm_of_z;
fprintf('Init: error: %e, timer: %f \n', init_err, init_timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


threshold = mu*c_s*r/n;

% iterations
for iter = 1:max_iter
    tic;
    
    % estimate the signal
    x = z-s;
    
    UtHx = zeros(r,q);
    HxV = zeros(p,r);
    Z = zeros(r,r); 
    
    for i = 1:r
        ui = U(:,i);
        UtHx(i,:) = (fhmvmultiply_1D(conj(x),ui))';
        vi = V(:,i);
        HxV(:,i) = fhmvmultiply_1D(x,vi);
    end
    
    C = UtHx*V;
    Xt = UtHx-C*(V');
    X = Xt';
    Y = HxV-U*C;
    [Q1,R1] = qr(X,0);
    [Q2,R2] = qr(Y,0);
    M = [C R1';R2 Z];
    
    % [Uc,SS,Vc] = svdecon(M);
    [Uc,SS,Vc] = svd(M,'econ');
    ss = diag(SS(1:r,1:r));
    Uc = Uc(:,1:r);
    Vc = Vc(:,1:r);
    
    U = [U Q2]*Uc;
    V = [V Q1]*Vc;
    
    x = zeros(n,1);
    for i = 1:r
        ui = U(:,i);
        vi = V(:,i);
        x = x+ss(i)*conv_fft(ui,conj(vi));
    end
    x = x./DD;
   
    

    % estimate the outlier
    eta = beta*(gamma^iter)*SS(1,1);
    s = wthresh( z - x ,'h', eta);
    
    
    err(iter) = norm (z-x-s)/norm_of_z;
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
