% Example 2
clear,clc, close all
fprintf('Example 2\n')

% Random seed for reproducibility
rng(1);

% choose problem: {'poisson', 'add32'}
problem = 'poisson';
if strcmp(problem, 'poisson')
    p = 3;
    fprintf('Poisson matrix, p = %d\n', p)
    n0 = 80;          % so n = 80^2 = 6400
    n = 6400;
    e = ones(n0,1);
    T = spdiags([e -2*e e], -1:1, n0, n0);
    I = speye(n0);  
    A = kron(I, T) + kron(T, I);
    h = 1/(n0+1);
    A = A/(h^2);
    B = rand(n,p);
    C = B;
    ts = [1,2];
    sigmamin = 1;
    sigmamax = 4/h^2;
elseif strcmp(problem, 'add32')
    n = 4960;
    %Problem = ssget(540);
    load("lib/add32.mat");
    A = Problem.A;
    p = 4;
    fprintf('add32 matrix, p = %d\n', p)
    rng(0);
    B = rand(n,p);
    C = B;
    %C = rand(n,p);
    ts = [1e-2, 1e-1];
    sigmamin = -5;
    sigmamax = -0.01;
else
    fprintf('Invalid problem selection.\n');
    return;
end

% RBLA
m = 40;
[blockV, blockW, H1_0, H, sigmas] = rbla(A, B, C, m, sigmamin, sigmamax);


eyemp = eye(m*p);
Am = H(1:m*p, 1:m*p);
dims = [10,20,30,40];

for j = 1:2
    t = ts(j);
    fprintf('t = %d \n', t)
    Xt = expm(t*A)*B;
    for i = 1:4
        k = dims(i);
        E1 = eyemp(1:k*p,1:p);
        Ak = Am(1:k*p,1:k*p);
        Xkt = blockV(:,1:k*p)*expm(t*Ak)*E1(1:k*p,:)*H1_0;
        fprintf('dim = %d, true error: %d \n', k, norm(Xt - Xkt, inf))
    end
end
