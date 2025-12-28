% Example 3
clear,clf
close all

% list of parameters t in which function is evaluated:
ts = [1e-2, 5*1e-2, 1e-1, 5*1e-1, 1, 10];

% choose problem: {'add32', 'rail1357', 'rail5177', 'rail20209'}
problem = 'add32';
if strcmp(problem, 'add32')
    fprintf('\n add32 matrix\n')
    n = 4960;
    p = 4;
    Problem = ssget(540);
    A = Problem.A;
    B = rand(n,p);
    C = B;
    m = 10;
    % manually set sigmas work better
    %spectrum = eigs(A,100);
    %sigmamin = -max(spectrum)*1.2;
    %sigmamax = -min(spectrum)*0.8;
    sigmamin = -0.7;
    sigmamax = -0.4;
elseif strcmp(problem, 'rail1357')
    fprintf('\n rail-1357 matrix\n')
    n = 1357;
    Problem = ssget(1442);
    A = Problem.A;
    p = 6;
    B = full(Problem.aux.B(:,1:p));
    C = B;
    m = 10;
    %spectrum = eigs(A,100);
    %sigmamin = -max(spectrum)*1.2;
    %sigmamax = -min(spectrum)*0.8;
    sigmamin = 1e-6;
    sigmamax = 1e-2;
elseif strcmp(problem, 'rail5177')
    fprintf('\n rail-5177 matrix\n')
    n = 5177;
    Problem = ssget(1444);
    A = Problem.A;
    p = 6;
    B = full(Problem.aux.B(:,1:p));
    C = B;
    m = 20;
    %spectrum = eigs(A,100);
    %sigmamin = -max(spectrum)*1.2;
    %sigmamax = -min(spectrum)*0.8;
    sigmamin = 1e-6;
    sigmamax = 1e-2;
elseif strcmp(problem, 'rail20209')
    fprintf('\n rail-20209\n')
    n = 20209;
    p = 6;
    Problem = ssget(1443);
    A = Problem.A;
    B = full(Problem.aux.B(:,1:6));
    C = B;
    m = 25;
    sigmamin = 1e-6;
    sigmamax = 1e-2;
end


[blockV, blockW, H1_0, H, sigmas, blockK] = rbla(A, B, C, m, sigmamin, sigmamax);
Vm = blockV(:, 1:m*p);
Wm = blockW(:, 1:m*p);
eyemp = eye(m*p);
E1 = eyemp(1:m*p,1:p);
Am = H(1:m*p, 1:m*p);

% for theorem 3.2 residuum
Hmp1_m = sigmas(m+1)*blockK(m*p+1:(m+1)*p,(m-1)*p+1:m*p);
Vmp1 = blockV(:, m*p + 1:(m+1)*p);
Km = blockK(1:m*p,:);
BlockEm = eye(m*p, m*p);
Em = BlockEm(:,(m-1)*p + 1:end);
for k = 1:6
    t = ts(k);
    Xmt = Vm*expm(t*Am)*E1*H1_0;
    Xmtprime = Vm*Am*expm(t*Am)*E1*H1_0;
    Rk = A*Xmt - Xmtprime;
    residua(k) = norm(Rk, inf);
    % Theorem 3.2 residuum
    temp = Km\expm(t*Am)*E1*H1_0;
    Rk2 = Vmp1*Hmp1_m*Em'*temp;
    residua2(k) = norm(Rk2,inf);
    fprintf('\nt = %d\n', t)
    fprintf('Norm of residuum: %d\n', norm(Rk, Inf))
    fprintf('Theorem 3.2 residuum for t = %g: %g\n', t, residua2(k));
end
tic
if strcmp(problem, 'rail20209')
    Xt = zeros(size(B));
    t = ts(6);
    koef = 1e4;
    for j = 1:size(B,2)
       Xt(:,j) = expmv(A, B(:,j),t);
    end
else
    Xt = expm(ts(6)*A)*B;
end
time = toc;
fprintf('time to evaluate six times expm(t*A)*B: %d seconds.\n', 6*time)


