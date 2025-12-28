clear,clf, close all
n = 2500;
m = 30;
p = 5;

rng(1);
B = rand(n,p);
C = rand(n,p);

e = ones(n,1);
A1 = spdiags([e 2*e e], -1:1, n, n);

c = 1/2;
A2 = zeros(n,n);

for i = 1:n/2
    ai = (2*i - 1)/(n + 1);
    block = [ai, c; c, ai];
    A2(2*i-1:2*i, 2*i-1:2*i) = block;
end

lambda = eig(A2);
lambda_min = -max(lambda)*1e+2;
lambda_max = -min(lambda)*1e-2;

%lambda_min = -max(lambda)*1e+4;
%lambda_max = -min(lambda)*1e-8;

f = logm(A2 + eye(n)) / A2;

tic

[blockV,blockW, H1_0, H, sigmas] = rbla(A2,B,C,m,lambda_min,lambda_max);
eyemp = eye((m+1)*p);
E1 = eyemp(1:(m+1)*p,1:p);
fm = logm(H + eye((m+1)*p)) / H;
Fm = blockV*fm*E1*H1_0;
%Am = blockW'*A1*blockV;
error = norm(f*B - Fm,Inf)

toc