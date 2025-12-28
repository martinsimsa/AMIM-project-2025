% example from the paper, page 9
clear,clf
n = 1600;
p = 3;
A1 = diag(log(linspace(0.2,0.99,n)));
c = arrayfun(@(i) [(2*i-1)/(n+1), 1/2; 1/2, (2*i-1)/(n+1)], 1:n/2, 'UniformOutput', false);
A2 = blkdiag(c{:});

A = A2;

rng(3);
B = rand(n,p);
C = rand(n,p);

mineig = min(eig(A))
maxeig = max(eig(A))
condA = cond(A)

% A2 spectrum crosses zero and is mostly positive - Place all of the shifts to the right of the entire spectrum
spectrum_width = maxeig - mineig;
margin = 0.2 * spectrum_width; 
sigmamin = maxeig + margin;     
sigmamax = maxeig + 2*margin;   

% error bound precomputation
eigs = eig((A + A')/2);
mu2 = max(eigs);

% RBLA
m = 30;
%[blockV,blockW, H1_0, sigmas] = rbla_reortho(A,B,C,m, 0.0101^2, 1.61^2);
[blockV,blockW, H1_0, H, sigmas] = rbla(A,B,C,m, sigmamin,sigmamax);

%=========% Choose t
t = 1e-2;
%t = 1;
%=========%


eyemp = eye((m+1)*p);
Xt = expm(t*A)*B;

% points for bound computation
len_points = 50;
points = linspace(0,t,len_points);

% computing approximations
for k = 1:m
    % approximate solution and get errors
    E1 = eyemp(1:(k+1)*p,1:p);
    Ak = blockW(:,1:(k+1)*p)'*A*blockV(:,1:(k+1)*p);
    Xkt = blockV(:,1:(k+1)*p)*expm(t*Ak)*E1*H1_0;
    Xktprime = blockV(:,1:(k+1)*p)*Ak*expm(t*Ak)*E1*H1_0;
    residua(k) = norm(A*Xkt - Xktprime);
    errors(k) = norm(Xt - Xkt);
    
    % compute error bound
    res_norm = @(s) norm(A*blockV(:,1:(k+1)*p)*expm(s*Ak)*E1*H1_0  - blockV(:,1:(k+1)*p)*Ak*expm(s*Ak)*E1*H1_0);
    for i = 1:len_points
        res_norms(i) = res_norm(points(i));
    end
    error_bound(k) = max(res_norms)* (exp(t*mu2) - 1)/mu2;
end

% plot
semilogy(errors)
hold on
semilogy(error_bound)
legend('errors', 'error bound')
title(['Fig 3.2, approximation results, e^t^AB, t = ',num2str(t)])
xlabel('Iteration') 
ylabel('Log10 of exact error and error bound') 
hold off

