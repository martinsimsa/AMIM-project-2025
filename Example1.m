clear; clc; close all;

% choose testCase = {'poisson, 'fdm'}
testCase = "poisson";
n0 = 40;
p  = 3;
m  = 10; % m = 10 or 15
h = 1/(n0+1);
tvals = logspace(-1, 0, 50);
rngSeed = 1;
makeA_sparse = true;

switch testCase
    case "fdm"
        A = build_fdm_matrix(n0);
        sigmamin = -8/h^2;
        sigmamax = -1e-2;
    case "poisson"
        A = build_poisson_matrix(n0);
        sigmamin = -4/h^2; 
        sigmamax = -1;
    otherwise
        error('Unknown testCase.');
end

if makeA_sparse
    A = sparse(A);
end

n = size(A,1);
rng(rngSeed);
B = rand(n,p);
C = rand(n,p);

[blockV, blockW, H1_0, Hfull, sigmas] = rbla(A, B, C, m, sigmamin, sigmamax);

% Use V_m and W_m 
mp = m*p;
Vm = blockV(:, 1:mp);
Wm = blockW(:, 1:mp);

% Projected matrix A_m = W_m^T A V_m
Am = Hfull(1:mp, 1:mp);

% E1 selector (first p columns of I_{mp})
E1 = [eye(p); zeros(mp-p,p)];



for it = 1:numel(tvals)
    t = tvals(it);

    Xm = Vm * (expm(t*Am) * (E1 * H1_0));
    Xexact = expm(t*A) * B;

    errVals(it) = log10(norm(Xexact - Xm, inf));
end




figure;
plot(tvals, errVals, '-*r');
grid on;
xlabel('time t');
ylabel('log_{10} of the error');
title('Approximation of exp(tA)B');



% Building matrices
function A = build_fdm_matrix(n0)
h = 1/(n0+1);
x = (1:n0)*h;
y = (1:n0)*h;

n = n0^2;
A = spalloc(n,n, 5*n);
idx = @(i,j) (j-1)*n0 + i;

for j = 1:n0
    for i = 1:n0
        k = idx(i,j);
        xi = x(i); yj = y(j);

        f = exp(xi*yj);
        g = sin(xi*yj);
        hh = yj^2 - xi^2;
        diagVal = -4/(h^2) - hh;
        if i < n0
            kp = idx(i+1,j);
            A(k,kp) = A(k,kp) + 1/(h^2) - f*(1/(2*h));
        end
        if i > 1
            km = idx(i-1,j);
            A(k,km) = A(k,km) + 1/(h^2) + f*(1/(2*h));
        end
        if j < n0
            kp = idx(i,j+1);
            A(k,kp) = A(k,kp) + 1/(h^2) - g*(1/(2*h));
        end
        if j > 1
            km = idx(i,j-1);
            A(k,km) = A(k,km) + 1/(h^2) + g*(1/(2*h));
        end

        A(k,k) = A(k,k) + diagVal;
    end
end
end

function A = build_poisson_matrix(n0)
h = 1/(n0+1);
n = n0^2;
A = spalloc(n,n, 5*n);

idx = @(i,j) (j-1)*n0 + i;

for j = 1:n0
    for i = 1:n0
        k = idx(i,j);
        A(k,k) = -4/(h^2);
        if i < n0, A(k, idx(i+1,j)) = 1/(h^2); end
        if i > 1,  A(k, idx(i-1,j)) = 1/(h^2); end
        if j < n0, A(k, idx(i,j+1)) = 1/(h^2); end
        if j > 1,  A(k, idx(i,j-1)) = 1/(h^2); end
    end
end
end

