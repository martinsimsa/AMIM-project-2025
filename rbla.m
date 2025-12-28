function [blockV, blockW, H1_0, H, sigmas, blockK] = rbla(A,B,C,m, sigmamin, sigmamax)
% RBLA with re-orthogonalization and robust SVD-normalization.

% ---- init ----
n = size(A,1);
p = size(B,2);

S0 = B;
R0 = C;

% initial biorthogonalization
[delta,beta] = qr(C'*B);
Vk = S0 / beta;
Wk = R0 * delta;
H1_0 = beta;

blockV = Vk;
blockW = Wk;

% --- CHANGE 1: allocate sigmas as vector (m+1) ---
sigmas = zeros(m+1,1);
sigmas(1) = sigmamin;

% --- CHANGE 2: etas must be sorted nodes for interval scanning ---
etas = sort([sigmamin, sigmamax]);  % start with two endpoints

% pre-allocate H small matrix (we'll fill at end by projection)
H = zeros((m+1)*p, (m+1)*p);

% threshold for tiny singular values (avoid divide-by-0)
tol_svd_factor = 1e3;
min_abs_tol = 1e-14;

% initialize block K matrix for computing residua
blockK = zeros((m+1)*p,m*p);
Em = eye((m+1)*p, m*p);

for k = 1:m
    % build small projected Ak for shift selection
    Ak = blockW(:,1:k*p)' * A * blockV(:,1:k*p);
    lambdas = eig(Ak);

    sigma_next = compute_next_sigma(lambdas, sigmas(1:k), etas);
    
    % --- CHANGE 3: stabilizační clamp ---
    sigma_next = min(max(sigma_next, sigmamin), sigmamax);
    if abs(sigma_next) < 1e-6
        sigma_next = sign(sigmamin)*1e-6;
    end
    
    sigma_next = sigma_next + sign(sigma_next)*1e-6;
    sigmas(k+1) = sigma_next;
    
    etas = sort(unique([etas, sigma_next]));

    % Solve shifted systems robustly
    sigma_k = sigmas(k+1);
    M = speye(n) - A / sigma_k;
    if condest(M) < 1e-14
        warning('rbla_reortho_fixed: M ill-conditioned at k=%d, perturbing sigma', k);
        sigma_k = sigma_k + sign(sigma_k)*1e-3;
        sigmas(k+1) = sigma_k;
        M = speye(n) - A / sigma_k;
    end

    Sk = M \ (A * Vk);
    Rk = M' \ (A' * Wk);

    % --------- biorthogonalization (1st pass) ---------
    Hk = blockW' * Sk;
    Gk = blockV' * Rk;
    Sk = Sk - blockV * Hk;
    Rk = Rk - blockW * Gk;

    % --------- re-biorthogonalization (2nd pass) ------
    Hk2 = blockW' * Sk;
    Gk2 = blockV' * Rk;
    Sk  = Sk - blockV * Hk2;
    Rk  = Rk - blockW * Gk2;

    % --------- orthonormalize new blocks robustly ----
    [Q, R_q] = qr(Sk, 0);
    [Q2, R_q2] = qr(Q, 0);
    Vk = Q2;
    R_total_V = R_q2 * R_q;

    [Qw, R_w] = qr(Rk, 0);
    [Qw2, R_w2] = qr(Qw, 0);
    Wk = Qw2;
    R_total_W = R_w2 * R_w;

    % --------- biorthogonal normalization via SVD ----
    Msmall = Wk' * Vk;
    [Pk, Dk, Qk] = svd(Msmall);
    d = diag(Dk);
    dmax = max(d);
    tolSmall = max(eps(dmax) * tol_svd_factor, min_abs_tol);
    d(d < tolSmall) = tolSmall;
    sqrtD = diag(sqrt(d));
    invSqrtD = diag(1./sqrt(d));

    Vk = Vk * Qk * invSqrtD;
    Wk = Wk * Pk * invSqrtD;

    % keep these (not used later, but preserve your computations)
    Hkp1_k = sqrtD * Qk' * R_total_V; %#ok<NASGU>
    Gkp1_k = sqrtD * Pk' * R_total_W; %#ok<NASGU>

    % append new blocks
    blockV = [blockV, Vk];
    blockW = [blockW, Wk];

    % Construct block K for residuum computation
    Hkboth = Hk + Hk2;
    Kkbar = 1/sigma_k*[Hkboth;Hkp1_k];
    blockK(1:(k+1)*p,(k-1)*p+1:k*p) = Kkbar;
end

% construct small projected matrix H
H = blockW' * (A * blockV);

% add identity part to the blockK
blockK = blockK + Em;

end
