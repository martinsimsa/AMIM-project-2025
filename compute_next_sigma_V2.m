function sigma_next = compute_next_sigma_V2(lambda, sigma, eta,p)
%COMPUTE_NEXT_SIGMA_V2 Summary of this function goes here
%   Detailed explanation goes here
eta=sort(eta);
m=length(sigma);
l=length(eta);
numPoints = 10^3;  % počet bodů na interval
mu = zeros(l-1,1);
fval = zeros(l-1,1);

for i=1:l-1
    x=linspace(eta(i),eta(i+1),numPoints);
    % spočítej hodnoty r_m(x)
    rm_values = ones(size(x));
    for k = 1:m
        for j=1:p
            rm_values = rm_values .* (x - lambda(j+(k-1)*p)) ./ (x - sigma(k));
        end
    end
    % spočítej 1/|r_m(x)|
    inv_abs_rm = 1 ./ abs(rm_values);
    % najdi maximum a jeho pozici
    [fval(i), idx] = max(inv_abs_rm);
    mu(i) = x(idx);
end
% najdi celkové maximum
[~, idx_max] = max(fval);
sigma_next = mu(idx_max);
end

