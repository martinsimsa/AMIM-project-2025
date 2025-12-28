function sigma_next = compute_next_sigma(lambda, sigma, eta)
% lambda ... vektor lambda_j (1×m)
% sigma  ... vektor sigma_j  (1×m)
% eta    ... uzly [eta_1, ..., eta_l]
% výstup: sigma_next ... nový posun sigma_{m+1}

m = length(sigma);
l = length(eta);

% Vytvoření husté mřížky pro hledání maxima
numPoints = 200;  % počet bodů na interval
mu = zeros(l-1,1);
fval = zeros(l-1,1);

for j = 1:l-1
    % definuj interval [eta_j, eta_{j+1}]
    x = linspace(eta(j), eta(j+1), numPoints);

    % spočítej hodnoty r_m(x)
    rm_values = ones(size(x));
    for k = 1:m
        rm_values = rm_values .* (x - lambda(k)) ./ (x - sigma(k));
    end

    % spočítej 1/|r_m(x)|
    inv_abs_rm = 1 ./ abs(rm_values);

    % najdi maximum a jeho pozici
    [fval(j), idx] = max(inv_abs_rm);
    mu(j) = x(idx);
end

% najdi celkové maximum
[~, idx_max] = max(fval);
sigma_next = mu(idx_max);

end

