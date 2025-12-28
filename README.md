
# Project for the article BarBenJbi23
A RATIONAL KRYLOV SUBSPACE METHOD FOR THE COMPUTATION OF THE MATRIX EXPONENTIAL OPERATOR

## Contains
Project contain the following files:
- lib - matrices for Example 2 and 3
    - add32.mat
    - rail_1357.mat
    - rail_5177.mat
    - rail_20209.mat
- rbla.m - main algorithm to compute rational block lanczos algorithm
- compute_next_sigma.m - function computing shifts in each iteration
- Numerical_example_A1.m - page 9 choose t = 1e-2 or t = 1 on line 26
- Numerical_example_A2.m - page 9 choose t = 1e-2 or t = 1 on line 34
- Example1.m - page 12, choose testCase from {'poisson', 'fdm'} on line 4 and m = 10 or 15 on line 7
- Example2.m - page 13-14, choose problem from {'poisson', 'add32'} on line 9
- Example3.m - page 14, choose problem from {'add32', 'rail1357', 'rail5177', 'rail20209'} on line 9
- Example4.m - page 14-15
- Example5.m - page 16
- Example6.m - page 16
- discretize_L1_operator.m - function used in Example4
- discretize_L2_operator.m - function used in Example4
- example4_precomputed.mat - precompution for Example4
- README.md

## Authors
Michaela Brzková, Jan Cirbus, Matouš Brodský, Svatopluk Vaňous, Martin Šimša

## Supervisor
Stefano Pozza, Dr., Ph.D.