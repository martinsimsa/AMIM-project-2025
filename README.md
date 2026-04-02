
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

## References
  title={A Rational Krylov Subspace Method for the Computation of the Matrix Exponential Operator}, 
  author={H. Barkouki and A. H. Bentbib and K. Jbilou},
  year={2023},
  eprint={2308.14639},
  archivePrefix={arXiv},
  primaryClass={math.NA},
  url={https://arxiv.org/abs/2308.14639}, 


V. Druskin, V. Simoncini,
Adaptive rational Krylov subspaces for large-scale dynamical systems,
Systems & Control Letters,
Volume 60, Issue 8,
2011,
Pages 546-560,
ISSN 0167-6911,
https://doi.org/10.1016/j.sysconle.2011.04.013.
(https://www.sciencedirect.com/science/article/pii/S0167691111000946)
Abstract: The rational Krylov space is recognized as a powerful tool within model order reduction techniques for linear dynamical systems. However, its success has been hindered by the lack of a parameter-free procedure, which would effectively generate the sequence of shifts used to build the space. In this paper we propose an adaptive computation of these shifts. The whole procedure only requires us to inject some initial rough estimate of the spectral region of the matrix, while further information is automatically generated during the process. The approach is a full generalization to the nonsymmetric case of the idea first proposed in Druskin et al. (2010) [18] and it is used for two important problems in control: the approximation of the transfer function and the numerical solution of large Lyapunov equations. The procedure can be naturally extended to other related problems, such as the solution of the Sylvester equation, and parametric or higher order systems. Several numerical experiments are proposed to assess the quality of the rational projection space over its most natural competitors.
Keywords: Lyapunov equation; Rational Krylov subspace; Model order reduction; Iterative methods; Transfer function

## Authors
Michaela Brzková, Jan Cirbus, Matouš Brodský, Svatopluk Vaňous, Martin Šimša

## Supervisor
Stefano Pozza, Dr., Ph.D.
