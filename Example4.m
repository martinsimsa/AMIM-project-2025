%EXAMPLE 4-using RBLA to aproximate Cauchy-Stieltjes funcion where matrix A
%comes from three point finite difference method used on two different
%differetial operators
clear ,clf
precomputed=load("example4_precomputed.mat");

%here you can choose parameters
use_precomputed=true;   %change to true/false to let code use precomputed variables/let it compute everything from scratch
n_0=100;                 %number of interval divisions in each coordinate, choose 60,80 or 100(in paper), precomputed variables for these choices
p=5;                    %number of columns of the block vectors B,C
m=70;                   %number of iterations in rbla algoritm, choose 20,30,40(in paper) or even bigger to see convergence of algorithm(I used 70),use 40 or more for table
rng(3)                  %sets seed for randomizer, for precomputed variables I used rng(3)


n=n_0^2;    %dimension of matrix
h=1/n_0;    %discretization parameter

A_1=discretize_L1_operator(n_0); %create matrix A_1 that we get from discretization of operator L_1 using finite difference method
A_2=discretize_L2_operator(n_0); %same for L_2


%using precomputed variables to save time(for n_0=100 it takes tens of minutes to compute the reference solution)
if(n_0==60&&use_precomputed)
    B=precomputed.B_60;
    C=precomputed.C_60;

    reference_solution_A_1=precomputed.reference_solution_L1_60; %see else case for explanation of variables
    reference_solution_A_2=precomputed.reference_solution_L2_60;

    sigma_min_A_1 = precomputed.sigma_min_L1_60;
    sigma_max_A_1 = precomputed.sigma_max_L1_60;
    
    sigma_min_A_2 = precomputed.sigma_min_L2_60;
    sigma_max_A_2 = precomputed.sigma_max_L2_60;

elseif(n_0==80&&use_precomputed)
    B=precomputed.B_80;
    C=precomputed.C_80;

    reference_solution_A_1=precomputed.reference_solution_L1_80;
    reference_solution_A_2=precomputed.reference_solution_L2_80;

    sigma_min_A_1 = precomputed.sigma_min_L1_80;
    sigma_max_A_1 = precomputed.sigma_max_L1_80;
    
    sigma_min_A_2 = precomputed.sigma_min_L2_80;
    sigma_max_A_2 = precomputed.sigma_max_L2_80;

elseif(n_0==100&&use_precomputed)
    B=precomputed.B_100;
    C=precomputed.C_100;

    reference_solution_A_1=precomputed.reference_solution_L1_100;
    reference_solution_A_2=precomputed.reference_solution_L2_100;

    sigma_min_A_1 = precomputed.sigma_min_L1_100;
    sigma_max_A_1 = precomputed.sigma_max_L1_100;
    
    sigma_min_A_2 = precomputed.sigma_min_L2_100;
    sigma_max_A_2 = precomputed.sigma_max_L2_100;

else
    B=rand(n,p);    %create random block vectors B and C   
    C=rand(n,p);

    %Computing spectrum of matrix A_1
    t1_start=tic;
    A_1=full(A_1);
    spectrum_A_1=eig(A_1);
    t1_end=toc(t1_start);
    fprintf("Computitation time of spectrum for A1: %.2f seconds\n",t1_end)
    %figure()                                           %optional to see graphed spectrum
    %plot(real(spectrum_A_1),imag(spectrum_A_1),"b*")
    A_1=sparse(A_1);
        
    %Computing spectrum of matrix A_2
    t2_start=tic;
    A_2=full(A_2);
    spectrum_A_2=eig(A_2);
    t2_end=toc(t2_start);
    fprintf("Computitation time of spectrum for A2: %.2f seconds\n",t2_end)
    %figure()
    %plot(real(spectrum_A_2),imag(spectrum_A_2),"b*")
    A_2=sparse(A_2);
        
    %Computation of reference solution (A_1)^(-1/2)*B we take as exact
    t5_start=tic;
    reference_solution_A_1=(A_1)^(-1/2)*B;
    t5_end=toc(t5_start);
    fprintf("Computation time of reference solution for A1: %.2f seconds\n",t5_end)
        
    %Computation of reference solution (A_2)^(-1/2)*B we take as exact
    t4_start=tic;
    reference_solution_A_2=(A_2)^(-1/2)*B;
    t4_end=toc(t4_start);
    fprintf("Computation time of reference solution for A2: %.2f seconds\n",t4_end)
        
    %Choosing initial interval for shifts
    sigma_min_A_1 = -max(spectrum_A_1);
    sigma_max_A_1 = -min(spectrum_A_1);
        
    sigma_min_A_2 = -max(spectrum_A_2);
    sigma_max_A_2 = -min(spectrum_A_2);
end

%RBLA for A_1
t5_start=tic;
[blockV_A_1,blockW_A_1, H1_0_A_1, H_A_1]=rbla(A_1,B,C,m,sigma_min_A_1,sigma_max_A_1); 
t5_end=toc(t5_start);
fprintf("Computation time of the RBLA algorithm for A1: %.2f seconds\n",t5_end)

%RBLA for A_2
t6_start=tic;
[blockV_A_2,blockW_A_2, H1_0_A_2, H_A_2]=rbla(A_2,B,C,m,sigma_min_A_2,sigma_max_A_2);
t6_end=toc(t6_start);
fprintf("Computation time of the RBLA algorithm for A2: %.2f seconds\n",t6_end)

%Computation of estimate using RBLA outputs of (A_1)^(-1/2)*B and (A_2)^(-1/2)*B
exact_error_vector_A_1=zeros(m,1);
exact_error_vector_A_2=zeros(m,1);
eyemp = eye((m+1)*p);
for k = 1:m
    E1 = eyemp(1:(k+1)*p,1:p);
    Ak_A_1 = blockW_A_1(:,1:(k+1)*p)'*A_1*blockV_A_1(:,1:(k+1)*p);
    Ak_A_2 = blockW_A_2(:,1:(k+1)*p)'*A_2*blockV_A_2(:,1:(k+1)*p);
    fk_A_1 = blockV_A_1(:,1:(k+1)*p)*Ak_A_1^(-1/2)*E1*H1_0_A_1;
    fk_A_2 = blockV_A_2(:,1:(k+1)*p)*Ak_A_2^(-1/2)*E1*H1_0_A_2;
    exact_error_vector_A_1(k)=norm(reference_solution_A_1-fk_A_1,"inf");
    exact_error_vector_A_2(k)=norm(reference_solution_A_2-fk_A_2,"inf");
end

%displays convergence of approximate solution to reference one for L1
figure()
semilogy(exact_error_vector_A_1,"--xr")
title(sprintf("Exact error of computed approximation for L1 operator with n_{0}=%d",n_0))
xlabel("Number of iterations")
ylabel("Error")

%displays convergence of approximate solution to reference one for L2
figure()
semilogy(exact_error_vector_A_2,"--xr")
title(sprintf("Exact error of computed approximation for L2 operator with n_{0}=%d",n_0))
xlabel("Number of iterations")
ylabel("Error")

fprintf("Exact error at last iteration: %d  (m=%d)\n",norm(reference_solution_A_1-fk_A_1,"inf"),m)
fprintf("Exact error at last iteration: %d  (m=%d)\n",norm(reference_solution_A_2-fk_A_2,"inf"),m)

if m>=40 %displays table with wanted quantities
    fprintf("\n")
    fprintf("---------------------------------\n")
    fprintf("Problem| Dim    | Exact error   |\n")
    fprintf("---------------------------------\n")
    fprintf("L1(u)  | m=20   | %.3e     |\n",exact_error_vector_A_1(20))
    fprintf("       | m=30   | %.3e     |\n",exact_error_vector_A_1(30))
    fprintf("       | m=40   | %.3e     |\n",exact_error_vector_A_1(40))
    fprintf("---------------------------------\n")
    fprintf("L2(u)  | m=20   | %.3e     |\n",exact_error_vector_A_2(20))
    fprintf("       | m=30   | %.3e     |\n",exact_error_vector_A_2(30))
    fprintf("       | m=40   | %.3e     |\n",exact_error_vector_A_2(40))
    fprintf("---------------------------------\n")
end
fprintf("\n")