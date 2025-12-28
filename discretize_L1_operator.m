function [A] = discretize_L1_operator(number_of_partitions)
%Creates matrix that coresponds to discretizing L1 operator using finite
%difference method on square grid with number of nodes taken as parameter
%in each coordinal direction

n_0=number_of_partitions;
dim_of_A=number_of_partitions^2;

%diagonal and superdiagonals of A
main_diagonal=4*ones(dim_of_A,1);
superdiagonals_of_A=-1*ones(dim_of_A-1,1);
for i=1:dim_of_A-1
    if mod(i,n_0)==0
        superdiagonals_of_A(i)=0;
    end
end
%subdiagonal coresponding to identity matrices in structure of A
subdiagonals_of_A=-1*ones(dim_of_A-n_0,1);

A=sparse(diag(main_diagonal)+diag(superdiagonals_of_A,1)+diag(superdiagonals_of_A,-1)+diag(subdiagonals_of_A,n_0)+diag(subdiagonals_of_A,-(n_0)));
A=dim_of_A*A;
end

