function [A] = discretize_L2_operator(number_of_partitions)
%Creates matrix that coresponds to discretizing L2 operator using finite
%difference method on square grid with number of nodes taken as parameter
%in each coordinal direction


n_0=number_of_partitions;
dim_of_A=number_of_partitions^2;

%diagonal and superdiagonals of A
main_diagonal=202*dim_of_A*ones(dim_of_A,1);
superdiagonals_of_A=-100*dim_of_A*ones(dim_of_A-1,1);
%subdiagonal coresponding to identity matrices in structure of A
subdiagonals_of_A=-1*dim_of_A*ones(dim_of_A-n_0,1);

%addition to superdiagonals corresponding to discretization of part containing x in operator
upper_super_diagonal=superdiagonals_of_A;
lower_super_diagonal=superdiagonals_of_A;
for i=1:(dim_of_A-1)
    if i<n_0
        upper_super_diagonal(i)=upper_super_diagonal(i)+5*i;
        lower_super_diagonal(i)=lower_super_diagonal(i)-5*(i+1);
    elseif mod(i,n_0)==0
        upper_super_diagonal(i)=0;
        lower_super_diagonal(i)=0;
    else
        upper_super_diagonal(i)=upper_super_diagonal(mod(i,n_0));
        lower_super_diagonal(i)=lower_super_diagonal(mod(i,n_0));
    end
end

A=sparse(diag(main_diagonal)+diag(upper_super_diagonal,1)+diag(lower_super_diagonal,-1)+diag(subdiagonals_of_A,n_0)+diag(subdiagonals_of_A,-(n_0)));
end


