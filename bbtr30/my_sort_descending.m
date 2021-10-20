function OutputMatrix = my_sort_descending(Matrix,n,m,c)
%
% OutputMatrix = my_sort_descending(Matrix,n,m,column)
%
% This function sorts in descending order a n*m matrix by a given column 
%
% Inputs & output of this function are:
%
% Matrix: the matrix to be sorted
% n: number of rows of matrix
% m: number of columns of matrix
% c: the column used as reference for the sort
%
% OutputMatrix: the sorted matrix
%

for i = 1 : n-1
    for j = i+1 : n
        
        if Matrix(j,c) > Matrix(i,c)
            temp = Matrix(j,1:m);
            Matrix(j,1:m) = Matrix(i,1:m);            
            Matrix(i,1:m) = temp;
        end
        
    end
end

OutputMatrix = Matrix;

return