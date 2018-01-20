function [conv_matrix] = conv(x,h)
[x_row, x_col]=size(x);%Size of x
[h_row, h_col]=size(h);%Size of h

size_x_col=x_col+h_col-1;%X column length

X=zeros(size_x_col,h_col);%Making X matrix full of zeros

k=0;%variable for shifting in the matrix
for i=1:h_col%Looping through column
   for j=1:x_col%Looping through the row
    X(j+k,i)=x(j);
   end
   k=k+1;
end
X;
Y = mtimes(X,transpose(h));%Matrix Multiplication
conv_matrix=transpose(Y);
end
