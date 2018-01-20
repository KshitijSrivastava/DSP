function [ X ] = X_xy( x,y )
%Finding the X matrix from x(input) and y (output)
%   Detailed explanation goes here
[x_row, x_col]=size(x);%Size of x
[y_row, y_col]=size(y);%Size of y

h_size=y_col-x_col+1;%Finding the size of impulse response
h_col=h_size;%impulse response column size

size_x_col=x_col+h_col-1;%X column length or simply y_col

X=zeros(size_x_col,h_col);%Making X matrix full of zeros

k=0;%variable for shifting in the matrix
for i=1:h_col%Looping through column
   for j=1:x_col%Looping through the row
    X(j+k,i)=x(j);
   end
   k=k+1;
end

end

