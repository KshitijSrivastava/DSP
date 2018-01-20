function [ conv_ans ] = cross_corr( x,y )
y=fliplr(y);%Flipped the y 

[x_row, x_col]=size(x);%Size of x
[y_row, y_col]=size(y);%Size of y

size_x_col=x_col+y_col-1;%X column length

X=zeros(size_x_col,y_col);%Making X matrix full of zeros

k=0;%variable for shifting in the matrix
for i=1:y_col%Looping through column
   for j=1:x_col%Looping through the row
    X(j+k,i)=x(j);
   end
   k=k+1;
end
X;
Y = mtimes(X,transpose(y));%Matrix Multiplication
conv_ans=transpose(Y);

end

