function [ H ] = H_hy( h,y )
%This function is used to find H matrix from input h(impulse response),
%y(output)
%   Detailed explanation goes here
[h_row, h_col]=size(h);%Size of h
[y_row, y_col]=size(y);%Size of y

x_col=y_col-h_col+1;
H=zeros(y_col,x_col);
k=0;%variable for shifting in the matrix
for i=1:x_col%Looping through column
    for j=1:h_col%Looping through the row
        H(j+k,i)=h(j);
    end
    k=k+1;
end

end

