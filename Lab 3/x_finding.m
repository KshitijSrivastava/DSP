function [ x_find ] = x_finding( H,y )
%This function is used to find input x from H matrix and y (output)
%   Detailed explanation goes here
hth=transpose(H)*H;%Finding the H*Htranspose
det(hth);%Determinant value of hth
x11=inv(hth);%Finding the inverse of hth matrix
x22=x11*transpose(H);%Multiplying the previous output by Htranspose
x_find=x22*transpose(y);%Mulatiplying the previous output by ytranspose

end

