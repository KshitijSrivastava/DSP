function [ h_find ] = h_finding( X,y )
%Finding the h(impulse response) from output y and X matrix
%   Detailed explanation goes here
xtx=transpose(X)*X;%Multiplication of X into Xtranspose
det(xtx);%finding the determinant value of xtx
h11=inv(xtx);%finding the inverse of xtx matrix
h22=h11*transpose(X);%Multiplying the previous output by Xtranspose 
h_find=h22*transpose(y);%Multiplying the previous output by ytranspose 
end

