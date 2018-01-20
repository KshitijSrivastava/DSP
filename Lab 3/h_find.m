function [ h_find ] = h_find( X,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
xtx=transpose(X)*X;
det(xtx)
h1=inv(xtx);
h2=h1*transpose(X);
h_find=h2*transpose(y);
end

