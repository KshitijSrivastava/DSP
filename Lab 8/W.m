function [ out ] = W( N,n )
%To find the twiddle function of a number
%  N is the total index present and n is the multiplication power 
%in twiddle factor
number=(2*pi*n)/N;
out=exp(  1j*number*-1);

end

