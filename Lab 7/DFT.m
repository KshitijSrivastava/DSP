function [ out ] = DFT( x )
%This function is used to find DFT of any sequence

N=length(x);
 D_mat_har=zeros(N,N);
 D_mat=zeros(N,N);
%Finding the D matrix harmesian for size N
for i=1:N 
    for j=1:N 
        number=((2*pi*(i-1)*(j-1))/N);
        D_mat(i,j)=exp(  1j*number );
        D_mat_har(i,j)=conj(D_mat(i,j));%Finding the conjugate
       
    end
end

out=D_mat_har*transpose(x);%Finding the DFT

end

