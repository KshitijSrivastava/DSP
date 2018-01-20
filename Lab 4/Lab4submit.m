%% Digital Signal Processing |[Lab-4]|                                       
% * Authors: _Kshitij Srivastava(1510110200)_ and _Nilambar Saha(partner)(1510110246)_
% * Lab Instructor: _Dr. Ravi Kant Saini_ 
% * Branch and Section: ECE32 Friday-1pm to 3pm
% * Email-Id: ks435@snu.edu.in
%% Objective: 
% Correlation (In this lab we computed correlation of signals. 
% we also computed the timeperiod of noisy period signal by computing the
% the auto-correlation of the signal and then finding period with its peaks) 
%% Program: 
clc;
clear all;
close all;

% * |*Matlab Commands Cross-Correlation of class example*|

x=[1,2,3,4];
y=[4,3,2,1];
out=zeros(1,4);%For reversing the y
for i=1:4
    out(i)=y(5-i);%Changing the poition of y and storing in out varible
end
y=out;%y is out variable (Reversing is done)

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
Y = mtimes(X,transpose(y));%Matrix Multiplication
conv_ans=transpose(Y);

%%
% * |*Matlab Commands for Cross-Correlation using function*|

% cross_corr function stored in cross_corr.m file

% function [ conv_ans ] = cross_corr( x,y )
% This function computes cross_correlation of input x,y

% y=fliplr(y);%Flipped the y  
% [x_row, x_col]=size(x);%Size of x
% [y_row, y_col]=size(y);%Size of y
% size_x_col=x_col+y_col-1;%X column length
% X=zeros(size_x_col,y_col);%Making X matrix full of zeros
% 
% k=0;%variable for shifting in the matrix
% for i=1:y_col%Looping through column
%    for j=1:x_col%Looping through the row
%     X(j+k,i)=x(j);
%    end
%    k=k+1;
% end
% Y = mtimes(X,transpose(y));%Matrix Multiplication
% conv_ans=transpose(Y);
%
% end

x=[1,2,3,4];
y=[4,3,2,1];

func_crossc=cross_corr(x,y);

%%
% * |*Matlab Commands for finding period in noise added periodic signal*|

load('noiseData.mat');
noise_autoc=cross_corr(noiseData,noiseData);%Computing cross-correlation
count=0;%Counter for counting the peaks of cross correlation
for i=2:98
    if noise_autoc(i)>noise_autoc(i-1)&&noise_autoc(i)>noise_autoc(i+1)
        %if the ith sample is more than its immediate neighbours
    count=count+1;%Increase the counter
    peak(count)=noise_autoc(i);%Storing the peak value in the peak vaiable
    index(count)=i;%Storing the peak value indexes
    end
end

[size_index_row, size_index_col]=size(index);%Finding the size of index
for i=1:size_index_col-1
    diff_index(i)=index(i+1)-index(i);%Taking the difference between indexes
end
[diff_index_row, diff_index_col]=size(diff_index);%Finding the size of 
% diff_index variable
avg_ans=sum(diff_index)/diff_index_col;%Finding the average difference in
%the peaks by computing the sum and divding by size of diff_index

%% Results:
% * |*Result for Cross-Correlation of class example *|
conv_ans

%%
% * |*Result for Cross-Correlation using function*|
func_crossc

%%
% * |*Result for finding period in noise added periodic signal*|
avg_ans

%%