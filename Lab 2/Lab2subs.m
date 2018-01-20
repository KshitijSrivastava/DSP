%% Digital Signal Processing |[Lab-2]|                                       
% * Authors: _Kshitij Srivastava(1510110200)_ and _Nilambar Saha(Partner)(1510110246)_
% * Lab Instructor: _Dr. Ravi Kant Saini_   
%% Objective: 
% Convolution and Block Convolution (In this lab, we did the simple convolution
% and for a bigger input block convolution was done so input was divided 
% into smaller blocks and convlution was individually done)
%% Program: 
clc;
clear all;
close all;

% * |*Matlab Commands for Convolution*|
x=input('enter the x');%Enter the value of x
[x_row, x_col]=size(x);%Size of x
h=input('enter the h');%Enter the value of h
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
conv_ans=transpose(Y)

%%
% * |*Matlab Commands for Block Convolution fo impulse response h1*|
[y, fs]=audioread('Signal_Processing_Audio.mp3');
y_n=y(:,1);
t=0:1/fs:5;% 5 sec of samples at fs sampling
size_y = size(y_n);
[size_t_row, size_t_col] = size(t);%Size of time samples
y_fivesec=y_n(1:size_t_col);%Taking 5 sec of audio samples

zeros_to_add=mod( size(y_fivesec) , 512 );
y_fivesec = vertcat(y_fivesec,zeros(171,1));%171 zeros added to the input
no_interations=size(y_fivesec)/512;%no_iterations=431

vector_y=transpose(y_fivesec);%Making to horizontal matrix
isvector(vector_y);

h=h1;
input_matrix=zeros(431,512);%Input matrix for making input to groups of 512
k=1;
for i=1:431%Iterating over 431 rows
    input=vector_y(k:k+511);
    input_matrix(i,:)=input;%Addign to the ith row
    k=k+512;
end

off_transient=zeros(1,60);%Initializing off transients
on_transient=zeros(1,60);%Initializing on transients
for i=1:431
conv_matrix=conv(input_matrix(i,:),h);%Selecting groups of 512 elements

if i==1%if first group of elements then add on-transients to the output
output1=conv_matrix(1:60);
else
on_transient=conv_matrix(1:60);%Finding on transients
output1=horzcat(output1,off_transient+on_transient);%Adding off and on transients
end

%Adding elements which are neither off nor on transients
output1=horzcat(output1,conv_matrix(61:512));

if i==431%if last group of elements then add off-transients to the output
output1=horzcat(output1,conv_matrix(513:572));
else
off_transient=conv_matrix(513:572);%Finding the off transients        
end  

end

%%
% * |*Matlab Commands for Block Convolution for h2*|
[y, fs]=audioread('Signal_Processing_Audio.mp3');
y_n=y(:,1);
t=0:1/fs:5;%Time period for taking 5sec of samples
size_y = size(y_n);
[size_t_row, size_t_col] = size(t);
y_fivesec=y_n(1:size_t_col);%Taking 5 sec of audio samples

zeros_to_add=mod( size(y_fivesec) , 512 );
y_fivesec = vertcat(y_fivesec,zeros(171,1));%171 zeros added to the input
no_interations=size(y_fivesec)/512;%no_iterations=431

vector_y=transpose(y_fivesec);%Making to horizontal matrix
isvector(vector_y);

h=h2;
input_matrix=zeros(431,512);%Input matrix for making input to groups of 512
k=1;
for i=1:431%Iterating over 431 rows
    input=vector_y(k:k+511);
    input_matrix(i,:)=input;%Addign to the ith row
    k=k+512;
end

off_transient=zeros(1,60);%Initializing off transients
on_transient=zeros(1,60);%Initializing on transients
for i=1:431
conv_matrix=conv(input_matrix(i,:),h);%Selecting groups of 512 elements

if i==1%if first group of elements then add on-transients to the output
output=conv_matrix(1:60);
else
on_transient=conv_matrix(1:60);%Finding on transients
output=horzcat(output,off_transient+on_transient);%Adding off and on transients
end

%Adding elements which are neither off nor on transients
output=horzcat(output,conv_matrix(61:512));

if i==431%if last group of elements then add off-transients to the output
output=horzcat(output,conv_matrix(513:572));
else
off_transient=conv_matrix(513:572);%Finding the off transients        
end  

end

%% Results:
% * |*Results for Convlution*|
conv_ans
%%
% * |*Results for Block Convlution for impulse response h1*|
output1
%%
% * |*Results for Block Convlution for impulse response h2*|
output
