%% Digital Signal Processing |[Lab-5]|                                       
% * Authors: _Kshitij Srivastava(1510110200)_ and _Nilambar Saha(partner)(1510110246)_
% * Lab Instructor: _Dr. Ravi Kant Saini_
% * Branch and Section: ECE32 Friday-1pm to 3pm
% * Email-Id: ks435@snu.edu.in
%% Objective: 
% Data Handling (In this lab we made a periodic signal from an existing 
% waveform. Also impulse response was calculated from Infinite impulse
% response for some threshold value and then finally convoluted to get the
% output)
%% Program: 
clc;
clear all;
close all;

% * |*Matlab Commands for generating the waveform by recursive equation*|
R=[1,2,3,4,5];% waveform for one period

%looping from 1 to 500 to make recursive equation for 100 cycles of size 5
for i=1:500%  looping throught the ith  position
    if i<=5
        y(i)=R(i);
    elseif mod(i,5)==0 %if modulus of i with 5 is zero
        pos=5;
        y(i)=R(pos);
    else
        pos=mod(i,5); %modulus of i with 5 
        y(i)=R(pos);
    end
end
 
%%
% * |*Matlab Commands for generating the the IIR impulse response*|
n_ref=log(0.01)/log(0.9);%Finding the n effective(Size of impulse response)
% for 1 percent of the initial value as reference
n_ref=ceil(n_ref);%Rounding off the value of n effective

for i=1:n_ref%Finding the nth term of impulse response
    h_one_per(i)=power(0.9,i-1);%0.9 to the power i-1 stored in ith index
end

%For 0.1 percent of initial value as a threshold
n_ref=log(0.001)/log(0.9);%Finding the n effective for 0.1% of intital value 
n_ref=ceil(n_ref);%Rounding off the value of impulse response

for i=1:n_ref%Finding the nth term of impulse response
    h_pointone_per(i)=power(0.9,i-1);%0.9 to the power i-1
end

% y(n)=x(n) + 0.9y(n-1)  Final equation of the impulse response

%%
% * |*Matlab Commands for convolution of h with audio signal sample*|

[audio,Fs] = audioread('Signal_Processing_Audio.mp3');
audio=audio(:,1);%Reading only the 1 column of the audio channel
audio=transpose(audio);%taking the transpose of previous ans
t=1:1/Fs:2;%Finding the sample between 0 to 2 sec at Fs frequency
[size_t_row, size_t_col]=size(t);%Finding the size of t matrix
audio=audio(1:size_t_col);%Making the audio sample for only 2 seconds

audio_out_pointone=conv(h_pointone_per,audio);%Finding the convolution 
%of audio with impulse response h_pointone_per

audio_out_one=conv( h_one_per,audio);%Finding the convolution of the audio 
%signal with impulse response h_pointone_per
%%
% * |*Matlab Commands for finding impluse response and convolution  of 
% audio with the impulse response *|

n_ref=log(0.01)/log(0.9);%Finding the n effective for 1 percent of the
% initial value
n_ref=ceil(n_ref);%Rounding off the n effective got from previous line

for i=1:n_ref %Finding the nth term of impulse response
    h_a(i)=power(-0.9,i-1);%-0.9 to the power i-1
end
% y(n) = -0.9y(n-1)+ x(n) is the recursion equation

%%
% 2*(0.9)^n=0.01  when n is even
n_ref=( log(0.01)-log(2) )/log(0.9);%Finding the n effective 
%for 1 percent of the initial value
n_ref=ceil(n_ref);%Rounding off the n effective got from previous line

for i=1:n_ref %Finding the nth term of impulse response
    h_b(i)=power(-0.9,i-1)+power(0.9,i-1);%-0.9 to the power i-1 + 0.9 
    %to the power i-1
end
% y(n)=0.81y(n-2) + 2x(n) is the recursion equation

%%
n_ref=log(0.01)/log(0.9);%Finding the n effective for 1 percent of the
% initial value
n_ref=ceil(n_ref);%Rounding off the n effective got from previous line

for i=1:n_ref %Finding the nth term of impulse response
    h_c(i)=power(0.5,i-1)+power(0.9,i-1);% 0.5 to the power i-1 + 0.9
    %to the power i-1
end
%y(n)=1.4y(n-1) -0.45y(n-2)+ 2x(n) -1.4x(n-1) is a recursion equation
%% Results:
% * |*Plot for the Periodic waveform of period 5*|

figure;plot(y);
title('Periodic waveform');xlabel('Index ');ylabel('Amplitude');

%%
% * |*Plot for the Infinite Impluse response till 1% of the initial value*|

figure;plot(h_one_per);
title('IIR till 1% of inital value');
xlabel('Index');ylabel('Amplitude');

% * |*Plot for the Infinite Impluse response till 0.1% of the initial value*|
figure;plot(h_pointone_per);
title('IIR till 0.1% of the inital value');
xlabel('Index');ylabel('Amplitude');

%%
% * |*Plot for the convolution of audio signal with IIR till 
% 1% of the initial value of impulse response*|
figure;plot(audio_out_one);
title('Convolution Output for IIR');xlabel('Index');ylabel('Amplitude');

%|*Plot for the convolution of audio signal with IIR till 
% 0.1% of the initial value of impulse response*|
figure;plot(audio_out_pointone);
title('Convolution Output for IIR');xlabel('Index');ylabel('Amplitude');

%%
%Plots for IIR response for 1% of initial value as thresholds

% * |*Plot for the Infinite Impluse response for Question 4(a)*|
figure;plot(h_a);
title('Impulse response for Ques 4(a)');xlabel('Index');ylabel('Amplitude');

% * |*Plot for the Infinite Impluse response for Question 4(b)*|
figure;plot(h_b);
title('Impulse response for Ques 4(b)');xlabel('Index');ylabel('Amplitude');

% * |*Plot for the Infinite Impluse response for Question 4(c)*|
figure;plot(h_c);
title('Impulse response for Ques 4(c)');xlabel('Index');ylabel('Amplitude');
%%
