%% Digital Signal Processing |[Lab-9]|                                       
% * Authors: _Kshitij Srivastava(1510110200)_ and _Nilambar Saha(partner)(1510110246)_
% * Lab Instructor: _Dr. Ravi Kant Saini_
% * Batch: ECE32, Friday 1 to 3 PM
% * Email: ks435@snu.edu.in
%% Objective: 
% Design of window based linear phase FIR filter  (In this we found the 
% freqency response of various kinds of window based linear phase
% FIR filters )
%% Program: 
clc;
clear all;
close all;

% * |*Matlab Commands for rectangular windowed ideal low pass filter *|

wc=pi/4; % Variable for cut-off frequency
h1=zeros(1,61); %Initializing the impulse response h1
for i=-30:30 %Looping through all the 61 elements of the impulse response
    if i==0
        h1(i+31)=wc/pi;
    else
        h1(i+31)=sin(wc*i)/(pi*i);
    end
end

%%
% * |*Matlab Commands to plot the 512 point frequency response *|
freq=zeros(1,1024); %For finding the 512 point frequency response
diff=(2*pi)/1024; %Finding the difference between 2 frequency from -pi to pi
%with 512 divisions
for i=1:1024 %Looping through the 512 impulse response
    w=(-1*pi)+((i-1)*diff); %Finding the angular frequency
    sum=0; %Variable  sum=0
    for k=0:(60/2) %Looping through the 30 elements index
        index=k+1;
        if k==0
            sum=sum+ ( h1(31)* cos(0*w) );
        else
            sum=sum+ ( 2*h1(30-k+1)* cos(k*w) );
        end
    end
    freq(i)=sum* exp((-1*1i*60*w)/2); %Finding the frequency response 
    %of each term
end
freq=abs(freq); %Finding the absolute value for 512 point freq response

%%
% * |*Matlab Commands for frequency response by using the built in function *|
rec_lp=freqz(h1);

%%
% * |*Matlab Commands for Rectangular window  *|
N=61; %Variable N to store the size N

for i=0:60 %Looping through 61 elements of the rectangular window
    if i>=0 && i<=60  %If 0<=n<=N-1, then
        w1(i+1)=1;
    else              %else 
        w1(i+1)=0;
    end
end
freq_rec=freqz(w1); %Finding the frequency response of rectangular window


%%
% * |*Matlab Commands for Hamming Window  *|
alpha=0.54; %Alpha value for generalized Hamming window
for i=0:60 %Looping through the 61 elements of the frequency response
    w_ham(i+1)=alpha- ( (1-alpha)*cos((2*pi*i)/60) )  ;
end
freq_hamm=freqz(w_ham); %Finding the frequency response of Hamming window

%%
% * |*Matlab Commands for Hanning window  *|
alpha=0.5;%Alpha value for generalized Hanning window
for i=0:60  %Looping through the 61 elements of the frequency response
    w_han(i+1)=alpha-( (1-alpha)*cos((2*pi*i)/60) );
end
freq_hann=freqz(w_han); %Finding the frequency response of Hanning window

%%
% * |*Matlab Commands for Bartlet window  *|
for i=1:60 %Looping through the 61 elements of the frequency response
    if i>=0 && i<=30 
        w3(i+1)=(2*i)/60;
    else
        w3(i+1)=2- ( (2*i)/60) ;
    end
end
freq_bar=freqz(w3); %Finding the frequency response of Bartlet window

%%
% * |*Matlab Commands for Blackman window  *|
for i=0:60 %Looping through the 61 elements of the frequency response
    w4(i+1)=0.42 - ( 0.5*cos((2*pi*i)/60) ) +  ( 0.08*cos((4*pi*i)/60) ) ;
end
freq_black=freqz(w4); %Finding the frequency response of Blackman window

%% Results:
% * |*Plot for the Question No 1(a)*|
figure;plot(h1);
title('Impulse response rectangular windowed ideal low pass filter ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 1(b)*|
figure;plot(freq);
title(' 512 point frequency response rectangular ideal low pass filter ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 1(c)*|
figure;freqz(h1);
title(' Frequency Response for Rectangular ideal low pass filter ');
xlabel('index');ylabel('Amplitude in dB');

%%
% * |*Plot for the Question No 2(i)*|
figure;plot(w1);
title(' Impulse response for Rectangular window ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 2(i)(a)*|
figure;freqz(w1);
title(' Frequency response for Rectangular window ');
xlabel('index');ylabel('Amplitude in dB');

%%
% * |*Plot for the Question No 2(iii)*|
figure;plot(w_ham);
title(' Impulse response for Hamming Window ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 2(iii)(a)*|
figure;freqz(w_ham);
title(' Frequency response for Hamming Window ');
xlabel('index');ylabel('Amplitude in dB');

%%
% * |*Plot for the Question No 2(iv)*|
figure;plot(w_han);
title(' Impulse response for Hanning Window ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 2(iv)(a)*|
figure;freqz(w_han);
title(' Frequency response for Hanning Window ');
xlabel('index');ylabel('Amplitude in dB');

%%
% * |*Plot for the Question No 2(v)*|
figure;plot(w3);
title(' Impulse response for Bartlet Window ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 2(v)(a)*|
figure;freqz(w3);
title(' Frequency response for Bartlet Window ');
xlabel('index');ylabel('Amplitude in dB');

%%
% * |*Plot for the Question No 2(vi)*|
figure;plot(w4);
title(' Impulse response for Blackman Window ');
xlabel('index');ylabel('Amplitude');

%%
% * |*Plot for the Question No 2(vi)(a)*|
figure;freqz(w4);
title(' Frequency response for Blackman Window ');
xlabel('index');ylabel('Amplitude in dB');

%%
