%% Digital Signal Processing |[Lab-1]| 
% * Authors: _Kshitij Srivastava(1510110200)_ and _Nilambar Saha(1510110246)_
% * Lab Instructor: _Dr. Ravi Kant Saini_ 
%% Objective: 
%Data Handling -In this experiment we tried to convert all types of files 
%into its matrix and then it was converted to a different file type using 
%transformation of the matrix.
%% Program: 
clc;
clear all;
close all;

% * |*Matlab Commands for Speech or Audio Reading and Writing*|
[y, fs]=audioread('Signal_Processing_Audio.mp3');
y_n=y(:,1);
t=0:1/fs:2;
sz = size(y_n);
sx = size(t);
isvector(y_n);
ismatrix(y_n);
y_t=y_n(1:88201);

[m, n]=size(y_t);
Y=reshape(y_t,[m,n]);
audiowrite('output_sound.wav',Y,fs);
sound(Y, fs);
%%
%|*Matlab Commands for Image Reading and Writing*|
IMG = imread('RGB_Image.jpg');
redImage = IMG(:,:,1);
greenImage=IMG(:,:,2);
blueImage=IMG(:,:,3);

I = rgb2gray(IMG);

MAT=reshape(IMG,1,[]);
isvector(MAT);
REIMG=imresize(IMG, 0.5);

%%
%|*Matlab Commands for Video File Reading and Writing*|
V = VideoReader('Signal_Processing_Video.mp4');
numFrames = V.NumberOfFrames;

NEWV = VideoWriter('new_video.avi');
open(NEWV)
for k = 1:60
   img = read(V,k);
   writeVideo(NEWV,img);
end
close(NEWV);

VMAT=zeros(518400,1);
for k = 1:60 %less iterations due to lots of Data
   img = read(V,k);
   B = reshape(img,[518400,1]);
   VMAT = horzcat(VMAT,B);
end

%%
%|*Matlab Commands for Excel File Reading and Writing*|
[num, str]=xlsread('Text_Data.xlsx');
STR_ARR=char(str);
ASCII_STR=double(STR_ARR);
[a, b]=size(ASCII_STR);
MAT_STR = reshape(ASCII_STR,[143,3]);

%%
%|*Matlab Commands for Audio files to ASCII value*|
[AUD, fs]=audioread('output_sound.wav');
AUD_ARR=num2str(AUD);
AUD_CHAR=char(AUD_ARR);
AUD_DOB=double(AUD_CHAR);

%%
%|*Matlab Commands for using FLOOR, CEIL, ROUND fuctions*|
[SAM, fs]=audioread('output_sound.wav');
SAM_FLOOR=floor(SAM);
SAM_CEIL=ceil(SAM);
SAM_ROUND=round(SAM);
FLOOR_ERROR=SAM-SAM_FLOOR;
CEIL_ERROR=SAM-SAM_CEIL;
ROUND_ERROR=SAM-SAM_ROUND;

%% Results:
%|*Plot for the Question No 1(a)(i)*|
figure;plot(t,y_t);
title('Speech Signal');xlabel('Time [Sec]');ylabel('Amplitude');

audiowrite('output_sound.wav',Y,fs);
sound(Y, fs);

%%
%|*Plot for the Question No 1(b)*|
figure;imshow(redImage);title('Red pallets from the image');
figure;imshow(greenImage);title('Green pallets from the image');
figure;imshow(blueImage);title('Blue pallets from the image');

figure;imshow(I);title('Color image to gray scale image');

figure;imshow(REIMG);title('Resizing the image');

%%
%|*Results for the Question No 1(c)*|

%the number of frames received is
numFrames

%%
%|*Plot for the Question No 1(d)*|
figure;image(MAT_STR);title('ASCII values from excel to image');

%%
%|*Plot for the Question No 1(e)*|
figure;image(AUD_DOB);title('Audio samples to image');
