clear
clc 
close all

IM = imread("meme.jpg");
A = cast(IM,"double")/255;
B = rot90(A);
C = rot90(B);
D = rot90(C);

dA = abs(diff(A));
dB = abs(diff(B));
dC = abs(diff(C));
dD = abs(diff(D));

dA1 = dA(:,1:end-1,:);
dB1 = rot90(dB(:,1:end-1,:),3);
dC1 = rot90(dC(:,1:end-1,:),2);
dD1 = rot90(dD(:,1:end-1,:),1);


size(dA1)
size(dB1)
size(dC1)

E = (dA1+dB1+dC1+dD1)/4;

%A = cast(A,"double")/255;
E = cast(E,"double")/255*1000;

figure
imshow(A)
figure
imshow(E)
figure
%%
clear
clc 

IM = imread("meme.jpg");
A = (cast(IM,"double")/255);

threshg = 0.61;
thresh1 = 0.8;
thresh2 = 0.6;
thresh3 = 0.3;

G = rgb2gray(A);
A1 = A(:,:,1);
A2 = A(:,:,2);
A3 = A(:,:,3);

wA1 = A1>thresh1;
wA2 = A2>thresh2;
wA3 = A3>thresh3;
wG = G>threshg;

wwG = medfilt2(wG,[5 5]);
wwA1 = medfilt2(wA1,[5 5]);
wwA2 = medfilt2(wA2,[5 5]);
wwA3 = medfilt2(wA3,[5 5]);

figure(1)
imshow(wwG);
%%
figure(1)
imshow(A1>thresh1);
figure(2)
imshow(A2>thresh2);
figure(3)
imshow(A3>thresh3);
figure(4)
imshow(G>threshg);
%%

% figure(2)
% imshow(wwA1);
% figure(3)
% imshow(wwA2);
% figure(4)
% imshow(wwA3);