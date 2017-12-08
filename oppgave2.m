clear all

lib = oblig2_lib;

%1. Lese inn av bildet
im = double(imread('uio.png'));
entro= lib.my_entropy(uint8(im));
[M,N]=size(im);

figure
imshow(uint8(im)); title('Original Image');

%Kvantisering Matrise 
Q=[16 11 10 16 24 40 51 61;
12 12 14 19 26 58 60 55;
14 13 16 24 40 57 69 56;
14 17 22 29 51 87 80 62;
18 22 37 56 68 109 103 77;
24 35 55 64 81 104 113 92;
49 64 78 87 103 121 120 101;
72 92 95 98 112 100 103 99];

%Kavantisering matrise med Faktor
q=1;
Q=Q*q;

% 2.Trekk 128 fra bildet
im= im-128;

% 3. 8x8 Blokkvis DCT transform av bildet
F = lib.DCT_2D(im,8,8);

% 5. Kvantisering av DCT koeffisientene
Fq = lib.quantizier(F,Q);

% 6. Beregning av det forventet entropien av DCT bildet
%    og kompresjons raten
entro_comp=lib.my_entropy(Fq);
cr=entro/entro_comp;

% 7.a Restoring DCT matrisen
Fdq = lib.dequantizier(Fq,Q);

% 7.b IDCT transform of legge 128 til resultaten
f = lib.IDCT_2D(Fdq,8,8);
f=floor(f+128);

% Klipping ut de null utvided verdiene og plotting
f=f(1:M,1:N);
str=strcat('Compressed image cr =  ', num2str(cr));
figure
imshow(uint8(f));title(str);

