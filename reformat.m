%Kod pobieraj¹cy dane binarne z USRP i przetwarzaj¹ce je na serie próbek
%   zespolonych do dalszej analizy
%
%Arkadiusz Majdecki
clear all
close all
sizeA=211040000;
fileID = fopen('final1.bin');
A = fread(fileID,sizeA,'uint16');
%Zakladam tu ¿e USRP formatuj¹cy dane w complex int16, przez pierwsze 16 bitów wystawia
%   wartoœæ rzeczywist¹ próbki, a przez nastêpne 16 urojon¹, po czym przechodzi
%   do nastêpnej próbki
Areal=A(1:2:sizeA-1);
Aimag=A(2:2:sizeA);
B=Areal(:)+1i.*Aimag(:);
save('B.mat','B');