%Kod pobieraj�cy dane binarne z USRP i przetwarzaj�ce je na serie pr�bek
%   zespolonych do dalszej analizy
%
%Arkadiusz Majdecki
clear all
close all
sizeA=211040000;
fileID = fopen('final1.bin');
A = fread(fileID,sizeA,'uint16');
%Zakladam tu �e USRP formatuj�cy dane w complex int16, przez pierwsze 16 bit�w wystawia
%   warto�� rzeczywist� pr�bki, a przez nast�pne 16 urojon�, po czym przechodzi
%   do nast�pnej pr�bki
Areal=A(1:2:sizeA-1);
Aimag=A(2:2:sizeA);
B=Areal(:)+1i.*Aimag(:);
save('B.mat','B');