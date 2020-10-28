%Algorytm tworzenia serii map range-doppler z danych pomiarowych 
%
%Arkadiusz Majdecki
clear all;
close all;
clc;
load('B.mat');%zaladowanie ciagu probek sformatowanych funkcja reformat
mtx=B;

f_samp=10e6;%czestotliwosc probkowania
t_samp=1/f_samp;%czas jednej probki
T_c=1e-4;%czas trwania chirpu
N_f=T_c/t_samp;%ilosc probek w fast timie
N_s=1000;%ilosc okresow branych pod analize vel ilosc probek na slow time
B=160e6;%Pasmo
f_0=9.6e9;%czestotliwosc nosna
c=3e8;%m/s
lambda_0=c/f_0;%dlugosc fali nosnej
S=B/T_c;%Nachylenie krzywej modulacji

Vmax=lambda_0/4/T_c;%predkosc maksymalna
Rmax=f_samp*c/2/S;%Odleglosc maksymalna


mov=VideoWriter('ifexp.avi'); 
open(mov);
for frame=1:100
    frmsz=N_f*N_s;%frame size: ilosc komorek na mape
    frmshft=frame*frmsz;%indeks ostatniej komorki poprzedniej mapy (klatki)
    C=reshape(mtx(1+frmshft:frmsz+frmshft),[1000,1000]);%WYMAGANA RECZNE WPROWADZANIE WIELKOSCI MACIERZY
    RTM=(fft(C',[],2));
    RDM=fftshift(fft(RTM,[],1),1);
    mesh(linspace(0,Rmax,N_f),linspace(-Vmax,Vmax,N_s),abs(RDM));
    F = getframe(gcf); 
    writeVideo(mov,F); 
end
close(mov); 
