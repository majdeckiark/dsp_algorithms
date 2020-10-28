%Algorytm symuluj¹cy sygna³ zdudnieñ przy zadanych parametrach i rozk³adzie
%   obiektów, po czym przeprowadza na nim analizê furierowsk¹ w celu
%   wyznaczenia mapy range-doppler
%
%Arkadiusz Majdecki

clear all
close all

%DANE PARAMETRYCZNE
F_if_max=10e6;%Maksymalna czestotlwiosc sygnalu zdudnien
N_f=1e3;%Ilosc komorek w fast time
N_s=1e2;%ilosc komorek w slow time 
T_c=1e-4;%Czas trwania chirpu
T_i=T_c*N_s;%Czas integracji
B=160e6;%Pasmo 
f_0=9.6e9;%czestotliwosc nosna
c=3e8;%m/s
lambda_0=c/f_0;%dlugosc fali nosnej
S=B/T_c;%Nachylenie krzywej modulacji

%WARTOSCI GRANICZNE
Vmax=lambda_0/4/T_c;%predkosc maksymalna
Rmax=F_if_max*c/2/S;%Odleglosc maksymalna
Vres=lambda_0/4/T_i;%rozdzielczosc predkosci
Rres=c/2/B;         %rozdzielczosc odlegosci


%WEKTOR CZASU ZDEFINIOWANY DLA POJEDYNCZEGO CHIRPU
time=linspace(0,T_c,N_f);

%ZADANE ODLEGLOSCI I PREDKOSCI OBIEKTOW
objrng=[300,300,150,150];
objspd=[-30,0,50,-45];

%GENERACJA SYGNALU ZDUDNIEN DLA DANEGO SCENARIUSZA
objshift=objspd.*T_c;%zmiana odleglosci po jednym okresie chirpu
ifs=zeros(N_s,N_f);%alokacji macierzy sygna³u
for k=1:N_s%slow time 
    for l=1:numel(objrng)%petla po obiektach
        fif=S*2*objrng(l)/c;%f_if=S2d/c dla aktualnie opracowywanego obiektu
        ifs(k,:)=ifs(k,:)+exp(1i*(2*pi*fif.*time+4*pi*k*objshift(l)/lambda_0));%dodanie do sygnalu skladowej wynikajacej z odbicia od aktualnie opracowywanego obiektu
    end                                                         
end  
%analiza czasowo-czêstotliwoœciowa
RTM=(fft(ifs,[],2));
RDM=fftshift(fft(RTM,[],1),1); 

%UTWORZENIE WEKTORÓW DZIEDZIN ODLEG£OŒCI, PRÊDKOŒCI
Rdomain=linspace(0,Rmax,N_f);%Wynik range-fft (czestosc ko³owa zmian czestotliwosc) jest liniowo zale¿na od odleg³osci R=F_if*c/2/S
Vdomain=linspace(-Vmax,Vmax,N_s);%Wynik time-fft (czestosc kolowa zmian fazy) jest liniowo zalezna od predkosci V=lambda_o*omega/4/pi/T_c

%RYSOWANIE WYKRESÓW
mesh(linspace(0,Rmax,N_f),linspace(-Vmax,Vmax,N_s),abs(RDM));
save('mtx.mat','Rdomain','Vdomain')