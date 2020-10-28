%Algorytm symuluje sygna³ zdudnieñ przy zadanych parametrach i rozk³adzie
%   obiektów, po czym przeprowadza na nim analizê furierowsk¹ w celu
%   wyznaczenia mapy range-doppler-azymut
%
%Arkadiusz Majdecki

clear all;
close all;

%DANE PARAMETRYCZNE
F_if_max=10e6;%Maksymalna czestotlwiosc sygnalu zdudnien
N_f=1e3;%Ilosc komorek w fast time
N_s=1e2;%ilosc komorek w slow tim
N_a=20;%Ilosc anten
T_c=1e-4;%Czas trwania chirpu
T_i=T_c*N_s;%Czas integracji
B=160e6;%Pasmo 
f_0=9.6e9;%czestotliwosc nosna
c=3e8;%m/s
lambda_0=c/f_0;%dlugosc fali nosnej
S=B/T_c;%Nachylenie krzywej modulacji
d=lambda_0/2;%rozstaw anten ZAK£ADAM TU MAKSYMALNY K¥T AMAX=90stp

%WARTOSCI GRANICZNE
Vmax=lambda_0/4/T_c;%predkosc maksymalna
Rmax=F_if_max*c/2/S;%Odleglosc maksymalna
Vres=lambda_0/4/T_i;%rozdzielczosc predkosci
Rres=c/2/B;         %rozdzielczosc odlegosci
Amax=asin(lambda_0/2/d);%azymut maksymalny


%WEKTOR CZASU ZDEFINIOWANY DLA POJEDYNCZEGO CHIRPU
time=linspace(0,T_c,N_f);

%ZADANE ODLEGLOSCI I PREDKOSCI OBIEKTOW
objrng=[400,700];
objspd=[30,-50];
objang=[-pi/4,pi/6];

%GENERACJA SYGNALU ZDUDNIEN DLA DANEGO SCENARIUSZA
objshift=objspd.*T_c;%zmiana odleglosci po czasie integracji 
ifs=zeros(N_s,N_f,N_a);%alokacji macierzy sygna³u

for ant=1:N_a%petla po antenach
    for k=1:N_s%slow time 
        for l=1:numel(objrng)%petla po obiektach
            fif=S*2*objrng(l)/c;%f_if=S2d/c dla aktualnie opracowywanego obiektu
            ifs(k,:,ant)=ifs(k,:,ant)+exp(1i*(2*pi*fif.*time+4*pi*k*objshift(l)/lambda_0+2*pi*d*sin(objang(l))*ant/lambda_0));%dodanie do sygnalu skladowej wynikajacej z odbicia od aktualnie opracowywanego obiektu
        end                                                        
    end
end    

%ANALIZA CZASOWO CZÊSTOTLIWOŒCIOWA
RAM=fftshift(fft(ifs,[],3),3);%fft shift bo dziedzina parzysta
RTM=(fft(RAM,[],2));
RDM=fftshift(fft(RTM,[],1),1); 

%UTWORZENIE WEKTORÓW DZIEDZIN ODLEG£OŒCI, PRÊDKOŒCI, AZYMUTU
Rdomain=linspace(0,Rmax,N_f);%Wynik range-fft (czestosc ko³owa zmian czestotliwosc) jest liniowo zale¿na od odleg³osci R=F_if*c/2/S
Vdomain=linspace(-Vmax,Vmax,N_s);%Wynik time-fft (czestosc kolowa zmian fazy) jest liniowo zalezna od predkosci V=lambda_o*omega/4/pi/T_c
Adomain=90./Amax.*asin(linspace(ceil(-N_a/2),ceil(N_a/2),N_a)/ceil(N_a/2));%Wynik anthen-fft (czestosc kolowa zmian fazy) jest sinusoidalnie zalezny od kata Theta=asin(lambda_0*omega/2/pi/d)


%RYSOWANIE WYKRESÓW
%title('Macierz Range-Doppler-Azymut po wyskalowaniu');
%   Wstawienie tytu³u korumpuje w jakiœ sposób isosurface
isosurface(Rdomain,Vdomain,Adomain,abs(RDM(:,:,:)));
view(0,90);
figure()
%title('Macierz Range-Doppler-Azymut');
isosurface(abs(RDM(:,:,:)));
view(0,90);

%PRZYK£ADOWY WYCINEK PRZEZ PRÊDKOŒÆ -50 [M/S] ~ 19 [KOMÓRKA W PRÊDKOŒCI]
RD(:,:)=RDM(19,:,:);
figure()
mesh(Adomain, Rdomain, abs(RD))