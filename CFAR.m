%Arkadiusz Majdecki

function det=CFAR(RDM)

%DANE PARAMETRYCZNE
P_fa=1e-130;
N_g=2;%ilosc komorek ochronnych (guard-cells)
N_t=5;%ilosc komorek branych do sumy (training cells)

%ZA£ADOWANIE DANYCH WEJŒCIOWYCH
load('mtx.mat');%Pobranie Rdomain i Vdomain do wyskalowania macierzy det

%GENERACJA MASEK (testowej i ochronnej)
%   +1 za komórke centraln¹ (cell under test)
Sum_mask=ones(N_t+1,N_t+1);
Dif_mask=ones(N_g+1,N_g+1);
Train_mask_size=(2*N_t+2*N_g+1)^2-(2*N_g+1)^2; %iloœæ komórek w masce "treningowej"

%POWIELENIE MATRYCY
%   Powielenie matrycy ma na celu oddanie jej 
%   okresowoœci podczas wyliczania sumy funkcj¹ filter2 
%   jednak¿e zwiêksza nak³ad obliczeniowy
mtx=[mtx mtx mtx; mtx mtx mtx; mtx mtx mtx];

%WYLICZENIE ŒREDNIEJ WARTOŒCI KOMÓREK TRENINGOWYCH
Sum=filter2(Sum_mask,mtx);
Dif=filter2(Dif_mask,mtx);
Avg_mtx=(Sum-Dif)/(Train_mask_size);

%PRZYCIÊCIE MACIERZY W CELU POWROTU DO ROZMIARU MACIERZY WEJŒCIOWEJ
mtx_sz=size(RDM);
Avg_mtx=Avg_mtx(mtx_sz(1)+1:2*mtx_sz(1),mtx_sz(2)+1:2*mtx_sz(2));

%WYLICZENIE WARTOSCI THARSHOLDU
%Nie mam pewnoœci co do linii 37, ze wzglêdu na kosmiczne wartoœci Pfa
%   wymaganej do uzysakania w miarê uzytecznego wyniku detekcji
alfa=Train_mask_size*(P_fa^(-1/Train_mask_size)-1);
Avg_mtx=alfa*Avg_mtx;

%PRZYROWNANIE MACIERZY RANGE-DOPPLER DO MACIERZY THRASHOLDU
det=abs(RDM)>abs(Avg_mtx);
mesh(Rdomain,Vdomain,det)