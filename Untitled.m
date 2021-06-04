clc

clear all

j=sqrt(-1);

Nt=16;

Nr=Nt/2;

%Generate the NrxNt MIMO channel matrix

H0=1/sqrt(2)*(randn(Nr,Nt)+j*randn(Nr,Nt));

dis=100;

v=3*10^8;

fc=1.2*10^9;

wavelength=v/fc;

%FSPL=(4*pi*dis/wavelength)^2

FSPL=(4*pi*dis/wavelength)^2

H1=FSPL^(-1/2).*H0;

noise=1;

MTPNR=40;

x=10^(MTPNR/10);

Pmax=noise*(x)

%Pmax=10^(40/10);

[V,D] = eig(H1'*H1);

x = eig(H1'*H1);
K1=sort(real(x),'descend'); 
figure, plot (K1,'o--');


title('PLOT OF EIGEN VALUES');

xlabel('INDEX');

ylabel('EIGEN VALUES');