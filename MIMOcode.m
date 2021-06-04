clear all
clc
Nt = 16;
Nr = 8;
fc=1.2*10^9; %Carrier Frequency
d=100;% distance
c=3*10^8;
BW=20*10^6;
l=c/fc; 

FSPL=(4*pi*d/l)^2;
Pmax = 10000; %40db=10log(Pmax/N) then we calculate the maximum transmitted power
j=sqrt(-1);



if 1
    H0=(1/sqrt(2))*(randn(Nr,Nt)+j*randn(Nr,Nt));

save Channel H0
else 
    load Channel H0
end;

H=H0*(FSPL)^(-1/2);  

%Plot eigen values
A=H'*H
E=eig(A)

K1=sort(real(E),'descend'); 
figure, plot(K1,'o-')
[U,D] = eig(A)

%nth diagonal elemetnt of matrix D and nth coulmn of matrix U
%reconstructing matrix A

B=zeros(Nt,Nt);
for n=1:Nt;
    B=B+D(n,n)*U(:,n)*U(:,n)';
end;
K=eig(B);
K2=sort(real(K),'descend');
hold on,plot(K2,'*--');

%Singular value decomposision for matrix H

[Us,Ds,Vs]=svd(H); %Singular Value Decomposition
product=real(Us'*H*Vs)%real part of complex number
SingularVal=sort(diag(real(Ds)),'ascend')

%Question 5
%data rate

Nall=find(Ds>1e-6);

P1=Pmax/1;
C1=log2(1+P1*Ds(1,1)^2);

P2=Pmax/2;
C2=log2(1+P2*Ds(1,1)^2)+log2(1+P2*Ds(2,2)^2);

P3=Pmax/3;
C3=log2(1+P3*Ds(1,1)^2)+log2(1+P3*Ds(2,2)^2)+log2(1+P3*Ds(3,3)^2);

P4=Pmax/4;
C4=log2(1+P4*Ds(1,1)^2)+log2(1+P4*Ds(2,2)^2)+log2(1+P4*Ds(3,3)^2)+log2(1+P4*Ds(4,4)^2);

P5=Pmax/5;
C5=log2(1+P5*Ds(1,1)^2)+log2(1+P5*Ds(2,2)^2)+log2(1+P5*Ds(3,3)^2)+log2(1+P5*Ds(4,4)^2)+log2(1+P5*Ds(5,5)^2);

P6=Pmax/6;
C6=log2(1+P6*Ds(1,1)^2)+log2(1+P6*Ds(2,2)^2)+log2(1+P6*Ds(3,3)^2)+log2(1+P6*Ds(4,4)^2)+log2(1+P6*Ds(5,5)^2)+log2(1+P6*Ds(6,6)^2);

P7=Pmax/7;
C7=log2(1+P7*Ds(1,1)^2)+log2(1+P7*Ds(2,2)^2)+log2(1+P7*Ds(3,3)^2)+log2(1+P7*Ds(4,4)^2)+log2(1+P7*Ds(5,5)^2)+log2(1+P7*Ds(6,6)^2)+log2(1+P7*Ds(7,7)^2);

P8=Pmax/8;
C8=log2(1+P8*Ds(1,1)^2)+log2(1+P8*Ds(2,2)^2)+log2(1+P8*Ds(3,3)^2)+log2(1+P8*Ds(4,4)^2)+log2(1+P8*Ds(5,5)^2)+log2(1+P8*Ds(6,6)^2)+log2(1+P8*Ds(7,7)^2)+log2(1+P8*Ds(8,8)^2);


C=[C1 C2 C3 C4 C5 C6 C7 C8 ];
DataRate=BW*C
x=1:8;
figure, plot(x,DataRate,'--');

%Question 6

SingularValues=find(Ds>1e-6);
mu=1000;
epsilon=1e-9;
 
Pi=subplus(mu-(1./(SingularValues.^2)));
if sum(Pi)>Pmax
    mu=mu-mu/2;
    Pi=subplus(mu-(1./(SingularValues.^2)));
end

if sum(Pi)<Pmax-ee
        mu=mu+mu/2;
        Pi=subplus(mu-(1./(SingularValues.^2)));
end
if Pmax>sum(Pi)&& sum(Pi)>Pmax-ee
    save mu
end
c=(log2(1+Pi.*(SingularValues.^2)));
r=BW*c;






