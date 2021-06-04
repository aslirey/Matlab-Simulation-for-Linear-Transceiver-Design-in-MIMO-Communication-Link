
clear all
clc
Nt = 16;
Nr = Nt/2;
fc=1.2e9; %Carrier Frequency
d=100;% distance
v=3e8;
BW=20*10^6;
wavelength=v/fc;
noise=1;
MTPNR=40;
 FSPL=(4*pi*d/wavelength)^2
 Pmax=noise*10^(MTPNR/10) 
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
 title('PLOT OF EIGEN VALUES');
 xlabel('INDEX');
 ylabel('EIGEN VALUES');
 [U,D] = eig(A)
 
%nth diagonal elemetnt of matrix D and nth coulmn of matrix U
 %reconstructing matrix A
 
 B=zeros(Nt,Nt);
 for n=1:Nt;
     B=B+D(n,n)*U(:,n)*U(:,n)';
 end;
 K=eig(B);
 K2=sort(real(K),'descend');
 hold on,plot(K2,'*--');%comparing eigen values of B with A
 
 %Singular value decomposision for matrix H
 
 [Us,Ds,Vs]=svd(H) %Singular Value Decomposition
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
 title('Data Rate');

 
 %Question 6
 
 %SingularValues=find(diag(Ds)>1e-6);
 SingularValues = diag(Ds);
 mu=1000;
 epsilon=1e-5;
  

step = mu/2;

m = 0;
flag = 0;
while flag == 0
    m = m +1;
    
 Pi=subplus(mu-(1./(SingularValues.^2)));
 if sum(Pi)>Pmax
      step = step/2;
     mu=mu-step;
     Pi=subplus(mu-(1./(SingularValues.^2)));
 end
 
 if sum(Pi)<Pmax-epsilon
          step = step;
         mu=mu+step;
         Pi=subplus(mu-(1./(SingularValues.^2)));
 end
% > if Pmax>sum(Pi)&& sum(Pi)>Pmax-ee
% >     save mu
% > end

muall(m) = mu;
Powerall(m) = sum(Pi);

if sum(Pi)<=Pmax && Pmax-sum(Pi)<= epsilon
    flag = 1;
end

end

Pi

 c=(log2(1+Pi.*(SingularValues.^2))) %%% sum of log
 r=BW*c
 figure,plot(r,'o--')
 title('Data Rate_question_six');

figure,plot(muall,'o-')
figure,plot(Powerall,'o-')
hold on,plot(Pmax.*ones(size(Powerall)),'*-')