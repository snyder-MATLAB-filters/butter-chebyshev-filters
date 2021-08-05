clear all;format long;clc;
Ap=2;As=31;
fp=1500;fs=2100;
fsampling=8000;

excess='stop';
% excess='pass';

%Step 1 Digital Frequencies (f ---> theta)
thetap=2*pi*fp/fsampling
thetas=2*pi*fs/fsampling

%Step 2 Prewarp (theta ---> omega)
wp=2*fsampling*tan(thetap/2)
ws=2*fsampling*tan(thetas/2)

%Step 3 Backward transformation
Wp=1
Ws=ws/wp

%Order of the Normalized Lowpass Filter (Step 4)
n=acosh(sqrt((10^(0.1*As)-1)/(10^(0.1*Ap)-1)))/acosh(Ws)
n=ceil(n)
e=sqrt(10^(0.1*Ap)-1);

%Recalculation of Ap, As, e
Ap1=Ap;
As1=As;
e1=e;
if excess == lower('pass')
    Ap1=10*log10(1+(10^(0.1*As)-1)/cosh(n*acosh(Ws))^2);
    e1=sqrt(10^(0.1*Ap1)-1);
else
    As1=10*log10(1+(10^(0.1*Ap)-1)*cosh(n*acosh(Ws))^2);
end   
Ap1
As1
e1

%Poles and Zeros of the Normalized Lowpass Filter (Step 5)
a=(1/n)*asinh(1/sqrt(10^(0.1*Ap1)-1))
k=0:ceil((n-3)/2);
Sa=-sinh(a)*sin(pi*(2*k+1)/(2*n))+j*cosh(a)*cos(pi*(2*k+1)/(2*n));
Sa=[Sa,conj(Sa)];
Sa=cplxpair(Sa);
Sb=reshape(Sa,[2,ceil((n-1)/2)]);
S=Sa;
if rem(n,2)==1 S=[-sinh(a),Sa]; end
S
% S=[-0.59,-0.3+1j,-0.3-1j]
% w2=10470
% w3=15870
% S=[-0.56,-0.28+1j,-0.28-1j]
% w2=11283
% w3=16919
% S=[-0.54,-0.27+1j,-0.27-1j]
% w2=12127
% w3=18027
% S=[-0.51,-0.26+1j,-0.26-1j]
% w2=12914
% w3=17913
% S=[-0.48,-0.24+1j,-0.24-1j]
% w2=13823
% w3=19079
SZ=realmax*ones(1,n)

%Poles and Zeros of the Frequency Transformed Lowpass Filter (Step 6)
s=wp*S
sz=SZ

%Step 7 Bilinear transformation
pa=(2*fsampling+s)./(2*fsampling-s)
za=(2*fsampling+sz)./(2*fsampling-sz)
pa=cplxpair(pa)
za=cplxpair(za)
if rem(n,2)==1 p=[pa,0];z=[za,0]; else p=pa;z=za; end
p
z

%Step 8 Second order sections and single section
p1=reshape(p,[2,floor((n+1)/2)]);
z1=reshape(z,[2,floor((n+1)/2)]);
DenSOS=[ones(1,floor((n+1)/2));-sum(p1);prod(p1)]
NumSOS=[ones(1,floor((n+1)/2));-sum(z1);prod(z1)]
Den=poly(p);
Num=poly(z);

%Finding the Constant Term b0 (Step 9)
%Step 9 Constant terms
%b0
b0=sum(DenSOS)./sum(NumSOS)
singleb0=prod(b0)

%Plots
x=-1:0.05:1;
ty=sqrt(1-x.^2);
figure(1);
plot(x,ty,':b',x,-ty,':b');
hold on;
plot(real(pa),imag(pa),'bx','MarkerSize',12);hold on;
plot(real(za),imag(za),'bo','MarkerSize',12);grid on;
axis([-1,1,-1,1]);axis square;
nstr=num2str(n);
text(-0.9,0.1,nstr,'FontSize',12);
title('Digital Pole-Zero Diagram');hold off;
theta=0:0.002:pi-0.03;
Q=exp(-j*[0:2]'*theta);
H=singleb0*Num*exp(-j*[0:length(Num)-1]'*theta) ...
./(Den*exp(-j*[0:length(Den)-1]'*theta));
figure(2);plot(theta/pi,abs(H),'LineWidth',2);grid on;
xlabel('\theta/\pi');ylabel('|H(\theta)|');
title('Magnitude response (linear)');
Mag=20*log10(abs(H));
Th=-100;
Mag=(Mag >= Th).*Mag + Th*(Mag < Th);
figure(3);plot(theta/pi,Mag,'LineWidth',2);grid on;
xlabel('\theta/\pi');ylabel('20*log_1_0(|H(\theta)|)');
title('Magnitude response (dB)');
figure(4);plot(theta/pi,angle(H)*180/pi,'LineWidth',2);grid on;
xlabel('\theta/\pi');ylabel('\angleH(\theta)');
title('Phase response');

% %Step 10 (verify)
% H0=abs(H(1))
% H1=abs(H(length(H)))
% thetao=2*atan(wo/(2*fsampling));
% thetac=[thetap,thetas,thetao];
% Qc=exp(-j*[0:2]'*thetac);
% Hc=singleb0*Num*exp(-j*[0:length(Num)-1]'*thetac) ...
% ./(Den*exp(-j*[0:length(Den)-1]'*thetac));
% Hc=20*log10(abs(Hc))
% Ap1
% As1
