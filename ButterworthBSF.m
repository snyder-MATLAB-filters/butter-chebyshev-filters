clear all;format compact;format long;clc;
Ap=1;As=24;f1=1100;f2=1600;f3=2100;f4=2700;fsampling=8000;excess='pass';
%f to theta
theta1=2*pi*f1/fsampling
theta2=2*pi*f2/fsampling
theta3=2*pi*f3/fsampling
theta4=2*pi*f4/fsampling
%Prewarp
w1=2*fsampling*tan(theta1/2)
w2=2*fsampling*tan(theta2/2)
w3=2*fsampling*tan(theta3/2)
w4=2*fsampling*tan(theta4/2)
%Backward transformation
Ws1=w2*(w4-w1)/(w1*w4-w2^2)
Ws2=w3*(w4-w1)/(w3^2-w1*w4)
Wp=1
Ws=min(abs(Ws1),abs(Ws2))
%Cutoff frequency adjustments to meet the excess tolerance requirements.
w2a=w2;w3a=w3;
if abs(Ws1)>abs(Ws2)
    w2a=w1*w4/w3;
else
    w3a=w1*w4/w2;
end
w2a
w3a
theta2a=2*atan(w2a/(2*fsampling))
theta3a=2*atan(w3a/(2*fsampling))

%Step 4 Order n of normalized LPF
n=log10((10^(0.1*As)-1)/(10^(0.1*Ap)-1))/(2*log10(Ws))
n=ceil(n)

%Recalculation of Ap, As, epsilon, 3dB cutoff frequency Wc
Ap1=Ap;
As1=As;
if excess == lower('stop')
    As1=10*log10(Ws^(2*n)*(10^(0.1*Ap1)-1)+1);
    epsilon=sqrt(10^(0.1*As1)-1)/Ws^n;
    Wc=1/epsilon^(1/n);
else
    Ap1=10*log10((10^(0.1*As1)-1)/Ws^(2*n)+1);
    epsilon=sqrt(10^(0.1*Ap1)-1);
    Wc=1/epsilon^(1/n);
end
Ap1
As1
epsilon

%Step 5 Normalized lowpass poles and zeros
k=0:ceil((n-3)/2);
S=exp(j*pi*(2*k+n+1)/(2*n));
S=[S,conj(S)];
S=cplxpair(S);

%check if odd
if rem(n,2)==1 S=[-1,S]; end
%poles
S=reshape(S,[1,n])
%zeros at inf
SZ=realmax*ones(1,n)
S=Wc*S

%
%Frequency transformation
wo=sqrt(w1*w4)
s1=(1./(2*S)).*(w4-w1+sqrt((w4-w1)^2-4*w1*w4*S.^2));
s2=(1./(2*S)).*(w4-w1-sqrt((w4-w1)^2-4*w1*w4*S.^2));
s=[s1,s2];
s=cplxpair(s)
sza1=j*sqrt(w1*w4)
sza2=-j*sqrt(w1*w4)
sza=[sza2,sza1,sza2,sza1,sza2,sza1];
sa2=reshape(s,[2,n]);
sz=sza;
% if rem(n,2)==1 sz=[-j*wo,j*wo,sza]; end
sz
%Digital poles and zeros
p=(2*fsampling+s)./(2*fsampling-s)
p1=reshape(p,[2,n]);
z=(2*fsampling+sz)./(2*fsampling-sz)
z1=reshape(z,[2,n]);
DenSOS=[ones(1,n);-sum(p1);prod(p1)]
NumSOS=[ones(1,n);-sum(z1);prod(z1)]

%b0
wp=sqrt(w1*w4)
thetap=2*atan(wp/(2*fsampling))
b0=sum(DenSOS)./sum(NumSOS)
b0b=prod(b0)
b01=ones(3,1)*b0
Den=poly(p)
Num=poly(z)
singleb0=b0b
Num1=singleb0*Num
NumSOS1=b01.*NumSOS
%Plots
x=-1:0.05:1;
ty=sqrt(1-x.^2);
figure(1);
plot(x,ty,':b',x,-ty,':b');
hold on;
plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
axis([-1,1,-1,1]);axis square;
nstr=num2str(n);
title('Digital Pole-Zero Diagram');hold off;
theta=0:0.002:pi;
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
xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
title('Magnitude response (dB)');
figure(4);plot(theta/pi,angle(H)*180/pi,'LineWidth',2);grid on;
xlabel('\theta/\pi');ylabel('\angle(H(\theta))');
title('Phase response');
H0=abs(H(1))
H1=abs(H(length(H)))
wo=sqrt(w1*w4)
thetao=2*atan(wo/(2*fsampling));
thetac=[theta1,theta2a,theta2,theta3,theta3a,theta4,thetao];
Qc=exp(-j*[0:2]'*thetac);
Hc=singleb0*Num*exp(-j*[0:length(Num)-1]'*thetac) ...
./(Den*exp(-j*[0:length(Den)-1]'*thetac));
Hc=20*log10(abs(Hc))
Ap1
As1