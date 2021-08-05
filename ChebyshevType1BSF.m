clear all;format long;clc;
Ap=2;As=45;
f1=1100;f2=1600;f3=2100;f4=2700;fsampling=8000;

excess='stop';
% excess='pass';

%Digital Cutoff Frequencies (Step 1)
theta1=2*pi*f1/fsampling
theta2=2*pi*f2/fsampling
theta3=2*pi*f3/fsampling
theta4=2*pi*f4/fsampling

%Prewarping (Step 2)
w1=2*fsampling*tan(theta1/2)
w2=2*fsampling*tan(theta2/2)
w3=2*fsampling*tan(theta3/2)
w4=2*fsampling*tan(theta4/2)

%Backward Transformation (Step 3)
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

%Order of the Normalized Lowpass Filter (Step 4)
n=acosh(sqrt((10^(0.1*As)-1)/(10^(0.1*Ap)-1)))/acosh(Ws)
n=ceil(n)
e=sqrt(10^(0.1*Ap)-1)

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

SZ=inf*ones(1,n)

%Poles and Zeros of the Frequency Transformed Bandpass Filter (Step 6)
s1=(1./(2*S)).*(w4-w1+sqrt((w4-w1)^2-4*w1*w4*S.^2));
s2=(1./(2*S)).*(w4-w1-sqrt((w4-w1)^2-4*w1*w4*S.^2));
s=[s1,s2];
s=cplxpair(s)

sza1=j*sqrt(w1*w4);
sza2=-j*sqrt(w1*w4);
sza=[sza2,sza1,sza2,sza1,sza2,sza1,sza2,sza1];
sa2=reshape(s,[2,n]);
sz=sza
%Poles and Zeros of the Digital Filter (Step 7)
p=(2*fsampling+s)./(2*fsampling-s)
z=(2*fsampling+sz)./(2*fsampling-sz)
p1=reshape(p,[2,n]);
z1=reshape(z,[2,n]);

%Second Order Sections (Step 8)
DenSOS=[ones(1,n);-sum(p1);prod(p1)]
NumSOS=[ones(1,n);-sum(z1);prod(z1)]
%Finding the Constant Term b0 (Step 9)
wp=sqrt(w2*w3)
thetap=2*atan(wp/(2*fsampling))
Q=exp(-1j*[0:2]'*thetap);
b0a=sum(DenSOS)./sum(NumSOS)
b0b=prod(b0a);
if rem(n,2)==0 b0=b0b*10^(-0.05*Ap1); else b0=b0b; end
%if rem(n,2)==0 & excess=='pass' b0=b0b/sqrt(1+epsilon^2); end
b0
Den=poly(p)
Num=poly(z)
singleb0=b0
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
text(0.85,0.1,nstr,'FontSize',12);
text(-0.9,0.1,nstr,'FontSize',12);
title('Digital Pole-Zero Diagram');hold off;
theta=0:0.002:pi;
Q=exp(-1j*[0:2]'*theta);
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
xlabel('\theta/\pi');ylabel('\angle(H(\theta))/\pi');
title('Phase response');
%Verify (Step 10)
H0=abs(H(1))
H1=abs(H(length(H)))
wo=sqrt(w2*w3)
thetao=2*atan(wo/(2*fsampling));
thetac=[theta1,theta2,theta2a,theta3,theta3a,theta4,thetao];
Qc=exp(-j*[0:2]'*thetac);
Hc=singleb0*Num*exp(-j*[0:length(Num)-1]'*thetac) ...
./(Den*exp(-j*[0:length(Den)-1]'*thetac));
Hc=20*log10(abs(Hc))
Ap1
As1
