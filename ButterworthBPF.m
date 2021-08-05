clear all;format long;clc;
Ap=2;As=29;
f1=1000;f2=1500;f3=2000;f4=2500;fsampling=8000;

% excess='stop';
excess='pass';

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
Wp=1
Ws1=(w2*w3-w1^2)/w1/(w3-w2)
Ws2=(w4^2-w2*w3)/w4/(w3-w2)
Ws=min(abs(Ws1),abs(Ws2))

%Cutoff Frequency Adjustments
w1a=w1;
w4a=w4;
if abs(Ws1) > Ws2
    w1a=w2*w3/w4;
else
    w4a=w2*w3/w1;
end
w1a
w4a
theta1a=2*atan(w1a/(2*fsampling))
theta4a=2*atan(w4a/(2*fsampling))

%Order of the Normalized Lowpass Filter (Step 4)
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

%Poles and Zeros of the Normalized Lowpass Filter (Step 5)
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

%Poles and Zeros of the Frequency Transformed Bandpass Filter (Step 6)
s1=(1/2)*((w3-w2)*S+sqrt(S.^2*(w3-w2)^2-4*w2*w3));
s2=(1/2)*((w3-w2)*S-sqrt(S.^2*(w3-w2)^2-4*w2*w3));
s=[s1,s2];
s=cplxpair(s)

sz=reshape([0,realmax]'*ones(1,n),[1,2*n])
%Poles and Zeros of the Digital Filter (Step 7)
p=(2*fsampling+s)./(2*fsampling-s)
z=(2*fsampling+sz)./(2*fsampling-sz)
p1=reshape(p,[2,n])
z1=reshape(z,[2,n])

%Second Order Sections (Step 8)
DenSOS=[ones(1,n);-sum(p1);prod(p1)]
NumSOS=[ones(1,n);-sum(z1);prod(z1)]
%Finding the Constant Term b0 (Step 9)
wp=sqrt(w2*w3)
thetap=2*atan(wp/(2*fsampling))
Q=exp(-1j*[0:2]'*thetap);
b0a=abs(DenSOS'*Q)./abs(NumSOS'*Q)
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
thetac=[theta1,theta1a,theta2,theta3,theta4a,theta4,thetao];
Qc=exp(-j*[0:2]'*thetac);
Hc=singleb0*Num*exp(-j*[0:length(Num)-1]'*thetac) ...
./(Den*exp(-j*[0:length(Den)-1]'*thetac));
Hc=20*log10(abs(Hc))
Ap1
As1
