clear all;format compact;format long;clc;
%close all;
Ap=1;As=43;
fp=2500;fs=1500;fsampling=8000;
excess='pass';

%f (Hz) to theta (rad)
thetap=2*pi*fp/fsampling
thetas=2*pi*fs/fsampling

%Prewarp (theta to w (rad/s))
wp=2*fsampling*tan(thetap/2)
ws=2*fsampling*tan(thetas/2)

%Backward Transformation
Wp=1
Ws=wp/ws

%Order n
n=acosh(sqrt((10^(0.1*As)-1)/(10^(0.1*Ap)-1)))/acosh(Ws)
n=ceil(n)

%Recalculation of Ap and As
Ap1=Ap;
As1=As;
if excess == lower('pass')
    Ap1=10*log10(1+(10^(0.1*As1)-1)/cosh(n*acosh(Ws))^2);
    epsilon1=1/sqrt(10^(0.1*As1)-1)    ;
else
    As1=10*log10(1+(10^(0.1*Ap1)-1)*cosh(n*acosh(Ws))^2);
end   
Ap1
As1
ep1=1/sqrt(10^(0.1*As1)-1)
%Normalized lowpass poles and zeros
%Parameter a
a=(1/n)*asinh(sqrt(10^(0.1*As1)-1)) 
k=0:ceil((n-3)/2);
Sa=-1./(sinh(a)*sin(pi*(2*k+1)/(2*n))+j*cosh(a)*cos(pi*(2*k+1)/(2*n)));
Sa=[Sa,conj(Sa)];
Sa=cplxpair(Sa);
Sb=reshape(Sa,[2,ceil((n-1)/2)]);
S=Sa;
if rem(n,2)==1 S=[-1/sinh(a),Sa]; end
S
k1=0:floor((n-2)/2);
SZa=j./cos(pi*(2*k1+1)/2/n);
SZa=[SZa,conj(SZa)];
SZa=cplxpair(SZa);
SZb=reshape(SZa,[2,ceil((n-1)/2)]);
SZ=SZa;
if rem(n,2)==1 SZ=[inf,SZa]; end
SZ

%Frequency transformation  
wo=ws;
s=ws./S
sz=ws./SZ

%Digital poles and zeros
pa=(2*fsampling+s)./(2*fsampling-s);
za=(2*fsampling+sz)./(2*fsampling-sz);
pa=cplxpair(pa);
za=cplxpair(za);
if rem(n,2)==1 p=[pa,0];z=[za,0]; else p=pa;z=za; end
p
z
p1=reshape(p,[2,floor((n+1)/2)])
z1=reshape(z,[2,floor((n+1)/2)])
DenSOS=[ones(1,floor((n+1)/2));-sum(p1);prod(p1)]
NumSOS=[ones(1,floor((n+1)/2));-sum(z1);prod(z1)]

%b0
b0=abs([1 -1 1]*DenSOS)./abs([1 -1 1]*NumSOS)
if rem(n,2)==0
    b0=10^(-0.05*Ap1/(n/2))*b0;
end
b0
Den=poly(p);
Num=poly(z);
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
text(0.85,0.1,nstr,'FontSize',12);
title('Digital Pole-Zero Diagram');hold off;
M=5000;
theta=0:pi/M:pi;
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
xlabel('\theta/\pi');ylabel('\angleH(\theta) (deg)');
title('Phase response');
H0=abs(H(1));
H1=abs(H(length(H)));
thetao=2*atan(wo/(2*fsampling));
thetac=[thetap,thetas,thetao];
Qc=exp(-j*[0:2]'*thetac);
Hc=singleb0*Num*exp(-j*[0:length(Num)-1]'*thetac) ...
./(Den*exp(-j*[0:length(Den)-1]'*thetac));
Hc=20*log10(abs(Hc));
Ap1
As1