clc;
clear all
syms X
z=0.475*1j;
Vs=1.05;
A=1;
a1=real(A); a2=imag(A);
A=a1+a2*1j;
B=z;
b1=real(B); b2=imag(B);
C=0;
D=A;
fi=acos(0.9);
K1=a1*(b2-b1*tan(fi))+a2*(b1+b2*tan(fi));
K2=a1*(b1+b2*tan(fi))+a2*(b1-b2*tan(fi));
deltarcrit=(pi/4)+0.5*atan(K2/-K1);
Vrcrit=Vs/(2*(a1*cos(deltarcrit)+a2*sin(deltarcrit)));
K3=b1*cos(deltarcrit)+b2*sin(deltarcrit);
K4=a1*cos(deltarcrit)+a2*sin(deltarcrit);
Prcrit=((Vs^2)*(2*K3*K4-(a1*b1+a2*b2)))/((b1^2+b2^2)*4*K4);
Vr=[];
for P=60:1:350
        Qr=P*tan(fi);
        P1=a1^2+a2^2;
        P2=2*P*(a1*b1+a2*b2)+2*Qr*(a1*b2+a2*b1)-Vs^2;
        P3=((b1+b2).^2)*(P^2+Qr^2);
        equation=P1*(X^2)+P2*X+P3;
        these_roots = roots([P1 P2 P3]);
        mask = any(imag(these_roots) ~= 0,2);
        these_roots(mask,:) = nan;
        Vr=[Vr these_roots];
  end
        Pr=(60:1:350);
%plot(Pr,Vr.')
display(Prcrit)
display(Vrcrit)