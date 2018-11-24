clear all


deltaAB=46;
gamma0=3100;
gamma1=390;
%gamma1prime = 320; %320;
gamma2=-28/2;
gamma3=315; % control the touching point. 
gamma4=41; %40;
gamma5=50/2; %60/2;

xi=1;

coolor='-k';

a=1.42*sqrt(3);

G1=2*pi/a*[1,-1/sqrt(3)];
G2=2*pi/a*[0,2/sqrt(3)];

GA=(-G1*(2/3)-G2*(1/3))/40;
K2=(-G1*(1/3)+G2*(1/3))/40;
K1=[0,0,0];

NK=80;

KX1=linspace(GA(1,1),K1(1,1),NK) ;
KY1=linspace(GA(1,2),K1(1,2),NK) ;

KX2=linspace(K1(1,1),K2(1,1),NK) ;
KY2=linspace(K1(1,2),K2(1,2),NK) ;


Kp=[KX1,KX2;KY1,KY2];


for kpp=1:size(Kp,2)

k=Kp(:,kpp)';

HG_single=[-gamma2,(sqrt(3)*a/2)*gamma0*(xi*k(1,1)-1i*k(1,2));(sqrt(3)*a/2)*gamma0*(xi*k(1,1)+1i*k(1,2)),-gamma5+deltaAB];

HG_bilayer=[gamma2,sqrt(2)*(sqrt(3)*a/2)*gamma3*(xi*k(1,1)+1i*k(1,2)),-sqrt(2)*(sqrt(3)*a/2)*gamma4*(xi*k(1,1)-1i*k(1,2)),(sqrt(3)*a/2)*gamma0*(xi*k(1,1)-1i*k(1,2));
    sqrt(2)*(sqrt(3)*a/2)*gamma3*(xi*k(1,1)-1i*k(1,2)),0,(sqrt(3)*a/2)*gamma0*(xi*k(1,1)+1i*k(1,2)),-sqrt(2)*(sqrt(3)*a/2)*gamma4*(xi*k(1,1)+1i*k(1,2));
    -sqrt(2)*(sqrt(3)*a/2)*gamma4*(xi*k(1,1)+1i*k(1,2)),(sqrt(3)*a/2)*gamma0*(xi*k(1,1)-1i*k(1,2)),deltaAB,sqrt(2)*gamma1;
    (sqrt(3)*a/2)*gamma0*(xi*k(1,1)+1i*k(1,2)),-sqrt(2)*(sqrt(3)*a/2)*gamma4*(xi*k(1,1)-1i*k(1,2)),sqrt(2)*gamma1,gamma5+deltaAB]

T2=[0, 0,0,0;0,0,0,0];

T2_t=[0,0;0,0;0,0;0,0];

HG3=[HG_single,T2;T2_t,HG_bilayer];

Egval_single(:,kpp)=sort(real(eig(HG_single)));

Egval_bilayer(:,kpp)=sort(real(eig(HG_bilayer)));

Egval(:,kpp)=sort(real(eig(HG3)));


%HG_trilayer= 


end

%kp_num=[linspace(0,1,NK),linspace(1,2,NK)];
kp_num=[linspace(-0.0426,0.0426,NK*2)];
figure;
% for i=3:4 % size(Egval,1)
%     plot(kp_num,Egval(i,:),coolor,'linewidth',[5]);
%     hold on
% end

% for i=3:6 %1:size(Egval,1)
%     plot(kp_num,Egval(i,:),coolor,'linewidth',[3]);
%     hold on
% end

for i=1:2 % size(Egval,1)
    plot(kp_num,Egval_single(i,:),coolor,'linewidth',[2],'Color','b');
    hold on
end

for i=1:4 % size(Egval,1)
    plot(kp_num,Egval_bilayer(i,:),coolor,'linewidth',[2],'Color','r');
    hold on
end


axis([-0.03  0.03 -600 700]);
%axis([kp_num(1,1)  kp_num(1,size(kp_num,2)) -20 20]);