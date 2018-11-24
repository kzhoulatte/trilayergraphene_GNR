clear all

N=30; %number of 

%Pablo&Abanin
% Gamma0=3100;
% Gamma1=390;
% Gamma3=315;
% Gamma4=41;
% Gamma2=-28;
% Gamma5=50;
% delta=46;
% Delta2=0;

%Koshino
%Gamma0=3000;
%Gamma1=400;
%Gamma3=300;
%Gamma4=40;
%Gamma2=-20;
%Gamma5=40; 
%delta=50;
%Delta2=0;

%PRL
Gamma0=3230;
Gamma1=310;
Gamma3=300;
Gamma4=40;
Gamma2=-9800/Gamma1;
Gamma5=10; 
delta=1.5+(Gamma5-Gamma2)/2;
Delta2=1.8;


v3=Gamma3/Gamma0;
v4=Gamma4/Gamma0;
Deltaext=0;
initial=0.05;
M = 4; inc = 0.01; %maximum B - field and step
G = zeros(M,6*N+9+1);
GKp = zeros(M,6*N+9+1);
counter = 0; 

for B = initial:inc:M %B field values
    H = zeros(6*N+9); %matrix initialization
    HKp = zeros(6*N+9); 
    vel= (3/2*Gamma0*0.142);
    E0 = ((sqrt(2)*vel)/25.66)*sqrt(B);
    alpha3 = E0*v3;
    alpha4 = E0*v4;
    U1=Gamma2/2+Delta2;
    U4=-2*Delta2;
    U3=delta-2*Delta2;
    U2=Gamma5/2+delta+Delta2;
    U5=Delta2-Gamma2/2;
    U6=-Gamma5/2+delta+Delta2;
    counter = counter + 1;
  
   % valley K 
   % Bilayer block
    
    for m = 3:1:N+3        
        H(m,m) = E0*sqrt(m-2)+(U1+U2)/2;
        H(m,N+4+(m-3))=(U1-U2)/2;
        H(N+4+(m-3),m)=(U1-U2)/2;
        H(m,4*N+8+(m-3))=Deltaext;
        H(4*N+8+(m-3),m)=Deltaext;
        if m < N+3
        H(m,2*N+5+(m-2))=sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-2)+sqrt(m-1)));
        H(2*N+5+(m-2),m)=sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-2)+sqrt(m-1)));
        H(m,3*N+6+(m-2))=sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-2)-sqrt(m-1)));
        H(3*N+6+(m-2),m)=sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-2)-sqrt(m-1)));
        end
        if m > 4
        H(m,2*N+5+(m-5))=sqrt(2)*((1/2)*sqrt((m-3))*alpha3);
        H(2*N+5+(m-5),m)=sqrt(2)*((1/2)*sqrt((m-3))*alpha3);
        H(m,3*N+6+(m-5))=-sqrt(2)*((1/2)*sqrt((m-3))*alpha3);
        H(3*N+6+(m-5),m)=-sqrt(2)*((1/2)*sqrt((m-3))*alpha3); 
        end
      end
    
    for m=N+4:1:2*N+4
        H(m,m) = -E0*sqrt(m-(N+3))+(U1+U2)/2;
        H(m,5*N+9+(m-(N+4)))=Deltaext;
        H(5*N+9+(m-(N+4)),m)=Deltaext;
        if m < 2*N+4
        H(m,2*N+5+(m-(N+3)))=-sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-(N+2))-sqrt(m-(N+3))));
        H(2*N+5+(m-(N+3)),m)=-sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-(N+2))-sqrt(m-(N+3))));
        H(m,3*N+6+(m-(N+3)))=-sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-(N+2))+sqrt(m-(N+3))));
        H(3*N+6+(m-(N+3)),m)=-sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-(N+2))+sqrt(m-(N+3))));
        end
        if m > N+5
        H(m,2*N+5+(m-(N+6)))=sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3);
        H(2*N+5+(m-(N+6)),m)=sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3);
        H(m,3*N+6+(m-(N+6)))=-sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3);
        H(3*N+6+(m-(N+6)),m)=-sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3); 
        end
    end
    
    for m = 2*N+5:1:3*N+5
        H(m,m) = E0*sqrt(m-(2*N+4))+(U3+U4)/2;
        H(m,3*N+6+(m-(2*N+5)))=(U3-U4)/2;
        H(3*N+6+(m-(2*N+5)),m)=(U3-U4)/2;
    end
    
    for m = 3*N+6:1:4*N+6    
        H(m, m) = -E0*sqrt(m-(3*N+5))+(U3+U4)/2;
    end
    
    H(2,2*N+5)=Gamma1-alpha4;
    H(2,3*N+6)=Gamma1+alpha4;
    H(2*N+5,2)=Gamma1-alpha4;
    H(3*N+6,2)=Gamma1+alpha4;
    H(1,4)=alpha3;
    H(4,1)=alpha3;
    H(1,N+5)=alpha3;
    H(N+5,1)=alpha3;
    H(1,1)=U4;
    H(2,2)=U2;
    H(2,4*N+7)=Deltaext;
    H(4*N+7,2)=Deltaext;

    %Monolayer block
    
    for m=4*N+8:1:5*N+8
        H(m,m)=E0*sqrt(m-(4*N+8)+1)+(U5+U6)/2;
        H(m,5*N+9+m-(4*N+8))=(U5-U6)/2;
        H(5*N+9+m-(4*N+8),m)=(U5-U6)/2;
    end 
    for m=5*N+9:1:6*N+9
        H(m,m)=-E0*sqrt(m-(5*N+9)+1)+(U5+U6)/2;
    end
    H(4*N+7,4*N+7)=U6;
    
   % valley Kp 
   % Bilayer block
    
    for m = 3:1:N+3        
        HKp(m, m) = -E0*sqrt(m -2)+(U4+U3)/2;
        HKp(m,N+4+(m-3))=(U4-U3)/2;
        HKp(N+4+(m-3),m)=(U4-U3)/2;
        if m < N+3
        HKp(m,2*N+5+(m-2))=sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-2)+sqrt(m-1)));
        HKp(2*N+5+(m-2),m)=sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-2)+sqrt(m-1)));
        HKp(m,3*N+6+(m-2))=sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-2)-sqrt(m-1)));
        HKp(3*N+6+(m-2),m)=sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-2)-sqrt(m-1)));
        end
        if m > 4
        HKp(m,2*N+5+(m-5))=sqrt(2)*((1/2)*sqrt((m-3))*alpha3);
        HKp(2*N+5+(m-5),m)=sqrt(2)*((1/2)*sqrt((m-3))*alpha3);
        HKp(m,3*N+6+(m-5))=-sqrt(2)*((1/2)*sqrt((m-3))*alpha3);
        HKp(3*N+6+(m-5),m)=-sqrt(2)*((1/2)*sqrt((m-3))*alpha3); 
        end
    end
    
    for m=N+4:1:2*N+4
        HKp(m, m) = E0*sqrt(m-(N+3))+(U4+U3)/2;
        if m < 2*N+4
        HKp(m,2*N+5+(m-(N+3)))=-sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-(N+2))-sqrt(m-(N+3))));
        HKp(2*N+5+(m-(N+3)),m)=-sqrt(2)*(Gamma1/2+(alpha4/2)*(sqrt(m-(N+2))-sqrt(m-(N+3))));
        HKp(m,3*N+6+(m-(N+3)))=-sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-(N+2))+sqrt(m-(N+3))));
        HKp(3*N+6+(m-(N+3)),m)=-sqrt(2)*(Gamma1/2-(alpha4/2)*(sqrt(m-(N+2))+sqrt(m-(N+3))));
        end
        if m > N+5
        HKp(m,2*N+5+(m-(N+6)))=sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3);
        HKp(2*N+5+(m-(N+6)),m)=sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3);
        HKp(m,3*N+6+(m-(N+6)))=-sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3);
        HKp(3*N+6+(m-(N+6)),m)=-sqrt(2)*((1/2)*sqrt((m-(N+4)))*alpha3); 
        end
    end
    
    for m = 2*N+5:1:3*N+5
        HKp(m,m) =  -E0*sqrt(m-(2*N+4))+(U2+U1)/2;
        HKp(m,3*N+6+(m-(2*N+5)))=(U2-U1)/2;
        HKp(3*N+6+(m-(2*N+5)),m)=(U2-U1)/2;
        HKp(m,4*N+8+(m-(2*N+5)))=Deltaext;
        HKp(4*N+8+(m-(2*N+5)),m)=Deltaext;
    end
    
    for m = 3*N+6:1:4*N+6    
        HKp(m, m) = E0*sqrt(m-(3*N+5))+(U2+U1)/2;
        HKp(m,5*N+9+(m-(3*N+6)))=Deltaext;
        HKp(5*N+9+(m-(3*N+6)),m)=Deltaext;
    end
    
    HKp(2,2*N+5)=-(Gamma1+alpha4);
    HKp(2,3*N+6)=-(Gamma1-alpha4);
    HKp(2*N+5,2)=-(Gamma1+alpha4);
    HKp(3*N+6,2)=-(Gamma1-alpha4);
    HKp(1,4)=alpha3;
    HKp(4,1)=alpha3;
    HKp(1,N+5)=alpha3;
    HKp(N+5,1)=alpha3;
    HKp(1,1)=U1;
    HKp(2,2)=U3;
    HKp(1,4*N+7)=Deltaext;
    HKp(4*N+7,1)=Deltaext;

    %Monolayer block
    
    for m=4*N+8:1:5*N+8
        HKp(m,m)=-E0*sqrt(m-(4*N+8)+1)+(U5+U6)/2;
        HKp(m,5*N+9+m-(4*N+8))=-(U5-U6)/2;
        HKp(5*N+9+m-(4*N+8),m)=-(U5-U6)/2;
    end 
    for m=5*N+9:1:6*N+9
        HKp(m,m)=E0*sqrt(m-(5*N+9)+1)+(U5+U6)/2;
    end
    HKp(4*N+7,4*N+7)=U5;
    
    %Diagnolize
    [Q, D] = eig(H);
    
    %Filling a column in the G (graph) matrix
    for m = 1:6*N+9
        evalue(m) =D(m,m);
        G(counter, m + 1) = D(m, m);
    end
    %Filling a B-field column in the G (graph) matrix
    G(counter, 1) = B;
    
    [QKp, DKp] = eig(HKp);
    
    %Filling a column in the G (graph) matrix
    for m = 1:6*N+9
        evalue(m) =DKp(m,m);
        GKp(counter, m + 1) = DKp(m, m);
    end
    %Filling a B-field column in the G (graph) matrix
    GKp(counter, 1) = B;
    
end
evalue;
figure
%plot(G(:,1), G(:,2), 'k', 'LineWidth', 0.5)
% Set the axis limits
axis([initial M -60 60])

hold on
    for i = 2:3*N+6
        plot(G(:,1), G(:,i), 'b', 'LineWidth', 0.5)
    end
%     for i = 3*N+5:3*N+7
%         plot(G(:,1), G(:,i), 'k', 'LineWidth', 0.5)
%     end
    for i = 3*N+7:6*N+10
        plot(G(:,1), G(:,i), 'b', 'LineWidth', 0.5)
    end
hold on
    for i = 2:3*N+6
        plot(GKp(:,1), GKp(:,i), 'b--', 'LineWidth', 1.5)
    end
%     for i = 3*N+5:3*N+7
%         plot(G(:,1), G(:,i), 'k', 'LineWidth', 0.5)
%     end
    for i = 3*N+7:6*N+10
        plot(GKp(:,1), GKp(:,i), 'b--', 'LineWidth', 1.5)
    end
    
hold off

% Turn on the grid
grid on

% Add title and axis labels
title('Trilayer LL spectrum')
xlabel('B Field, T')
ylabel('Energy, meV')