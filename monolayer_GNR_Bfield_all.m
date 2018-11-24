% The band structure and NEGF for monolayer graphene Nanoribbons
% The major difference between graphene with infinite boundry-peoridic
% condictions: the unit cell is rectangular.

%% Device Setup
clear all; clc;

%define constants
a = 1.42;
q = 1.6e-19; 
hbar = 1.06e-34;
h = 2*pi*hbar;
m = 9.1e-31;
zplus=i*1e-3; 
%t has been given as a parameter
t = -3;

AtoM=1e-10;


f1 =1; f2 = 0;
fprintf('for GNR\n');
device = 1;

fprintf(1,'Enter the configuration of device you want to simulate Armchair/Zigzag\n')
fprintf(1,'1 for Armchair\n');
fprintf(1,'2 for Zigzag\n');
type = input('Valid inputs are 1 or 2 : ');

if (~((type == 1)||(type == 2)))
    type = 1;
end

%Points in GNR width
NW = input('Enter number of Unit Cells in the width : ');
if isempty(NW)
    NW = 12;
end

%Points in GNR length
NL = 1;
%Points in GNR unit cell
NU = 4;

%Define any scatterer as a delta potential barrier
%No scatterer

%Define phase breaking potential
D = 1e-12;
%Convergence criteria
errmax = 1e-3;
%Total H will be of LxWxU size
NT = NL*NW*NU;

fprintf(1,'Simulating...\n');

if (type == 1)
    fprintf(1,'Type: Armchair\n');
else
    fprintf(1,'Type: Zigzag\n');
end

%% NEGF Setup

%Define Hamiltonians and alpha, beta matrices
fprintf(1,'Hamiltonian Setup...please wait\n');
alphau = t*diag(ones(1,NU-1),1)+t*diag(ones(1,NU-1),-1); 
betaw = zeros(NU);
betal = zeros(NU);

%Armchair Structure for GNR
if ((type == 1)&&(device == 1))
    betaw(2,1) = t; betaw(3,4) = t;
    betal(4,1) = t;
elseif ((type == 2)&&(device == 1))

%Zigzag Structure for GNR
    betaw(4,1) = t;
    betal(2,1) = t; betal(3,4) = t;
elseif ((type == 1)&&(device == 2))
end


Bz = input('Magnetic field Bz: ');

% ignore z positions for calculation
A1 = [0,0];
B1 = [1/2,-sqrt(3)/2]*a;
A1_p = [3/2,-sqrt(3)/2]*a;
B1_p = [2,0]*a;

Xlist = [0,1/2,3/2,2]*a;
Ylist = [0,-sqrt(3)/2,-sqrt(3)/2,0]*a;

% dcell_alpha = [0,0];
% dcell_beta1 = [3,0]*a;
% dcell_beta2 = [0,sqrt(3)]*a;

alpha = kron(diag(ones(1,NW)),alphau)+kron(diag(ones(1,NW- 1),1),betaw)+kron(diag(ones(1,NW-1),-1),betaw');
beta = kron(diag(ones(1,NW)),betal);

l1 = sqrt(3)*a;
l2 = 3*a;

if (type == 1)
    
    fprintf('using gauge1 \n');
    
    % Modify alpha
    for i = 1:NU*NW
        for j = 1:NU*NW
            ii = floor((i-0.5)/NU);
            jj = floor((j-0.5)/NU);
            y_i = Ylist(i-ii*NU)+ii*l1;
            y_j = Ylist(j-jj*NU)+jj*l1;
            x_i = Xlist(i-ii*NU);
            x_j = Xlist(j-jj*NU);
            %dx = dcell_alpha(1);
            alpha(i,j) = alpha(i,j)*exp(1i*q*(-Bz)*(y_i+y_j)/2*(x_i-x_j)/hbar*AtoM^2);
        end
    end
    
    % Modify beta
    for i = 1:NU*NW
        for j = 1:NU*NW
            ii = floor((i-0.5)/NU);
            jj = floor((j-0.5)/NU);
            y_i = Ylist(i-ii*NU)+ii*l1;
            y_j = Ylist(j-jj*NU)+jj*l1;
            x_i = Xlist(i-ii*NU);
            x_j = Xlist(j-jj*NU);
            %dx = dcell_beta1(1);
            beta(i,j) = beta(i,j)*exp(1i*q*(-Bz)*(y_i+y_j)/2*(x_i-x_j+l2)/hbar*AtoM^2);
        end
    end
    
else
    
    fprintf('using gauge2 \n');
    
    % Modify alpha
    for i = 1:NU*NW
        for j = 1:NU*NW
            ii = floor((i-0.5)/NU);
            jj = floor((j-0.5)/NU);
            y_i = Ylist(i-ii*NU);
            y_j = Ylist(j-jj*NU);
            x_i = Xlist(i-ii*NU)+ii*l2;
            x_j = Xlist(j-jj*NU)+jj*l2;
            %dy = dcell_alpha(2);
            alpha(i,j) = alpha(i,j)*exp(1i*q*(-Bz)*(x_i+x_j)/2*(y_i-y_j)/hbar*AtoM^2);
        end
    end
    
    % Modify beta
    for i = 1:NU*NW
        for j = 1:NU*NW
            ii = floor((i-0.5)/NU);
            jj = floor((j-0.5)/NU);
            y_i = Ylist(i-ii*NU);
            y_j = Ylist(j-jj*NU);
            x_i = Xlist(i-ii*NU)+ii*l2;
            x_j = Xlist(j-jj*NU)+jj*l2;
            %dy = dcell_beta2(2);
            beta(i,j) = beta(i,j)*exp(1i*q*(-Bz)*(x_i+x_j)/2*(y_i-y_j+l1)/hbar*AtoM^2);
        end
    end

end


% Hamiltonian: 
H = kron(diag(ones(1,NL)),alpha)+kron(diag(ones(1,NL- 1),1),beta)+kron(diag(ones(1,NL-1),-1),beta');



%Define Energy grid for calculation of Transmission
E = linspace(-1,1,51);
%Define the matrices for NEGF
T = zeros(1,length(E));
g1 = inv(E(1)*eye(NW*NU)-alpha);
g2 = inv(E(1)*eye(NW*NU)-alpha);
Es = zeros(NT);
Esin = zeros(NT);
fprintf(1,'Done\n');

%% NEGF Calculations
fprintf(1,'Running NEGF...please wait\n')
fprintf(1,'Energy range = %d eV to %d eV\n',min(E),max(E)); 
fprintf(1,'%s\n',['[',blanks(length(E)),']']);
fprintf(1,' ');

%Run for each energy
for k = 1:length(E)
%Calculate surface GFs self consistently
%Calculate for beta'
    err = 100;
    while(err>errmax)
        g1new = inv((E(k)+zplus)*eye(NW*NU)-alpha-beta'*g1*beta); 
        err = (sum(sum(abs(g1new-g1))))/(sum(sum(abs(g1new+g1)))); 
        g1 = g1new;
    end
    sigma1  = beta'*g1*beta;
%Calculate for beta
    err = 100;
    while(err>errmax)
        g2new = inv((E(k)+zplus)*eye(NW*NU)-alpha-beta*g2*beta'); 
        err = (sum(sum(abs(g2new-g2))))/(sum(sum(abs(g2new+g2)))); 
        g2 = g2new;
    end
    sigma2  = beta*g2*beta';
%Calculate self energy matrices
    E1 = kron(diag([1 zeros(1,NL-1)]),sigma1);
    E2 = kron(diag([zeros(1,NL-1) 1]),sigma2);
%Calculate broadening
    G1 = 1i*(E1-E1'); 
    G2 = 1i*(E2-E2');
%Calculate G, Gn, A, T for coherent transport
%     G = inv((E(k)+1i*etaplus)*eye(NT)-H-E1-E2);
%     T(k) = real(trace(G1*G*G2*G'));
%     Gn = f1*(G*G1*G')+f2*(G*G2*G');
%     A = i*(G-G');
%Calculate G, Gn, A, T self consistently including the phase breaking processes
    err = 100;
    while(err>errmax)
        G = inv((E(k)+zplus)*eye(NT)-H-E1-E2-Es);
        Esnew = D*G;
        err = sum(sum(abs(Esnew-Es)))/sum(sum(abs(Esnew+Es))); 
        Es = Esnew;
    end
    % % % % % % % %
%      err = 100;
%      while(err>errmax)
%          Gn = f1*(G*G1*G')+f2*(G*G2*G')+(G*Esin*G');
%          Esinnew = D*Gn;
%          err = sum(sum(abs(Esinnew-Esin)))/sum(sum(abs(Esinnew+Esin))); Esin = Esinnew;
%      end
%      A=i*(G-G');
    T(k) = real(trace(G1*G*G2*G'));
    fprintf(1,'|');
end

fprintf(1,'\n');
fprintf(1,'Done\n');
%Calculate Subbands
ka = linspace(-pi,pi,501);
eigE = zeros(length(ka),NW*NU); 
fprintf(1,'Calculating Sub-bands...please wait\n'); 
for k = 1:length(ka)
    [V,D] = eig(alpha + beta*exp(1i*ka(k))+ beta'*exp(-1i*ka(k)));
    eigE(k,:) = diag(D);
end

fprintf(1,'Done\n');


%% Plot Figures
%Figure 1: Transmission and sub-band comparison figure(1);
subplot(1,2,1);
plot(T,E);
xlabel('T(E)');
ylabel('E[eV]');
title('Energy v Transmission');
xlim([0 max(T)+0.5]);

ylim([min(E) max(E)]);
subplot(1,2,2);
plot(ka/pi,eigE);
xlabel('ka/\pi');
ylabel('E[eV]');
title('Sub-bands');
ylim([min(E) max(E)]);

%Figure 2: Sub-bands
figure(2);
plot(ka/pi,eigE,'b');
xlabel('ka/\pi');
ylabel('E[eV]');
title('Sub-bands');