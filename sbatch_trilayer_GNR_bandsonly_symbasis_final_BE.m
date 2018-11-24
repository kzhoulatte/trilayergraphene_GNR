% The band structure and NEGF for monolayer graphene Nanoribbons
% The major difference between graphene with infinite boundry-peoridic
% condictions: the unit cell is rectangular.

%% Device Setup
clear all; clc;
warning('off','all');

device = 1;
type = 1;
NW = 10;
dE = 0;
Bz = 10;
Basis = 1;

file_out = fopen('TNR_output.txt','w');
kps_out = fopen('TNR_kps.txt','w');

if Basis==1
    eigEBL_out = fopen('TNR_eigEBL.txt','w');
    eigESL_out = fopen('TNR_eigESL.txt','w');
else
    eigE_out = fopen('TNR_eigE.txt','w');
end

tic;

a = 1.42;
%define constants
q = 1.6e-19; 
hbar = 1.06e-34;
h = 2*pi*hbar;
m = 9.1e-31;
zplus=i*1e-3; 
%t has been given as a parameter
Ea = 0;
Eb = 0;

AtoM = 1e-10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARAMs from Supeng Thesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% delta1 = 0;
% delta2 = 0.00451;
% deltaE = 0.040;
% t = 2.550;
% %gamma0 
% g0= t;
% %gamma1 
% g1 = 0.350;
% %gamma2 
% g2 = -0.0234/2;
% %gamma3 
% g3 = 0;%0.280;
% %gamma4 
% g4 = 0.130;
% %gamma5 
% g5 = 0.034/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARAMs from Paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta1 = 0;
delta2 = 0.0;
deltaE = 0.046;

%t = 3.1;
%gamma0 
g0= 3.1;
%gamma1 
g1 = 0.39;
%gamma2 
g2 = -0.028/2;
%gamma3 
g3 = 0.315;
%g3_prime = 0;%.315;
%gamma4 
g4 = 0.041;
%gamma5 
g5 = 0.05/2;

ymin = -10;
ymax = 10;

f1 =1; f2 = 0;
fprintf(file_out,'for GNR\n');


if (~((type == 1)||(type == 2)))
    type = 1;
end

%Points in GNR width

if isempty(NW)
    NW = 3;
end

%Points in GNR length
NL = 1;
%Points in GNR unit cell
NU = 12;

%Define any scatterer as a delta potential barrier
%No scatterer

%Define phase breaking potential
D = 1e-12;
%Convergence criteria
errmax = 1e-3;
%Total H will be of LxWxU size
NT = NL*NW*NU;

fprintf(file_out,'Simulating...\n');

if (type == 1)
    fprintf(file_out,'Type: Armchair\n');
else
    fprintf(file_out,'Type: Zigzag\n');
end



%% NEGF Setup

%Define Hamiltonians and alpha, beta matrices
fprintf(file_out,'Hamiltonian Setup...please wait\n');

%%%%%%%%%%%
%Calculate Alpha, Beta from alpha, beta etc. 
%%%%%%%%%%%

alpha= [Ea+delta1+delta2+dE,    g0, 0,   g0,    -g4,g3, -g4,g3,    g2, 0, 0, 0;
        g0,Eb+deltaE+delta1+delta2+dE, 0, 0,     g1,-g4,   0,-g4,    0, g5, 0, 0;
         0, 0,       Ea+delta1+delta2+dE,g0,     0,  0, -g4,g3,    0,  0,g2, 0;
        g0, 0,g0,Eb+deltaE+delta1+delta2+dE,      0, 0,  g1,-g4,    0,  0, 0,g5;
       -g4, g1,0, 0,                        Ea+deltaE-2*delta2, g0,0, g0,-g4, g1,0, 0;
       g3,-g4, 0,0,                        g0, Eb-2*delta2,   0,  0,g3,-g4, 0,0;
       -g4,  0,-g4,g1,0,0,Ea+deltaE-2*delta2,g0,                    -g4,0,-g4,g1;
       g3,-g4,g3,-g4,g0,0,g0,Eb-2*delta2,                           g3,-g4,g3,-g4;
        g2, 0,0,0,-g4,g3,-g4,g3,Ea-delta1+delta2-dE,g0,0,g0;
         0, g5,0,0,g1,-g4,0,-g4,g0,Eb+deltaE-delta1+delta2-dE,0,0;
         0, 0,g2,0,0,0,-g4,g3,0,0,Ea-delta1+delta2-dE,g0;
         0, 0,0,g5,0,0,g1,-g4,g0,0,g0,Eb+deltaE-delta1+delta2-dE];
    
    
%t = zeros(12);
t = [dE,0,0,0, 0,0,0,0, 0,0,0,0;
     0,dE,0,0, 0,0,0,0, 0,0,0,0;
     0,g0,dE,0, -g4,g3,0,0, 0,0,0,0;
     0,0,0,dE, 0,-g4,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,g0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0,  0,0,0,0, -dE,0,0,0;
     0,0,0,0, 0,0,0,0, 0,-dE,0,0;
     0,0,0,0, -g4,g3,0,0, 0,g0,-dE,0;
     0,0,0,0, 0,-g4,0,0, 0,0,0,-dE];

beta0 = [dE,0,0,0, 0,0,0,0, 0,0,0,0;
     g0,dE,0,0, 0,-g4,0,0, 0,0,0,0;
     0,0,dE,g0, 0,0,-g4,0, 0,0,0,0;
     0,0,0,dE, 0,0,0,0, 0,0,0,0;
     -g4,0,0,0, 0,g0,0,0, -g4,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     g3,0,0,-g4, 0,0,g0,0, g3,0,0,-g4;
     0,0,0,0, 0,0,0,0, -dE,0,0,0;
     0,0,0,0, 0,-g4,0,0, g0,-dE,0,0;
     0,0,0,0, 0,0,-g4,0, 0,0,-dE,g0;
     0,0,0,0, 0,0,0,0, 0,0,0,-dE];
 
beta1 = [dE,0,0,0, 0,0,0,0, 0,0,0,0;
     0,dE,0,0, 0,0,0,0, 0,0,0,0;
     0,0,dE,0, 0,0,0,0, 0,0,0,0;
     0,0,0,dE, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, -dE,0,0,0;
     0,0,0,0, 0,0,0,0, 0,-dE,0,0;
     0,0,0,0, 0,0,0,0, 0,0,-dE,0;
     0,0,0,0, 0,0,0,0, 0,0,0,-dE];
 
beta2 = [dE,0,0,0, 0,0,0,0, 0,0,0,0;
     0,dE,0,0, 0,0,0,0, 0,0,0,0;
     0,0,dE,0, 0,g3,0,0, 0,0,0,0;
     0,0,0,dE, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, 0,0,0,0;
     0,0,0,0, 0,0,0,0, -dE,0,0,0;
     0,0,0,0, 0,0,0,0, 0,-dE,0,0;
     0,0,0,0, 0,g3,0,0, 0,0,-dE,0;
     0,0,0,0, 0,0,0,0, 0,0,0,-dE];
 

M = [1/sqrt(2),0,0,0, 0,0,0,0, -1/sqrt(2),0,0,0;
    0,0,1/sqrt(2),0, 0,0,0,0, 0,0,-1/sqrt(2),0;
    0,1/sqrt(2),0,0, 0,0,0,0, 0,-1/sqrt(2),0,0;
    0,0,0,1/sqrt(2), 0,0,0,0, 0,0,0,-1/sqrt(2);
    1/sqrt(2),0,0,0, 0,0,0,0, 1/sqrt(2),0,0,0;
    0,0,1/sqrt(2),0, 0,0,0,0, 0,0,1/sqrt(2),0;
    0,0,0,0, 0,1,0,0, 0,0,0,0;
    0,0,0,0, 0,0,0,1, 0,0,0,0;
    0,0,0,0, 1,0,0,0, 0,0,0,0;
    0,0,0,0, 0,0,1,0, 0,0,0,0;
    0,1/sqrt(2),0,0, 0,0,0,0, 0,1/sqrt(2),0,0;
    0,0,0,1/sqrt(2), 0,0,0,0,0,0,0, 1/sqrt(2)];

%alpha = M*alpha*M';
%t = M*t*M';
%beta0 = M*beta0*M';
%beta1 = M*beta1*M';
%beta2 = M*beta2*M';


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% If added constant Bz magnetic field, change hamiltonian components
% accordingly: 



% ignore z positions for calculation
A1 = [0,0];
B1 = [-1/2,sqrt(3)/2]*a;
A1_p = [3/2,sqrt(3)/2]*a;
B1_p = [1,0]*a;
A2 = [-1/2,sqrt(3)/2]*a;
B2 = [-1,0]*a;
A2_p = [1,0]*a;
B2_p = [1/2,sqrt(3)/2]*a;
A3 = [0,0];
B3 = [-1/2,sqrt(3)/2]*a;
A3_p = [3/2,sqrt(3)/2]*a;
B3_p = [1,0]*a;

Xlist = [0,-1/2,3/2,1,-1/2,-1,1,1/2,0,-1/2,3/2,1]*a;
Ylist = [0,sqrt(3)/2,sqrt(3)/2,0,sqrt(3)/2,0,0,sqrt(3)/2,0,sqrt(3)/2,sqrt(3)/2,0]*a;

% dcell_alpha = [0,0];
% dcell_t = [3,0]*a;
% dcell_beta0 = [0,sqrt(3)]*a;
% dcell_beta1 = [-3,sqrt(3)]*a;
% dcell_beta2 = [3,sqrt(3)]*a;

% Alpha_1 = [alpha,beta0',zeros(12);
%            beta0,alpha,beta0';
%            zeros(12),beta0,alpha]
       
Alpha_1 = kron(diag(ones(1,NW)),alpha)+kron(diag(ones(1,NW- 1),1),beta0')+kron(diag(ones(1,NW-1),-1),beta0);

% Beta_plus_1 = [t',beta2',zeros(12);
%                beta1,t',beta2';
%                zeros(12),beta1,t']

Beta_plus_1 = kron(diag(ones(1,NW)),t')+kron(diag(ones(1,NW- 1),1),beta2')+kron(diag(ones(1,NW-1),-1),beta1);

% Beta_minus_1 = [t',beta2',zeros(12);
%                 beta1,t',beta2';
%                 zeros(12),beta1,t'] 



% Alpha_2 = [alpha,t,zeros(12);
%     transpose(t),alpha,t;
%     zeros(12),transpose(t),alpha]
% 
% Beta_plus_2 = [beta0',beta1',zeros(12);
%         beta2',beta0',beta1';
%         zeros(12),beta2',beta0']
% 
% Beta_minus_2 = [beta0',beta1',zeros(12);
%         beta2',beta0',beta1';
%         zeros(12),beta2',beta0']

Alpha_2 = kron(diag(ones(1,NW)),alpha)+kron(diag(ones(1,NW- 1),1),t)+kron(diag(ones(1,NW-1),-1),t');


Beta_plus_2 = kron(diag(ones(1,NW)),beta0')+kron(diag(ones(1,NW- 1),1),beta1')+kron(diag(ones(1,NW-1),-1),beta2');


% alphau = t*diag(ones(1,NU-1),1)+t*diag(ones(1,NU-1),-1); 
% betaw = zeros(NU);
% betal = zeros(NU);

%Armchair Structure for GNR
if ((type == 1)&&(device == 1))
    Alpha = Alpha_1;
    Beta = Beta_plus_1;
elseif ((type == 2)&&(device == 1))

%Zigzag Structure for GNR
    Alpha = Alpha_2;
    Beta = Beta_plus_2;
elseif ((type == 1)&&(device == 2))
end


l1 = sqrt(3)*a;
l2 = 3*a;

if (type == 1)
    
    fprintf(file_out,'using gauge1 \n');
    
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
            Alpha(i,j) = Alpha(i,j)*exp(1i*q*(-Bz)*(y_i+y_j)/2*(x_i-x_j)/hbar*AtoM^2);
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
            Beta(i,j) = Beta(i,j)*exp(1i*q*(-Bz)*(y_i+y_j)/2*(x_i-x_j+l2)/hbar*AtoM^2);
        end
    end
    
else
    
    fprintf(file_out,'using gauge2 \n');
    
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
            Alpha(i,j) = Alpha(i,j)*exp(1i*q*(-Bz)*(x_i+x_j)/2*(y_i-y_j)/hbar*AtoM^2);
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
            Beta(i,j) = Beta(i,j)*exp(1i*q*(-Bz)*(x_i+x_j)/2*(y_i-y_j+l1)/hbar*AtoM^2);
        end
    end

end



% Hamiltonian: 
H = kron(diag(ones(1,NL)),Alpha)+kron(diag(ones(1,NL- 1),1),Beta)+kron(diag(ones(1,NL-1),-1),Beta');


% %Define Energy grid for calculation of Transmission
% E = linspace(-1,1,51);
% %Define the matrices for NEGF
% T = zeros(1,length(E));
% green1 = inv(E(1)*eye(NW*NU)-Alpha);
% green2 = inv(E(1)*eye(NW*NU)-Alpha);
% Es = zeros(NT);
% Esin = zeros(NT);
% fprintf(1,'Done\n');
% 
% %% NEGF Calculations
% fprintf(1,'Running NEGF...please wait\n')
% fprintf(1,'Energy range = %d eV to %d eV\n',min(E),max(E)); 
% fprintf(1,'%s\n',['[',blanks(length(E)),']']);
% fprintf(1,' ');
% 
% %Run for each energy
% for k = 1:length(E)
% %Calculate surface GFs self consistently
% %Calculate for beta'
%     err = 100;
%     while(err>errmax)
%         g1new = inv((E(k)+zplus)*eye(NW*NU)-Alpha-Beta'*green1*Beta); 
%         err = (sum(sum(abs(g1new-green1))))/(sum(sum(abs(g1new+green1)))); 
%         green1 = g1new;
%     end
%     sigma1  = Beta'*green1*Beta;
% %Calculate for beta
%     err = 100;
%     while(err>errmax)
%         g2new = inv((E(k)+zplus)*eye(NW*NU)-Alpha-Beta*green2*Beta'); 
%         err = (sum(sum(abs(g2new-green2))))/(sum(sum(abs(g2new+green2)))); 
%         green2 = g2new;
%     end
%     sigma2  = Beta*green2*Beta';
% %Calculate self energy matrices
%     E1 = kron(diag([1 zeros(1,NL-1)]),sigma1);
%     E2 = kron(diag([zeros(1,NL-1) 1]),sigma2);
% %Calculate broadening
%     G1 = 1i*(E1-E1'); 
%     G2 = 1i*(E2-E2');
% %Calculate G, Gn, A, T for coherent transport
% %     G = inv((E(k)+1i*etaplus)*eye(NT)-H-E1-E2);
% %     T(k) = real(trace(G1*G*G2*G'));
% %     Gn = f1*(G*G1*G')+f2*(G*G2*G');
% %     A = i*(G-G');
% %Calculate G, Gn, A, T self consistently including the phase breaking processes
%     err = 100;
%     while(err>errmax)
%         G = inv((E(k)+zplus)*eye(NT)-H-E1-E2-Es);
%         Esnew = D*G;
%         err = sum(sum(abs(Esnew-Es)))/sum(sum(abs(Esnew+Es))); 
%         Es = Esnew;
%     end
%     % % % % % % % %
% %      err = 100;
% %      while(err>errmax)
% %          Gn = f1*(G*G1*G')+f2*(G*G2*G')+(G*Esin*G');
% %          Esinnew = D*Gn;
% %          err = sum(sum(abs(Esinnew-Esin)))/sum(sum(abs(Esinnew+Esin))); Esin = Esinnew;
% %      end
% %      A=i*(G-G');
%     T(k) = real(trace(G1*G*G2*G'));
%     fprintf(1,'|');
% end
% 
% fprintf(1,'\n');
% fprintf(1,'Done\n');

%Calculate Subbands
ka = linspace(-pi,pi,5001);
eigESL = zeros(length(ka),NW*NU/3); 
eigEBL = zeros(length(ka),NW*NU*2/3); 
eigE = zeros(length(ka),NW*NU); 
fprintf(file_out,'Calculating Sub-bands...please wait\n');

dims = NW*NU;
MM = kron(diag(ones(1,NW)),M);

% C1 = [1,0,0,0;
%     0,1,0,0;
%     0,0,1,0;
%     0,0,0,1];
% C0 = [0,0,0,0;
%     0,0,0,0;
%     0,0,0,0;
%     0,0,0,0];
% CC = [C1,C0,C0,C0,C0,C0,C0,C0,C0,C0,C0,C0;
%       C0,C0,C0,C1,C0,C0,C0,C0,C0,C0,C0,C0;
%       C0,C0,C0,C0,C0,C0,C1,C0,C0,C0,C0,C0;
%       C0,C0,C0,C0,C0,C0,C0,C0,C0,C1,C0,C0;
%       C0,C1,C0,C0,C0,C0,C0,C0,C0,C0,C0,C0;
%       C0,C0,C1,C0,C0,C0,C0,C0,C0,C0,C0,C0;
%       C0,C0,C0,C0,C1,C0,C0,C0,C0,C0,C0,C0;
%       C0,C0,C0,C0,C0,C1,C0,C0,C0,C0,C0,C0;
%       C0,C0,C0,C0,C0,C0,C0,C1,C0,C0,C0,C0;
%       C0,C0,C0,C0,C0,C0,C0,C0,C1,C0,C0,C0;
%       C0,C0,C0,C0,C0,C0,C0,C0,C0,C0,C1,C0;
%       C0,C0,C0,C0,C0,C0,C0,C0,C0,C0,C0,C1];

CC = zeros(NW*NU);
  

for i = 1:NW*NU/3
    for j = 1:NW*NU
        if (j == i+8*floor((i-0.5)/4))
            CC(i,j) = 1;            
        end
    end    
end  

for i = NW*NU/3+1:NW*NU
    for j = 1:NW*NU
        if (j== i - 4*(NW-1-floor((i-NW*NU/3-0.5)/8)))
            CC(i,j) = 1;            
        end
    end
end

toc;

for k = 1:length(ka)
    
    HH = Alpha + Beta*exp(1i*ka(k))+ Beta'*exp(-1i*ka(k));
    
    if Basis==1
    
        HHM = MM*HH*MM';
        HHMC = CC*HHM*CC';
    
%     HHM_SL= zeros(dims/3);
%     HHM_BL= zeros(2*dims/3);

        HHM_SL= HHMC(1:dims/3,1:dims/3);
        HHM_BL= HHMC(1+dims/3:dims,1+dims/3:dims);
    
        [VSL,DSL] = eig(HHM_SL);
        [VBL,DBL] = eig(HHM_BL);
    
        eigESL(k,:) = diag(DSL);
        eigEBL(k,:) = diag(DBL);
    end
        
    if Basis==2
        [V,D] = eig(HH);
        eigE(k,:) = diag(D);
    end
end

fprintf(file_out,'Done\n');

% %% Plot Figures
% %Figure 1: Transmission and sub-band comparison figure(1);
% subplot(1,2,1);
% plot(T,E);
% xlabel('T(E)');
% ylabel('E[eV]');
% title('Energy v Transmission');
% xlim([0 max(T)+0.5]);
% 
% ylim([min(E) max(E)]);
% subplot(1,2,2);
% plot(ka/pi,eigE);
% xlabel('ka/\pi');
% ylabel('E[eV]');
% title('Sub-bands');
% ylim([min(E) max(E)]);

if Basis==1
    figure;

    %plot(ka/pi,eigE,'b');

    plot(ka/pi,eigESL,'b');
    hold on; 
    plot(ka/pi,eigEBL,'r');

    xlabel('ka/\pi');
    ylabel('E[eV]');
    title('Bands');
    
    xlim([-0.3,0.6]);
    ylim([-0.1,0.1]);
end

if Basis==2
    figure;
    plot(ka/pi,eigE,'b');
    xlabel('ka/\pi');
    ylabel('E[eV]');
    title('Mixed Bands');
    
    xlim([-0.3,0.6]);
    ylim([-0.1,0.1]);
    
end

%save(['bands_width100_100T_',num2str(dE*1000),'meV.mat'])

fprintf(kps_out,'%f\n',ka/pi);

if Basis==1
    fprintf(eigEBL_out,'%f\n',eigEBL);
    fprintf(eigESL_out,'%f\n',eigESL);
else
    fprintf(eigE_out,'%f\n',eigE);
end


toc;

%saveas(gcf,['bands_type=',num2str(type),'_Bz=',num2str(Bz),'.png'])