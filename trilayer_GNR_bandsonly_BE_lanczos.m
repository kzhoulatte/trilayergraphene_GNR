% The band structure and NEGF for monolayer graphene Nanoribbons
% The major difference between graphene with infinite boundry-peoridic
% condictions: the unit cell is rectangular.

%% Device Setup

load('bands_width100_100T_100meV.mat')
tic;

Basis = input('Basis used(1 for odd&even;2 for mixed): ');
ka = linspace(-pi,pi,51);
ktop = 30;
eigESL = zeros(length(ka),ktop/3); 
eigEBL = zeros(length(ka),ktop*2/3); 
eigE = zeros(length(ka),ktop); 

for k = 1:length(ka)
    
    HH = Alpha + Beta*exp(1i*ka(k))+ Beta'*exp(-1i*ka(k));
    
    if Basis==1
    
        HHM = MM*HH*MM';
        HHMC = CC*HHM*CC';
    
%     HHM_SL= zeros(dims/3);
%     HHM_BL= zeros(2*dims/3);

        HHM_SL= HHMC(1:dims/3,1:dims/3);
        HHM_BL= HHMC(1+dims/3:dims,1+dims/3:dims);
        
        
        %TSL = myLanczos(HHM_SL,ktop/3);
        [VSL,DSL] = eigs(HHM_SL,ktop/3,'SM');
        %TBL = myLanczos(HHM_BL,ktop*2/3);
        [VBL,DBL] = eigs(HHM_BL,ktop*2/3,'SM');
    
        eigESL(k,:) = diag(DSL);
        eigEBL(k,:) = diag(DBL);
    end
        
    if Basis==2
        
        %T = myLanczos(HH,ktop);
        
        [V,D] = eigs(HH,ktop,'SM');
        eigE(k,:) = diag(D);
    end
end

fprintf(1,'Done\n');

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

% if Basis==1
%     figure;
% 
%     %plot(ka/pi,eigE,'b');
% 
%     plot(ka/pi,eigESL,'b');
%     hold on; 
%     plot(ka/pi,eigEBL,'r');
% 
%     xlabel('ka/\pi');
%     ylabel('E[eV]');
%     title('Bands');
%     
%     xlim([-0.3,0.6]);
%     ylim([-0.1,0.1]);
% end
% 
% if Basis==2
%     figure;
%     plot(ka/pi,eigE,'b');
%     xlabel('ka/\pi');
%     ylabel('E[eV]');
%     title('Mixed Bands');
%     
% %     xlim([-0.3,0.6]);
% %     ylim([-0.1,0.1]);
%     
% end

% Scatter plots: 

if Basis==1
    figure;

    %plot(ka/pi,eigE,'b');
    for ii=1:length(ka)
        kpsSL = ka(ii)/pi*ones(ktop/3,1);
        scatter(kpsSL,eigESL(ii,:),'filled');
        hold on; 
        kpsBL = ka(ii)/pi*ones(ktop*2/3,1);
        scatter(kpsBL,eigEBL(ii,:),'filled');
    end
    xlabel('ka/\pi');
    ylabel('E[eV]');
    title('Bands');
    
    xlim([-0.3,0.6]);
    ylim([-0.1,0.1]);
end

if Basis==2
    figure;
    
    for ii=1:length(ka)
        kps = ka(ii)/pi*ones(ktop,1);
        scatter(kpsSL,eigE(ii,:),'filled');
    end
    
    xlabel('ka/\pi');
    ylabel('E[eV]');
    title('Mixed Bands');
    
%     xlim([-0.3,0.6]);
%     ylim([-0.1,0.1]);
    
end


toc;
%save(['bands_width100_100T_',num2str(dE*1000),'meV.mat'])

%toc;

%saveas(gcf,['bands_type=',num2str(type),'_Bz=',num2str(Bz),'.png'])