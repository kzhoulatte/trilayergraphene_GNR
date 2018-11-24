
clear all; clc;

NW = input('Width of ribbon: ');
kps = fopen('TNR_kps.txt','r');
k = fscanf(kps,'%f');

Basis = input('Basis used(1 for odd&even;2 for mixed): ');
[ksize,m] = size(k);
klim=0.2;

if Basis==1
    eigEBL_out = fopen('TNR_eigEBL.txt','r');
    eigESL_out = fopen('TNR_eigESL.txt','r');
    eigESL_r = fscanf(eigESL_out,'%f');
    eigEBL_r = fscanf(eigEBL_out,'%f');
    
    eigESL = zeros(ksize,NW*4);
    eigEBL = zeros(ksize,NW*8);
    
    for i=1:NW*4
        for j = 1:ksize
            eigESL(j,i) = eigESL_r((i-1)*ksize+j);
        end
    end
    
    for i=1:NW*8
        for j = 1:ksize
            eigEBL(j,i) = eigEBL_r((i-1)*ksize+j);
        end
    end
    
    
else
    eigE_out = fopen('TNR_eigE.txt','r');
    eigE_r = fscanf(eigE_out,'%f');
    eigE = zeros(ksize,NW*12);
    
     for i=1:NW*12
        for j = 1:ksize
            eigE(j,i) = eigE_r((i-1)*ksize+j);
        end
    end
    
end

if Basis==1
    figure;

    %plot(ka/pi,eigE,'b');

    plot(k,eigESL,'b');
    hold on; 
    plot(k,eigEBL,'r');

    xlabel('ka/\pi');
    ylabel('E[eV]');
    title('Bands');
    
    xlim([-klim,klim]);
    ylim([-0.1,0.1]);
end

if Basis==2
    figure;
    
    plot(k,eigE,'b');
    xlabel('ka/\pi');
    ylabel('E[eV]');
    title('Mixed Bands');
    
    xlim([-klim,klim]);
    ylim([-0.1,0.1]);
    
end