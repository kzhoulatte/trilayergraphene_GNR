% This is the script I need to plot bands change vs (B,E). The band
% positions E can extracted from output files of calculations. 

clear all;
clc;

BE_file = fopen('./BEbands.txt','r');
BE_result = fscanf(BE_file,'%f');

B = 'B';
E = 'E';


NB = 7;
NE = 4;
Nbands = 8;

BE0 = reshape(BE_result,[NB*NE*2,Nbands/2]);
BE1 = zeros(NB*NE,Nbands);

for i = 1:NB*NE
    BE1(i,:) = reshape(BE0(i*2-1:i*2,:),[1,Nbands]);
end

%BE = zeros(NB,NE,Nbands);
BE = reshape(BE1,[NB,NE,Nbands]);
disp('B = 2:8,E = 0:10:30. ');

fixed = input('Whose effect you want see: ');


if fixed == 'E'
    Bvalue = input('Input B field to be fixed on(2:1:8): ');
    
    figure;
    for i = 1:8
        if i==3 || i==4
            %scatter(0:10:30,BE(Bvalue-1,:,i),'r','.');
            plot(0:10:30,BE(Bvalue-1,:,i),'r');
            hold on;
        else
            %scatter(0:10:30,BE(Bvalue-1,:,i),'b','.');
            plot(0:10:30,BE(Bvalue-1,:,i),'b');
            hold on;
        end
    end
    
else
    Evalue = input('Input E field to be fixed on(0,10,20,30): ');
    
    figure;
    for i = 1:8
        if i==3 || i==4
            %scatter(2:1:8,BE(:,Evalue/10+1,i),'r','.');
            plot(2:1:8,BE(:,Evalue/10+1,i),'r');
            hold on;
        else
            %scatter(2:1:8,BE(:,Evalue/10+1,i),'b','.');
            plot(2:1:8,BE(:,Evalue/10+1,i),'b');
            hold on;
        end
    end
    
end


% [X,Y,Z] = ndgrid(1:NB,1:NE,1:Nbands);
% pointsize =30;
% %scatter3(X(:),Y(:),Z(:),pointsize,BE(:));
% scatter3(X(:),Y(:),BE(:));
% hold on; 
% 
% BE_single = BE(:,:,3:4);
% [X,Y,Z] = ndgrid(1:NB,1:NE,1:Nbands/4);
% scatter3(X(:),Y(:),BE_single(:),'r');
% 
% view(0,0);
% 
% view(90,0);



