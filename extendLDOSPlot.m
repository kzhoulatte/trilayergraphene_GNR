
% clear all; clc;
% warning('off','all');
% 
% fprintf('for GNR\n');
% device = 1;
% 
% fprintf(1,'Enter the configuration of device you want to simulate Armchair/Zigzag\n')
% fprintf(1,'1 for Armchair\n');
% fprintf(1,'2 for Zigzag\n');
% type = input('Valid inputs are 1 or 2 : ');
% 
% if (~((type == 1)||(type == 2)))
%     type = 1;
% end
% 
% NW = input('Enter number of Unit Cells in the width : ');
% if isempty(NW)
%     NW = 3;
% end

%Points in GNR length
NL = 1;
%Points in GNR unit cell
NU = 12;
NL_extended=input('Enter number of Unit Cells to extend in L: ');

a = 1.42;
c = 6.71;

A1 = [0,0,c];
B1 = [-1/2,sqrt(3)/2,c/a]*a;
A1_p = [3/2,sqrt(3)/2,c/a]*a;
B1_p = [1,0,c/a]*a;
A2 = [-1/2,sqrt(3)/2,0]*a;
B2 = [-1,0,0]*a;
A2_p = [1,0,0]*a;
B2_p = [1/2,sqrt(3)/2,0]*a;
A3 = [0,0,-c];
B3 = [-1/2,sqrt(3)/2,-c/a]*a;
A3_p = [3/2,sqrt(3)/2,-c/a]*a;
B3_p = [1,0,-c/a]*a;

Xlist = [0,-1/2,3/2,1,-1/2,-1,1,1/2,0,-1/2,3/2,1]*a;
Ylist = [0,sqrt(3)/2,sqrt(3)/2,0,sqrt(3)/2,0,0,sqrt(3)/2,0,sqrt(3)/2,sqrt(3)/2,0]*a;
Zlist = [c,c,c,c,0,0,0,0,-c,-c,-c,-c];

l1 = sqrt(3)*a;
l2 = 3*a;


%colormap hot
DOS_k_norm = DOS_k/max(DOS_k);
size = 5;
shrink=1;
shrinkA = 1;

map = [
%     1,0.9,0.9
    1,0.8,0.8
    1,0.7,0.7
    1,0.6,0.6
    1,0.5,0.5
    1,0.4,0.4
    1,0.3,0.3
    1,0.2,0.2
    1,0.1,0.1
    1,0,0
]

% if type ==1
%     set(gcf,'units','points','position',[-5 -5 10 10+l1*NW]);
% else 
%     set(gcf,'units','points','position',[-5 -5  10+l2*NW 10])
% end

figure;
if (type ==1)
    for k=1:NL_extended
        for i=1:NW
            scatter3(Xlist(1:4)+[l2*k,l2*k,l2*k,l2*k],Ylist(1:4)+[l1*i,l1*i,l1*i,l1*i],Zlist(1:4),10,DOS_k_norm(1+(i-1)*12:4+(i-1)*12),'filled');
            hold on;
            scatter3(Xlist(5:8)+[l2*k,l2*k,l2*k,l2*k],Ylist(5:8)+[l1*i,l1*i,l1*i,l1*i],Zlist(5:8),30,DOS_k_norm(5+(i-1)*12:8+(i-1)*12),'filled');
            hold on;
            scatter3(Xlist(9:12)+[l2*k,l2*k,l2*k,l2*k],Ylist(9:12)+[l1*i,l1*i,l1*i,l1*i],Zlist(9:12),50,DOS_k_norm(9+(i-1)*12:12+(i-1)*12),'filled');
            hold on;
        end
    end
    
else
    for k=1:NL_extended
      for i=1:NW
        scatter3(Xlist(1:4)+[l2*i,l2*i,l2*i,l2*i],Ylist(1:4)+[l1*k,l1*k,l1*k,l1*k],Zlist(1:4),10,DOS_k_norm(1+(i-1)*12:4+(i-1)*12),'filled');
        hold on;
        scatter3(Xlist(5:8)+[l2*i,l2*i,l2*i,l2*i],Ylist(5:8)+[l1*k,l1*k,l1*k,l1*k],Zlist(5:8),30,DOS_k_norm(5+(i-1)*12:8+(i-1)*12),'filled');
        hold on;
        scatter3(Xlist(9:12)+[l2*i,l2*i,l2*i,l2*i],Ylist(9:12)+[l1*k,l1*k,l1*k,l1*k],Zlist(9:12),50,DOS_k_norm(9+(i-1)*12:12+(i-1)*12),'filled');
        hold on;
      end
    end
end

if type==1
    xmin = -4;
    xmax = 4+l1*NW;;
    ymin = -4;
    ymax = 4+l1*NW;
    
    xlim([xmin,xmax*shrinkA]);
    ylim([ymin,ymax]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3)*shrink a(4)]);
    
else
    xmin = -5;
    xmax = 5+l2*NW;
    
    ymin = -5;
    ymax = 5+l2*NW;
    
    xlim([xmin,xmax]);
    ylim([ymin,ymax*shrinkA]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3) a(4)*shrink]);
    
end 
    

grid off; 
view(0,90);
colormap(map);
title('All three layers');

figure;
if (type ==1)
    for k=1:NL_extended
    for i=1:NW
        scatter3(Xlist(1:4)+[l2*k,l2*k,l2*k,l2*k],Ylist(1:4)+[l1*i,l1*i,l1*i,l1*i],Zlist(1:4),size,DOS_k_norm(1+(i-1)*12:4+(i-1)*12),'filled');
        hold on;
    end
    end
    
else
    for k=1:NL_extended
    for i=1:NW
        scatter3(Xlist(1:4)+[l2*i,l2*i,l2*i,l2*i],Ylist(1:4)+[l1*k,l1*k,l1*k,l1*k],Zlist(1:4),size,DOS_k_norm(1+(i-1)*12:4+(i-1)*12),'filled');
        hold on;
    end
    end
end

if type==1
    xmin = -4;
    xmax = 4+l1*NW;;
    ymin = -4;
    ymax = 4+l1*NW;
    
    xlim([xmin,xmax*shrinkA]);
    ylim([ymin,ymax]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3)*shrink a(4)]);
    
else
    xmin = -5;
    xmax = 5+l2*NW;
    
    ymin = -5;
    ymax = 5+l2*NW;
    
    xlim([xmin,xmax]);
    ylim([ymin,ymax*shrinkA]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3) a(4)*shrink]);
    
end 

grid off; 
view(0,90);
colormap(map);
title('Top layer');

figure;
if (type ==1)
     for k=1:NL_extended
    for i=1:NW
        
        scatter3(Xlist(5:8)+[l2*k,l2*k,l2*k,l2*k],Ylist(5:8)+[l1*i,l1*i,l1*i,l1*i],Zlist(5:8),size,DOS_k_norm(5+(i-1)*12:8+(i-1)*12),'filled');
        hold on;
        
    end
     end
    
else
     for k=1:NL_extended
    for i=1:NW
        
        scatter3(Xlist(5:8)+[l2*i,l2*i,l2*i,l2*i],Ylist(5:8)+[l1*k,l1*k,l1*k,l1*k],Zlist(5:8),size,DOS_k_norm(5+(i-1)*12:8+(i-1)*12),'filled');
        hold on;
       
    end
     end
end


if type==1
    xmin = -4;
    xmax = 4+l1*NW;;
    ymin = -4;
    ymax = 4+l1*NW;
    
    xlim([xmin,xmax*shrinkA]);
    ylim([ymin,ymax]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3)*shrink a(4)]);
    
else
    xmin = -5;
    xmax = 5+l2*NW;
    
    ymin = -5;
    ymax = 5+l2*NW;
    
    xlim([xmin,xmax]);
    ylim([ymin,ymax*shrinkA]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3) a(4)*shrink]);
    
end 

grid off; 
view(0,90);
colormap(map);
title('Middle layer');

figure;
if (type ==1)
     for k=1:NL_extended
    for i=1:NW
        
        scatter3(Xlist(9:12)+[l2*k,l2*k,l2*k,l2*k],Ylist(9:12)+[l1*i,l1*i,l1*i,l1*i],Zlist(9:12),size,DOS_k_norm(9+(i-1)*12:12+(i-1)*12),'filled');
        hold on;
    end
     end
    
else
     for k=1:NL_extended
    for i=1:NW
        
        scatter3(Xlist(9:12)+[l2*i,l2*i,l2*i,l2*i],Ylist(9:12)+[l1*k,l1*k,l1*k,l1*k],Zlist(9:12),size,DOS_k_norm(9+(i-1)*12:12+(i-1)*12),'filled');
        hold on;
    end
     end
end

if type==1
    xmin = -4;
    xmax = 4+l1*NW;;
    ymin = -4;
    ymax = 4+l1*NW;
    
    xlim([xmin,xmax*shrinkA]);
    ylim([ymin,ymax]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3)*shrink a(4)]);
    
else
    xmin = -5;
    xmax = 5+l2*NW;
    
    ymin = -5;
    ymax = 5+l2*NW;
    
    xlim([xmin,xmax]);
    ylim([ymin,ymax*shrinkA]);
    
    a = get(gca,'position');
    set(gca, 'position',[a(1) a(2) a(3) a(4)*shrink]);
    
end 

grid off; 
view(0,90);
colormap(map);
title('Bottom layer');





