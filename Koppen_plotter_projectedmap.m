
%This script plots maps of the simplified Koppen class for different scenarios. 
%It uses the (simplified) Koppen classes input data with names like
%Koppendata_xxxx.mat (where xxxx is the scanario name, see below). 
%It is assumed that the intermediate data is in the same folder as the
%script. It is also assumed that the matlab function Koppen_formatter.m is
%in the same folder. 


scenchoi=1:7; 
  %1 REF, 2 CONTe, 3 CONTm, 4 CONTl, 5 S20e, 6 S20l,  7 S80l
nscen=length(scenchoi);

% scenario nmubers
 %1: REF (present-day control, i.e.. 2020-2035)
 %2: CONT_early (control 2050-2065; control means no SAI applied)
 %3: CONT_medium (control 2065-2080, i.e. efore SAI2080 starts)
 %4: CONT_late (control 2085-2100)
 %5: S20_early (SAI_2020 2050-2065 i.e. SAI starts in 2020)
 %6: S20_late  (SAI_2020 2085-2100 )
 %7: S80_late  (SAI_2080 2085-2100 i.e. SAI starts in 2080)
  
Nx=288;
Ny=192;

delx=360/(Nx-1);
dely=180/(Ny-1);

xvec=-180+delx/2:delx:180-delx;
yvec=-90+dely/2:dely:90-dely;

for k=scenchoi
    

    if scenchoi(k)==1
        load('Koppendata_REF.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,1);
        titlestring='Present-day';
    end
    
    if scenchoi(k)==2
        load('Koppendata_CONTe.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,1);
        titlestring='RCP8.5 2050-2065';
    end
    
    if scenchoi(k)==3
        load('Koppendata_CONTm.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,0);
        titlestring='RCP8.5 2065-2080'; 
    end
    
    if scenchoi(k)==4
        load('Koppendata_CONTl.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,1);
        titlestring='RCP8.5 2085-2100';
    end
    
    if scenchoi(k)==5
        load('Koppendata_S20e.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,1);
        titlestring='SAI2020 2050-2065';
    end
    
    if scenchoi(k)==6
        load('Koppendata_S20l.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,1);
        titlestring='SAI2020 2085-2100';
    end
    
    if scenchoi(k)==7
        load('Koppendata_S80l.mat')
        Koppen_aux=Koppen_formatter(Koppenraw,Nx,Ny,1);
        titlestring='SAI2080 2085-2100';
    end
    
    
 Koppen=Koppen_aux;   


% % % % % Koppen1=reshape(Koppenraw,Nx,Ny)';
% % % % % 
% % % % % 
% % % % % load('coastlines_0_360.mat');
% % % % % ff=find(coastlon>180);
% % % % % coastlon(ff)=coastlon(ff)-360;
% % % % % 
% % % % % Koppen2=Koppen1(Ny:-1:1,:);
% % % % % Koppen=Koppen2(:,[Nx/2+1:Nx,1:Nx/2]);
% % % % % 
% % % % % clear Koppen1; clear Koppen2
% % % % % %Re-organise some climate zones
% % % % % %NEW name: -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13
% % % % % %OLD name:  0   1 2 7 8 9 6 4 5 3 10 11 12 13
% % % % % 
% % % % % ff=find(Koppen==-1); Koppen(ff)= 0.1;
% % % % % 
% % % % % % ff=find(Koppen== 7); Koppen(ff)= 3.1;
% % % % % % ff=find(Koppen== 8); Koppen(ff)= 4.1;
% % % % % % ff=find(Koppen== 9); Koppen(ff)= 5.1;
% % % % % % ff=find(Koppen== 4); Koppen(ff)= 7.1;
% % % % % % ff=find(Koppen== 5); Koppen(ff)= 8.1;
% % % % % % ff=find(Koppen== 3); Koppen(ff)= 9.1;
% % % % % Koppen=round(Koppen);
    
    %manually make some big lakes white (zero)
    Koppen(47:57,185:186)=0; %caspian
    Koppen(56:57,187)=0;
    Koppen(48:49,184)=0;
    
    Koppen(48:49,192:193)=0; %Aral 
    
    Koppen(97:99,171:172)=0; %Victoria
    
    Koppen(46:46,73:77)=0; %great lakes
    Koppen(45:45,75:76)=0;
    Koppen(47:47,76:76)=0;
    Koppen(48:52,75:75)=0;
    Koppen(48:50,78:79)=0;

%colors based on fig 12 a in: https://personal.sron.nl/~pault/ 
%see also

colmat=[...
 255 255 255      %sea
 160 24  19     %tropical rainforest
 233 76 31     %subtropical forest
 253 154 68    %tropical grassland and savannah
 249 213 118     %deserts and xeric shrublands
 240 230 178     %mediterranean shrublands
 18 90 86     %warm temperate and coastal forests
 0 118 123     %temperate wet forests (small group)
 66 167 198    %temperate grassland and savannah
 96 188 233     %various forests (fairly temperate and continental)
 157 204 239     %montane and arctic grass and shrubland
 198 219 237     %boreal coniferous forest taiga
 236 234 218     %boreal and mountain tundra
 210 210 210     %ice cap
    ]/255;
    

% colmat=[...
%     1    1   1
%     0    0   0
%     0    0   1
%     0.5 0.5  1
%     0    0.5   0.5
%     0   0.5  0
%     0.3  0.8  0
%     0   0.6   0.3
%     1  0.6  0
%     0.9 0.8 0
%     0.5  0.2  0
%     0.5 0  0.5
%     0.8  0  0.8
%     0.4 0.6 0.8
%     0.6 0.8 1];

h=figure; hold on
worldmap('World')
load coastlines
surfm(yvec,xvec,Koppen(end:-1:1,:))
plotm(coastlat,coastlon,'Color','0.5, 0.5, 0.5]','linewidth',0.5)
set(gca,'FontSize',15)
% % % scatter(coastlon,coastlat,'k.') %just doing this to enfore that north is up 
% % % imagesc(xvec,-yvec,Koppen)
% % % scatter(coastlon,coastlat,0.1,'k.') %now actually adding the coast line
colormap(colmat)
hcb = colorbar;
ticks=(0:16)*0.93+0.5;
hcb.Ticks=ticks;
hcb.TickLength = 0;
%Koppenlabels={'sea','trop forest','subtrop forest','savannah','deserts','mediterr','warm temp forest','wet temp forest','temp grass','contin forest','subarctic shrub','boreal forest','boreal grass','ice cap','   '};
Koppenlabels={'sea','trop for','strop for','savannah','desert','medit','temp for','temp for,m','temp grass','cont for','arc shrub','boreal for','bor grass','ice cap','   ','      '};
%Koppenlabels={'sea','tropical forest','subtropical forest','savannah','desert','mediterranean','temperate forest','temperate forest, moist','temperate grassland','continental forest','arctic shrubland','boreal forest','boreal grassland','ice cap','   ','      '};

hcb.TickLabels=Koppenlabels;
%h.Position(1) = h.Position(1)-0.1;
%axis([-180,180,-90,90])
title(titlestring)
set(gca,'FontSize',15)
set(findall(h,'Tag','PLabel'),'visible','off')
set(findall(h,'Tag','MLabel'),'visible','off')

%colmat  (-1 = sea; 13 = antarctica)

end %looop over all scenarios to be plotted