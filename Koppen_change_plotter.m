
%This script explores how land grid cells change their simplified Koppen
%class between different moments of different scenarios. 
%It uses the (simplified) Koppen classes input data with names like
%Koppendata_xxxx.mat (where xxxx is the scanario name, see below). 
%It is assumed that the intermediate data is in the same folder as the
%script. It is also assumed that the matlab function Koppen_formatter.m is
%in the same folder. 


%Chose scenarios 

scenchoi=[3 7 7]; %which scenarios to compare 
% picking 3 different ones: 3-way comparison (e.g. [1 3 7])
% se the last two equal for a 2-way comparison (e.g. [1 3 3])
% scenario nmubers
 %1: REF (present-day control, i.e.. 2020-2035)
 %2: CONT_early (control 2050-2065; control means no SAI applied)
 %3: CONT_medium (control 2065-2080, i.e. efore SAI2080 starts)
 %4: CONT_late (control 2085-2100)
 %5: S20_early (SAI_2020 2050-2065 i.e. SAI starts in 2020)
 %6: S20_late  (SAI_2020 2085-2100 )
 %7: S80_late  (SAI_2080 2085-2100 i.e. SAI starts in 2080)
 
% The script first plots some maps of the changes in Koppen zones (2- or 3-way comparison) 
% next, it makes a matrix saying how much land area changed between which
% Köppen zones.
% Note that change matrix plots are always done between the first and second
% scenario chosen; if you pick [1 3 7], the matrix compared 1 (REF) and 3 (Control medium)
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No changes below this line unless you know what you're doing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=288;
Ny=192;


for k=1:3

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
    
    %manually make some big lakes white (zero)
    Koppen_aux(47:57,185:186)=0; %caspian
    Koppen_aux(56:57,187)=0;
    Koppen_aux(48:49,184)=0;
    
    Koppen_aux(48:49,192:193)=0; %Aral 
    
    Koppen_aux(97:99,171:172)=0; %Victoria
    
    Koppen_aux(46:46,73:77)=0; %great lakes
    Koppen_aux(45:45,75:76)=0;
    Koppen_aux(47:47,76:76)=0;
    Koppen_aux(48:52,75:75)=0;
    Koppen_aux(48:50,78:79)=0;

    if k==1
        Koppen1=Koppen_aux;
        tstring1=titlestring; 
    elseif k==2
        Koppen2=Koppen_aux; 
        tstring2=titlestring; 
    elseif k==3
        Koppen3=Koppen_aux; 
        tstring3=titlestring;
    end
end 

delx=360/(Nx-1);
dely=180/(Ny-1);

xvec=-180+delx/2:delx:180-delx;
yvec=-90:dely:90;


load('coastlines_0_360.mat');
ff=find(coastlon>180);
coastlon(ff)=coastlon(ff)-360;


%compute Koppen difference 

Koppendiff=Koppen1*0; 

ff=find(Koppen1>0); %land, general -- unchanged
Koppendiff(ff)=1;

ff=find(abs(Koppen1-Koppen2)>0 & abs(Koppen1-Koppen3)==0 ); %restoration
Koppendiff(ff)=2; 

ff=find(abs(Koppen1-Koppen2)>0 & abs(Koppen2-Koppen3)==0 ); %change by first switch
Koppendiff(ff)=3; 
    
ff=find(abs(Koppen1-Koppen2)==0 & abs(Koppen2-Koppen3)>0 ); %change by second switch
Koppendiff(ff)=4; 

ff=find(abs(Koppen1-Koppen2)>0 & abs(Koppen2-Koppen3)>0 & abs(Koppen1-Koppen3)>0 ); %change and change
Koppendiff(ff)=5; 

% areas 
cosvec=cos(yvec*pi/180)';

area_same_ind=Koppendiff;  %normalised area that is NOT changed 
ff=find(area_same_ind>1.5);  %eliminate changed and restored
area_same_ind(ff)=0; 
area_same_ind=sum(area_same_ind,2); %sum over longitudes;
area_same_ind=sum(area_same_ind.*cosvec);

area_rest_ind=Koppendiff;  %normalised area that is restored
ff=find(abs(area_rest_ind-2)>0.1);  %eliminate unchanged and permanently changed
area_rest_ind(ff)=0; 
area_rest_ind=area_rest_ind./(area_rest_ind+0.00000000001);
area_rest_ind=sum(area_rest_ind,2); %sum over longitudes;
area_rest_ind=sum(area_rest_ind.*cosvec);

area_chan_ind=Koppendiff;  %normalised area that IS changed
ff=find(area_chan_ind<2.5);  %eliminate non-changed and restores area
area_chan_ind(ff)=0; 
area_chan_ind=area_chan_ind./(area_chan_ind+0.00000000001);
area_chan_ind=sum(area_chan_ind,2); %sum over longitudes;
area_chan_ind=sum(area_chan_ind.*cosvec);

area_land=area_chan_ind + area_rest_ind+ area_same_ind;
area_frac_same= area_same_ind / (area_land) * 100; %unchanged area in % 
area_frac_rest= area_rest_ind / (area_land) * 100; %restored area in % 
area_frac_chan= area_chan_ind / (area_land) * 100; %changed area in % 


%colors based on fig 12 a in: https://personal.sron.nl/~pault/ 
%see also

colmat3=[...   %for 3 scenarios
 255 255 255      %sea
 198 219 237       %no change
 18 90 86         %restored
 249 213 118      % change 1st shift
 253 154 68       % change 2nd shift
 233 76 31        %two changes
    ]/255;
    
colmat2=[...   %for 2 scenarios
 255 255 255      %sea
 198 219 237        %no change
 253 154 68      % change 
    ]/255;
    

h=figure; hold on
worldmap('World')
load coastlines
surfm(yvec,xvec,Koppendiff(end:-1:1,:))
plotm(coastlat,coastlon,'Color','0.5, 0.5, 0.5]','linewidth',0.5)
hcb = colorbar;
if scenchoi(2)==scenchoi(3) %only one comparison
 colormap(colmat2)
 ticks=(0:2)+0.5;
 caxis([0 3])
 hcb.Ticks=ticks;
 hcb.TickLength = 0;
 Koppenlabels={'sea','no chan','chan'};
 hcb.TickLabels=Koppenlabels;
 title([tstring1 '  vs  ' tstring2 ' ;  land area no change: ' num2str(round(area_frac_same))  '%, change: ' num2str(round(area_frac_chan)) '%'])
else
 colormap(colmat3)
 ticks=(0:5)+0.5;
 caxis([0 6])
 hcb.Ticks=ticks;
 hcb.TickLength = 0;
 Koppenlabels={'sea','no chan','restored', 'chan 1','chan 2','chan 1+2'};
 hcb.TickLabels=Koppenlabels;
 title([tstring1 '  vs  ' tstring2 '  vs  ' tstring3 ' ;  land area no change: ' num2str(round(area_frac_same)) '%, restored: ' num2str(round(area_frac_rest)) '%, change: ' num2str(round(area_frac_chan)) '%'])
end
%axis([-180,180,-90,90])
set(gca,'FontSize',15)
set(findall(h,'Tag','PLabel'),'visible','off')
set(findall(h,'Tag','MLabel'),'visible','off')



% now we also make tables wit changes between Koppen1 and Koppen2
% entry i,j is: how much land was class i in Koppen1 and is now class j in Koppen2

Kchanmat=zeros(13,13); % first index: orig

for iy=1:Ny
for ix=1:Nx
    
    roundK1=round(Koppen1(iy,ix));
    roundK2=round(Koppen2(iy,ix));
    
    if roundK1>0.5 && roundK2>0.5 
        Kchanmat(roundK1,roundK2)=Kchanmat(roundK1,roundK2)+cosvec(iy);
    end
end
end

Kchanmat=Kchanmat/area_land;


landlabels={'trop F','subtrop F','savannah','desert','medit','temp F','temp wet F','temp G','cont F','arctic S','boreal F','boreal G','ice'};
nl=length(landlabels);

landlabels1=landlabels; 
landlabels2=landlabels; 
for k=1:nl,
    landlabels1{k} = [landlabels1{k} ' (' num2str(round(sum(Kchanmat(k,:))*1000)/10) ')'];
    landlabels2{k} = [landlabels2{k} ' (' num2str(round(sum(Kchanmat(:,k))*1000)/10) ')'];
end

% %from Iceland paper
%   grouptickarray={['vagabond A (' num2str(SD(1)) ')'],['vagabond C (' num2str(SD(2)) ')'],['pauper A (' num2str(SD(3)) ')']...
%       ,['pauper C (' num2str(SD(4)) ')'],['relative (' num2str(SD(5)) ')'],['foster C (' num2str(SD(6)) ')']...
%       ,['servant (' num2str(SD(7)) ')'],['cottar  (' num2str(SD(8)) ')'],['C of cottar (' num2str(SD(9)) ')'],['farmer (' num2str(SD(10)) ')']...
%       ,['C of farmer (' num2str(SD(11)) ')'],['elite (' num2str(SD(12)) ')'],['C of elite (' num2str(SD(13)) ')'],['widowed (' num2str(SD(14)) ')']};

figure; 
set(gca,'fontsize',14)
imagesc(Kchanmat*100) %plot as %... 

% colourmat
vec=round([1 sqrt(0.002:0.002:1)*256]);
clear hot
hot2=hot;
hot3=hot2(vec,:);
colormap(hot3); 
colorbar
%label stuff
title(['land area changes (%) from ' tstring1 ' to ' tstring2])
xticks(1:nl)
xticklabels(landlabels2)
xlabel({[tstring2 ],[ ' Note: F= forest, S=shrub, G=grass']})
xline([2 5 8 11]+0.5,'w--')
yticks(1:nl)
yticklabels(landlabels1)
ylabel(tstring1)
yline([2 5 8 11]+0.5,'w--')
set(gca,'fontsize',14)
caxis([0 24])


