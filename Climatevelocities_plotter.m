
%This script plots climate velocities (for temperauture-related climate
%velocities (T) or precipitation-related (prect) velocities) or ratios
%thereof, depending on the "choice" parameter below. 

%it is assumed that the velocity intermediate data (v_X_Y.nc where X is "T"
%ot "prect" and Y is the scenario (con, sai2020, sai2080) is in the same
%folder. Also, thefile coastlines.mat should be in the same folder. 

choice=13; 



if choice==1; filename='v_T_con.nc';     titlestring = 'Velocity of 2m Temp for RCP8.5, 2020-2100 (in km/year)'; end
if choice==2; filename='v_T_sai2020.nc'; titlestring = 'Velocity of 2m Temp for SAI2020, 2020-2100 (in km/year)'; end
if choice==3; filename='v_T_sai2080.nc'; titlestring = 'Velocity of 2m Temp for SAI2080, 2080-2090 (in km/year)'; end
if choice==4; filename='v_prect_con.nc';     titlestring = 'Velocity of Precip for RCP8.5, 2020-2100 (in km/year)'; end
if choice==5; filename='v_prect_sai2020.nc'; titlestring = 'Velocity of Precip for SAI2020, 2020-2100 (in km/year)'; end
if choice==6; filename='v_prect_sai2080.nc'; titlestring = 'Velocity of Precip for SAI2080, 2080-2090 (in km/year)'; end

if choice==11; filename1='v_T_sai2020.nc';   filename2='v_T_con.nc';     titlestring = 'Ratio |velocity of 2m Temp|  SAI2020, 2020-2100 vs RCP8.5, 2020-2100'; end
if choice==12; filename1='v_T_sai2080.nc';   filename2='v_T_con.nc';     titlestring = 'Ratio |velocity of 2m Temp|  SAI2080, 2080-2090 vs RCP8.5, 2020-2100'; end
if choice==13; filename1='v_T_sai2080.nc';   filename2='v_T_sai2020.nc'; titlestring = 'Ratio |velocity of 2m Temp|  SAI2080, 2080-2090 vs SAI2020, 2020-2100'; end


% work on the variables to make them less ugly 
if choice<7
vel=ncread(filename,'__xarray_dataarray_variable__');
end

if choice>10
    vel1=ncread(filename1,'__xarray_dataarray_variable__');
    vel2=ncread(filename2,'__xarray_dataarray_variable__');
    vel=vel1./vel2; clear vel1, clear vel2
end 


Nx=288;
Ny=192;

delx=360/(Nx);
dely=180/(Ny-1);

velaux=[vel(Nx/2+1:end,:); vel(1:Nx/2,:)];
vel=velaux; clear velaux

xvec=-180:delx:180-delx;
yvec=-90:dely:90;



load('coastlines_0_360.mat');
ff=find(coastlon>180);
coastlon(ff)=coastlon(ff)-360;



% % % %colors based on fig 12 a in: https://personal.sron.nl/~pault/ 
% % % %see also
% % % 
% % % colmat=[...
% % %  255 255 255      %sea
% % %  160 24  19     %tropical rainforest
% % %  233 76 31     %subtropical forest
% % %  253 154 68    %tropical grassland and savannah
% % %  249 213 118     %deserts and xeric shrublands
% % %  240 230 178     %mediterranean shrublands
% % %  18 90 86     %warm temperate and coastal forests
% % %  0 118 123     %temperate wet forests (small group)
% % %  66 167 198    %temperate grassland and savannah
% % %  96 188 233     %various forests (fairly temperate and continental)
% % %  157 204 239     %montane and arctic grass and shrubland
% % %  198 219 237     %boreal coniferous forest taiga
% % %  236 234 218     %boreal and mountain tundra
% % %  210 210 210     %ice cap
% % %     ]/255;
% % %     

%choices for colour bar... 
nraw=7; %nr of entries in clomat*_aux below
if choice<7
maxpot=9; %plot from 2^-maxpot to 2^maxpot. 
delcol=(nraw-1)/(maxpot-1);
end
if choice>10
maxpot=7; %plot from 2^-maxpot to 2^maxpot. 
delcol=(nraw-1)/(2*maxpot-1);
end

%colmatpos_aux=[1 1 1; 1 0 0.6; 0.6 0.2 0; 1 0.9 0]; 
colmatpos_aux=[255 255 255 ;236 234 218  ;240 230 178 ;249 213 118 ;253 154 68 ;233 76 31 ;160 24  19;    ]/255; 
colmatpos=[interp1(1:nraw,colmatpos_aux(:,1),1:delcol:nraw); interp1(1:nraw,colmatpos_aux(:,2),1:delcol:nraw); interp1(1:nraw,colmatpos_aux(:,3),1:delcol:nraw)]';

colmatneg_aux= [255 255 255 ; 198 219 237 ; 157 204 239 ; 96 188 233  ; 66 167 198 ; 0 118 123 ; 18 90 86 ;   ]/255; 
colmatneg=[interp1(1:nraw,colmatneg_aux(:,1),1:delcol:nraw); interp1(1:nraw,colmatneg_aux(:,2),1:delcol:nraw); interp1(1:nraw,colmatneg_aux(:,3),1:delcol:nraw)]';

ticks=[-10:2:-2 2:2:10] ;

if choice < 7
velaux=vel; 
ff=find(abs(velaux)<1);
velaux(ff)=1; 
plotvar = log(abs(velaux')).*sign(vel')/log(2);
tickvec =  2.^abs(ticks).*sign(ticks);
end

if choice > 10
plotvar = log(abs(vel'))/log(2);
tickvec= {'1/1024','1/256','1/64','1/16','1/4','4','16','64','256','1024'};
end



h=figure; hold on
worldmap('World')
load coastlines
surfm(yvec,xvec,plotvar)
plotm(coastlat,coastlon,'Color',[0.1, 0.1, 0.1],'linewidth',1.5)
set(gca,'FontSize',15)

colormap([colmatneg(end:-1:1,:); colmatpos])
hcb = colorbar;
hcb.Ticks=ticks;
hcb.TickLabels= tickvec;
hcb.TickLength = 0.02;
caxis([-maxpot maxpot])
xlabel('longitude')
ylabel('latitude')
% set(gca,'ColorScale','log')

%h.Position(1) = h.Position(1)-0.1;
title(titlestring)
set(gca,'FontSize',15)
set(findall(h,'Tag','PLabel'),'visible','off')
set(findall(h,'Tag','MLabel'),'visible','off')


%colmat  (-1 = sea; 13 = antarctica)