function [Koppen_new] = Koppen_formatter(Koppenraw,Nx,Ny,doland)
%re-orders Koppen such as to obtain the order: 

 % 0 sea
 % 1 tropical rainforest
 % 2 subtropical forest
 % 3 tropical grassland and savannah
 % 4 deserts and xeric shrublands
 % 5 mediterranean shrublands
 % 6 warm temperate and coastal forests
 % 7 temperate wet forests (small group)
 % 8 temperate grassland and savannah
 % 9 various forests (fairly temperate and continental)
 % 10 montane and arctic grass and shrubland
 % 11 boreal coniferous forest taiga
 % 12 boreal and mountain tundra
 % 13 ice cap




Koppen1=reshape(Koppenraw,Nx,Ny)';

Koppen2=Koppen1(Ny:-1:1,:);
Koppen=Koppen2(:,[Nx/2+1:Nx,1:Nx/2]);

clear Koppen1; clear Koppen2
%Re-organise some climate zones
%NEW name: -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13
%OLD name:  0   1 2 7 8 9 6 4 5 3 10 11 12 13

ff=find(Koppen==-1); Koppen(ff)= 0.1;
  if doland==1
ff=find(Koppen== 7); Koppen(ff)= 3.1;
ff=find(Koppen== 8); Koppen(ff)= 4.1;
ff=find(Koppen== 9); Koppen(ff)= 5.1;
ff=find(Koppen== 4); Koppen(ff)= 7.1;
ff=find(Koppen== 5); Koppen(ff)= 8.1;
ff=find(Koppen== 3); Koppen(ff)= 9.1;
  end 
Koppen_new=round(Koppen);

end