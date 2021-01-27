
close all;
clear all;

addpath(genpath('./'))
path=['data'];
pathtextures=[path '/Textures'];

file=['/S07_A20_C05_H04_LAT0361_LON0720'];
filename=[path file 'TOPO.fits']


%Charger la topo
topo=fitsread(filename);
topo=flip(topo);%because fits record have by convention rows increasing upward and png have rows increasing downward
nb_lat=size(topo,1);
nb_long=size(topo,2);
mintopo = min(min(topo));
maxtopo = max(max(topo));
topo = (topo - mintopo)/(maxtopo - mintopo); % rescale between 0 and 1


disp('Earth like');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intermediane=reshape(topo,[1 nb_lat*nb_long]);
mediane=median(intermediane)

cmap1=cmocean_modified('to2','pivot',mediane);
cmap=imresize(cmap1,[1000,3]);
cmap(cmap<0)=0.;

RGB_topoflat=colorflat(topo,cmap);
image(RGB_topoflat);
imwrite(RGB_topoflat,[pathtextures file  'COLOFLAT.png']);

%generer texture
[texture,topo_out]=create_spheric_hillshade(topo,'reduction',4,...
    'radius',6400,'altitude_range',170,'longitude_rotation',0,...
'exagerate_shading',7,'colormap',cmap,...
'shadowlevel',3);

imwrite(texture,[pathtextures file  'COLOSHAD.png']);
image(texture);

%plot surface 3D
create_3D_shape(topo_out,texture,'reduction_3D',0.16)
view(-128,39)   
set(gca, 'Visible', 'off')
saveas(gcf,[pathtextures file 'COLOTOPO.png']);          



disp('Small Body');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RGB_grayflat=zeros([size(topo) 3])+round(0.8*256);
RGB_grayflat=uint8(RGB_grayflat);
imwrite(RGB_grayflat,[pathtextures file 'GRAYFLAT.png']);

cmap=[.8,.8,.8];
%generer texture
[texture,topo_out]=create_spheric_hillshade(topo,'reduction',4,...
    'radius',80,'altitude_range',40,'longitude_rotation',0,...
'exagerate_shading',2,'colormap',cmap,'shadowlevel',3);

imwrite(texture,[pathtextures file 'GRAYSHAD.png']);
image(texture);

%plot surface 3D
create_3D_shape(topo_out,texture,'reduction_3D',0.16)
view(-128,39) 
set(gca, 'Visible', 'off')
saveas(gcf,[pathtextures file 'GRAYTOPO.png']);
  

    
function RGB=colorflat(topo,cmap)
    sz=size(cmap,1);
    cmap = round(cmap*256);
    %topo must be previously resized between 0 and 1
    dz=1/sz;
    RGB=zeros([size(topo) 3]);

    for i=1:sz
        mask1 = topo<i*dz ;
        mask2= topo<(i-1)*dz;
        mask = mask1-mask2;
        RGB(:,:,1)=RGB(:,:,1)+mask*cmap(i,1);
        RGB(:,:,2)=RGB(:,:,2)+mask*cmap(i,2);
        RGB(:,:,3)=RGB(:,:,3)+mask*cmap(i,3);
        clearvars -except i RGB dz cmap sz topo
    end  
    RGB=uint8(RGB);
end