function [H_shadow,H_light] = shade3D(DEM,param)


%extract lat and long
[nlat,nlong]=size(DEM);

%matrix of theta and phi value acrosse 2D grid
[theta,phi]=meshgrid(linspace(pi/nlat,pi-pi/nlat,nlat), linspace(0,2*pi,nlong));


%transpose
theta=theta';
phi=phi';

m=zeros(nlat,nlong,3);
m(:,:,1)=theta;
m(:,:,2)=phi;
m(:,:,3)=DEM;


res=reshape(m,[nlat*nlong,3]);
%Rdtheta acrosse DEM
Rdtheta=DEM*pi/nlat;

%Rdphicostheta theta acrosse DEM
Rdphicostheta=DEM.*2*pi/nlong.*sin(theta);
%

%calculate surface normals
[Nx,Ny,Nz] = surfnorm(DEM*param.exagerate_shading);

%transform Nx et Ny to take into account the size of each cell
Ny=Ny./Rdtheta;
Nx=Nx./Rdphicostheta;


%write normal vector in one mat
Nall=zeros([size(Nx) 3]);
Nall(:,:,1)=Nx;
Nall(:,:,2)=Ny;
Nall(:,:,3)=Nz;

%normalisation of Nall
norm=sqrt(sum(Nall.*Nall,3));
for i=1:3
    Nall(:,:,i)=Nall(:,:,i)./norm;
end


%perform rotation
Vall=normal_rotation(Nall);

%force solar inclination
sx=0 ;
sy=0;
sz=1;
 
%perform scalar product
H =Vall(:,:,1)*sx + Vall(:,:,2)*sy + Vall(:,:,3)*sz;

%dark and light side
H_light= uint8(H*255);
H_shadow=uint8(-H*255);

end