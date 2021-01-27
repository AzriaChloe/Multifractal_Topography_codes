%rotation of normal vector for projected hemisphere
function [Vall]=normal_rotation(Uall)

[N M b]=size(Uall);

%result
Vall=zeros(size(Uall));

long=linspace(-pi,pi,M);
lat=linspace(-pi/2,pi/2,N);

%meshgrid
[latmat,longmat]=meshgrid(lat,long);

%adjust
latmat=-flipud(latmat');
longmat=-longmat';

%1st rotation along phi (long)
%noral is uy
Vall(:,:,2)=Uall(:,:,2);
Vall(:,:,1)=cos(longmat).*Uall(:,:,1)-sin(longmat).*Uall(:,:,3);
Vall(:,:,3)=sin(longmat).*Uall(:,:,1)+cos(longmat).*Uall(:,:,3);


%2nd rotation along theta (lat) from z to y
%normal is ux
Vall2(:,:,1)=Vall(:,:,1);

Vall2(:,:,2)=sin(latmat).*Vall(:,:,3)+cos(latmat).*Vall(:,:,2);

Vall2(:,:,3)=cos(latmat).*Vall(:,:,3)-sin(latmat).*Vall(:,:,2);

Vall=Vall2;





    

%normalisation of Nall
norm=sqrt(sum(Vall.*Vall,3));
for i=1:3
    Vall(:,:,i)=Vall(:,:,i)./norm;
end


