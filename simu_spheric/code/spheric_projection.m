%convert field lat long matrix to a spherical projection
function DEM=spheric_projection(field,longmin,wantedsquaresize)
[A B]=size(field);

%radius of planet
radius=B/2;

%select only one visible half of the sphere
field1=[field field];
half=field1(:,longmin:longmin+radius-1);

%size of initial matrix;
[nlat nlong]=size(half);

%ini frame

%frame=zeros(wantedsquaresize,wantedsquaresize,90);

%length of original equator
Lequ=length(half(1,:));

%define grid cartesian coordinates
x=linspace(-1,1,wantedsquaresize);
y=linspace(1,-1,wantedsquaresize);

%meshgrid
x=meshgrid(x);
y=meshgrid(y)';
z=sqrt(1-x.^2-y.^2);
z(~imag(z)==0)=NaN;

%create theta and phi mat
thetamat=acos(y);
phimat=asin(x./sin(thetamat))+pi/2;
phimat(~imag(phimat)==0)=NaN;

%adpat to origial grid
coordlat=floor(thetamat/pi*Lequ);
coordlong=floor(phimat/pi*Lequ);


%replace zero value by 1
coordlat(coordlat==0)=1;
coordlong(coordlong==0)=1;

%reshape
rcoordlat=reshape(coordlat,[wantedsquaresize^2,1]);
rcoordlong=reshape(coordlong,[wantedsquaresize^2,1]);
rc=(rcoordlong-1)*nlat + rcoordlat;

%select altitude values
rc(isnan(rc))=1;
selectal=half(rc);

%reshape alitude matrix
proj=reshape(selectal,[wantedsquaresize,wantedsquaresize]);

dem=proj;
X=1:length(dem(1,:));
Y=1:length(dem(:,1));
% % 
DEM = GRIDobj(X,Y,dem);




        
        
        
