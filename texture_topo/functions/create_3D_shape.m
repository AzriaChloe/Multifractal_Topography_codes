function create_3D_shape(topo_out,texture,varargin)
 
%gestion des inputs
%-----------------------------------------------------------------------
p = inputParser;p.StructExpand  = true;p.KeepUnmatched = false; p.FunctionName = 'create_3D_shape'; 

%reduce field for 3D shape
addParameter(p,'reduction_3D',0.02),@(x) isscalar(x) && x>= 0;

%radius of the body 
addParameter(p,'radius',0 ,@(x) isscalar(x) && x>= 0);

%exagerate hillshade (default=1 for realist hillshade)
addParameter(p,'plot_body',1,@(x) isscalar(x) && x>= 0);


%prendre en compte les parametres utilisateurs
parse(p,varargin{:})

%objet contenant tous les parametres
param=p.Results;

%-----------------------------------------------------------------------

%duplique
topo_work=topo_out;

%reduce the dimension of topo
topo_work=imresize(topo_work,param.reduction_3D);

%resolution in latitude and longitude
[nlat,nlong]=size(topo_work);

%matrix of theta and phi value acrosse 2D grid
[theta,phi]=meshgrid(linspace(0,pi,nlat), linspace(0,2*pi,nlong));

%transpose
theta=theta';
phi=phi';

%construct shape in cartesian coordinates
x=topo_work.*cos(phi).*sin(theta);
y=topo_work.*sin(phi).*sin(theta);
z=topo_work.*cos(theta);

figure

% create a spherical shape from cartesiant coordinates
s = surface(x,y,z); 

% add surface texture
s.CData =texture;             
s.FaceColor = 'texturemap';    
s.EdgeColor = 'none';

set(gca,'color','black')
set(gcf,'color','black')
view(15,-28)

axis image

end
    
     

