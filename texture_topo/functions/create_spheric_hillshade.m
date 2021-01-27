function [RGB,topo_work]=create_spheric_hillshade(topo,varargin)
%genere la texture (RGB)
 
%gestion des inputs
%-----------------------------------------------------------------------
p = inputParser;p.StructExpand  = true;p.KeepUnmatched = false; p.FunctionName = 'create_spheric_hillshade'; 

%rotation to apply in longitude degrees
addParameter(p,'longitude_rotation',0,@(x) isscalar(x) && x>= 0 && x<=360);

%reduce field resolution for texture processing (1 for no reduction, 2 for reduce by 2 etc..)
addParameter(p,'reduction',8),@(x) isscalar(x) && x>= 0;

%radius of the body (default is earth in km)
addParameter(p,'radius',6400,@(x) isscalar(x) && x>= 0);

%altitude range (default is earth)
addParameter(p,'altitude_range',20,@(x) isscalar(x) && x>= 0);

%zero level between 0 and 1 (default is 0.5)
addParameter(p,'zero_level',0.5,@(x) isscalar(x) && x>= 0 && x<=1);

%exagerate hillshade (default=1 for realist hillshade)
addParameter(p,'exagerate_shading',1,@(x) isscalar(x) && x>= 0);

%shadow level for darkside
addParameter(p,'shadowlevel',2);

%colormap to apply
addParameter(p,'colormap',jet);

%prendre en compte les parametres utilisateurs
parse(p,varargin{:})

%objet contenant tous les parametres
param=p.Results;

%-----------------------------------------------------------------------

%duplique
topo_work=topo;

%extract lat and long
[nlat,nlong]=size(topo_work);

%rotate longitude
topo_work=circshift(topo_work,[0 round(param.longitude_rotation/360*nlong,0)]);

%set the maximum of topographic fluctation
topo_work=topo_work-min(topo_work(:));
topo_work=topo_work/max(topo_work(:))*param.altitude_range;

%add the radius of the body
topo_work=topo_work+param.radius;

%apply reduction of source field resolution
%topo_work=imresize(topo_work,1/param.reduction);

[RGB]=texture3D(topo_work,param);


end
    
     

