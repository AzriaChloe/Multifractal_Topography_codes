function [RGB] = texture3D(DEM,param)
%texture and colormap

%construct shading
[H_shadow,H_light] = shade3D(DEM,param);

%colormap
colmapfun=param.colormap;

% nr of colors
nhs = 256;

H_light = gray2ind(H_light,nhs);
H_shadow=gray2ind(H_shadow,nhs);


ncolors = size(colmapfun,1);
cmap = colmapfun;
alims = [min(DEM(:)) max(DEM(:))];
A = gray2ind(mat2gray(DEM,double(alims)),ncolors);

% create colormap for indexing
cmap = cmap(:);
cmap = bsxfun(@times,cmap,linspace(0,1,nhs));
cmap = reshape(cmap,[ncolors 3 nhs]);
cmap = permute(cmap,[3 1 2]);
cmap = reshape(cmap,[ncolors*nhs 3]);

A=imresize(A,size(DEM));
IND  = uint32(H_light+1) + nhs*uint32(A) + 1;
INDshadow  = uint32(H_shadow+1) + nhs*uint32(A) + 1;


cmapUINT8 = uint8(round(cmap*256));
RGB=reshape(cmapUINT8(IND(:),:),[size(IND),3]);
RGBshadow=reshape(cmapUINT8(INDshadow(:),:),[size(INDshadow),3]);

RGBshadow1=RGBshadow;
RGBshadow1(:,:,1)=RGBshadow1(:,:,1)/param.shadowlevel;
RGBshadow1(:,:,2)=RGBshadow1(:,:,2)/param.shadowlevel;
RGBshadow1(:,:,3)=RGBshadow1(:,:,3)/param.shadowlevel;

%recombination of dark and light
RGB=RGBshadow1+RGB;

end
