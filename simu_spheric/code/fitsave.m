function fitsave(path, height, idSeed, Seed, alpha, H, C1)
import matlab.io.* ;

%Fits convention: row number increases upward

nb_lat=size(height,1);
nb_long=size(height,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%location of the file
filename=[path File_name(idSeed,nb_lat,nb_long,alpha,C1,H) '.fits'];
delete(filename);%delete file if exists

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calc for further parameter allocations
increment_lat=180./nb_lat;
increment_long=360./nb_long;

%coordinates of repoint 
%(in pixel)
reflong=(nb_long)/2;
reflat=(nb_lat)/2;
%(in degree)
crval1=reflong*(360./nb_long)-180.;
crval2=reflat*(180./nb_lat)-90.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%create and allocate fits file
fptr=fits.createFile(filename);
fits.createImg(fptr,'double_img',[nb_lat nb_long]);
fits.writeImg(fptr,height);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Keywords
fits.writeKey(fptr,'DATAMAX',max(max(height)));
fits.writeKey(fptr,'DATAMIN',min(min(height)));


%Parameters
fits.writeComment(fptr,'           PARAMETERS');
fits.writeComment(fptr,['nb_lat = ' num2str(nb_lat)  ' : Latitude sampling']);
fits.writeComment(fptr,['nb_long = ' num2str(nb_long) ' : Longitude sampling']);
fits.writeComment(fptr,['idSeed = ' num2str(idSeed) ' : Random Seed identifier']);
fits.writeComment(fptr,['Seed = ' num2str(Seed) ' : Random Seed Value']);
fits.writeComment(fptr,['C1 = ' num2str(C1) ' : degree of intermittency']);
fits.writeComment(fptr,['alpha = ' num2str(alpha) ' : degree of multifractality']);
fits.writeComment(fptr,['H = ' num2str(H) ' : degree of smoothness']);


%Coordinates description
fits.writeKey(fptr,'WCSNAME','PLANETOCENTRIC COORDINATES');
fits.writeKey(fptr,'WCSAXES',int32(2));

fits.writeKey(fptr,'CUNIT1','deg     ');
fits.writeKey(fptr,'CUNIT2','deg     ');

fits.writeKey(fptr,'CTYPE1','RA---CAR','Longitude-cylindrical proj -180W +180E ');
fits.writeKey(fptr,'CTYPE2','DEC--CAR','Latitude-cylindrical proj -90S +90N ');


%Reference point
fits.writeKey(fptr,'CRPIX1',reflong+0.5 ,'horizontal position of reference point(pixel)');
fits.writeKey(fptr,'CRPIX2',reflat+0.5 ,'vertical position of reference point(pixel)');

fits.writeKey(fptr,'CRVAL1',crval1,'longitude of reference point (deg)');
fits.writeKey(fptr,'CRVAL2',crval2,'latitude of reference point (deg)');


%Conversion Matrix 
fits.writeKey(fptr,'CD1_1', increment_long );
fits.writeKey(fptr,'CD1_2', 0);
fits.writeKey(fptr,'CD2_1', 0);
fits.writeKey(fptr,'CD2_2',increment_lat);

%Other method to define Conversion: rotation PCi_j and increment CDELTi
% fits.writeKey(fptr ,'CDELT1',increment_long,'longitude increment/pixel'); 
% fits.writeKey(fptr ,'CDELT2',increment_lat,'latitude increment/pixel'); 
% fits.writeKey(fptr,'PC1_1',1 );
% fits.writeKey(fptr,'PC1_2',0);
% fits.writeKey(fptr,'PC2_1',0);
% fits.writeKey(fptr,'PC2_2',1);

fits.closeFile(fptr);
%fitsdisp(filename);



