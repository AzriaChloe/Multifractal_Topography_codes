function filename=File_name(idSeed,nb_lat,nb_long,alpha,C1,H)
FS2 = '%02.0f'; %format special
FS4 = '%04.0f';
FS5 = '%05.0f';


SidSeed=num2str(idSeed,FS2);
Snb_lat=num2str(nb_lat,FS4);

if nb_lat==5001
   Snb_long=num2str(nb_long,FS5);
   Salpha=num2str(alpha*100,FS2); 
else   
   Salpha=num2str(alpha*10,FS2);
   Snb_long=num2str(nb_long,FS4);
end
SC1=num2str(C1*100,FS2);
SH=num2str(H*10,FS2);

filename=[ '/S' SidSeed '_A' Salpha '_C' SC1 '_H' SH '_LAT' Snb_lat '_LON' Snb_long];
