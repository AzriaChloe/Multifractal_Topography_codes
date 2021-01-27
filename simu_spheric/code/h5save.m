function h5save(path, nb_lat, nb_long, height , lat, long, idSeed, Seed, alpha, H, C1, radius, phi, expp )

% metarndfilename=[path meta_File_name(idSeed,nb_lat,nb_long)]
% delete (metarndfilename);
% h5create(metarndfilename,'/expp', [nb_lat nb_long]);
% h5create(metarndfilename,'/phi', [nb_lat nb_long]);
% h5write(metarndfilename,'/expp', expp);
% h5write(metarndfilename,'/phi', phi);
% h5writeatt(metarndfilename,'/expp','idSeed',idSeed);
% h5writeatt(metarndfilename,'/expp','Seed',Seed);


filename=[path File_name(idSeed,nb_lat,nb_long,alpha,C1,H) '.hdf5'];

%supprime les fichiers h5 si ils existent déjà
delete (filename) ; 

h5create(filename ,'/data',[nb_lat nb_long]);
h5create(filename,'/lat',[nb_lat]);
h5create(filename, '/long',[nb_long]);

h5write(filename,'/data', height );
h5write(filename,'/lat',lat);
h5write(filename,'/long',long);

h5writeatt(filename,'/data','C1',C1);
h5writeatt(filename,'/data','alpha',alpha);
h5writeatt(filename,'/data','H',H);
h5writeatt(filename,'/data','radius',radius);
h5writeatt(filename,'/data','idSeed',idSeed);

h5writeatt(filename,'/data','nb_lat',nb_lat);
h5writeatt(filename,'/data','nb_long',nb_long);
h5writeatt(filename,'/data','Seed',Seed);



