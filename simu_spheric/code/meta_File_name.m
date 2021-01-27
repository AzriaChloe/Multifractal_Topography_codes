function metarndfilename=meta_File_name(idSeed,nb_lat,nb_long)
    FS2 = '%02.0f'; %format special
    FS4 = '%04.0f';
    FS5 = '%05.0f';
    SidSeed=num2str(idSeed,FS2);
    Snb_lat=num2str(nb_lat,FS4);
    Snb_long=num2str(nb_long,FS5);
    metarndfilename=['/metadata/S' SidSeed '_LAT' Snb_lat '_LON' Snb_long 'METARND.hdf5'];
end    