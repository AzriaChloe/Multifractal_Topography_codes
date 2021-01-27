%exponentialisation du bruit de levy integre (factor est determine empiriquement)
function out=expmulti(bruit,C1,factor,alpha)
    bruit_multifract=exp(factor*(C1*2/pi)^(1/alpha)*bruit);
    ff=mean(mean(bruit_multifract));
    out=bruit_multifract/ff;