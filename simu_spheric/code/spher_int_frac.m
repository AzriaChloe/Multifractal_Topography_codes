function [resultat,C]=spher_int_frac(champ,ordre)
% ---------------------------
%
% fonction: spher_int_frac
%
% But de la fonction : integration fractionnaire spherique
% 
% Input : champ= matrice lat/long, ordre= ordre de l'integration fractionnaire
%
% Output : resultat = champ spherique integre, lat long
%
% ---------------------------  

%dimension du champ initial
dim=length(champ(:,1));

%nombre de divisions en latitude
nb_lat=length(champ(:,1));

%nombre de divisions en longitude
nb_long=length(champ(1,:));

%----------------------------------------------------------
% A singularite

% definition de la singularite (vecteur lat)
xx=2*(1:nb_lat)-1;
singularity_vecteur=xx.^(-(2-ordre));

%definition de la singularite (matrice lat long)	
singularity_matrice=zeros(size(champ));
for j=1:length(champ(1,:)) 
    singularity_matrice(:,j)=singularity_vecteur;
end

%----------------------------------------------------------
% Decomposition de la singularite en harmoniques spherique

%decomposition
option_nosave = 1; % 0 = precompute and save the Legendre polynomial

[sigma_l_m]=xyz2plm_modified(singularity_matrice,dim-2,'gl', option_nosave);

%Duplication des sigma_l_m pour toutes les valeurs de m (dans le filtre, la valeur du couple(l,m=0) sera associee e tous les couples (l,m) du champ
for k=1:length(sigma_l_m(:,1))

	% si la valeur de m est nulle
    if sigma_l_m(k,2)>0
	
		% dupliquer le coefficient m=0 pour tout m different de 0
		sigma_l_m(k,3)=sigma_l_m((k-1),3);
		sigma_l_m(k,4)=sigma_l_m((k-1),4);
    end
end

%l
l=sigma_l_m(:,1);

%m
m=sigma_l_m(:,2);

%coefficients de la decomposition de la singularite en harmoniques circulaires (+ les duplications)
sigma_l=[sigma_l_m(:,3)+1i*sigma_l_m(:,4)];

%----------------------------------------------------------
% Decomposition du champ en harmoniques spherique

%decomposition
ulm=xyz2plm_modified(champ,dim-2,'gl', option_nosave);
ulm=ulm(:,3)+1i*ulm(:,4);

%----------------------------------------------------------
%application du filtre dans l'espace des harmoniques spheriques

%filtre
filtrage=sigma_l.*ulm.*(sqrt((4*pi)./(2*l+1))); 

%Remise en dimension nlignes*4colonnes pour l'utilisation de plm2xyz
C=[l,m,real(filtrage),imag(filtrage)];


%----------------------------------------------------------
%retour dans l'espace reel
resultat=plm2xyz_modified(C,option_nosave);


%----------------------------------------------------------
%application d'un coefficient empirique
resultat=0.4*dim*dim*2*resultat;




