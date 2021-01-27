function resultat=spher_sim(alpha,C1,H, phi, expp)
% ---------------------------
%
% fonction: spheric_simulation
%
% But de la fonction : genere une carte lat long d'un champ FIF spherique
%
% Note :
% 
% Input : nb_lat  (nb de pas en latitude
% Input : nb_long  (nb de pas en longitude
% Input : alpha  degre de multifractalite
% Input : C1, degre d'intermittence
% Input : H, degre de lissage
%
% Output : resultat = FIF spherique
%
% ---------------------------  

%generation du bruit

bruit=levy(alpha,phi,expp);

%1ere integration fractionnaire
[bruit_integrated,C]=spher_int_frac(bruit,2*(1-1/alpha));


%exponentialisation
epsed=expmulti(bruit_integrated,C1,0.75/sqrt(2),alpha);

%2eme integration fractionnaire
resultat=spher_int_frac(epsed,H);




    
	


