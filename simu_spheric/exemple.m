
%Exemple de generation et representation d'un champ FIF spherique
close all
clear all

%inclure les fonctions des sous-repertoires
addpath(genpath('./'))

%chemin pour l'enregistrement des données de sortie
mkdir('data');
path='data'; 

%----Parametres------------------------------------------

% latitude
nb_lat=361;
% longitude
nb_long=720;

%degre de multifractalite
alpha=2.0;

%degre de lissage
H=0.4;

%degre d'intermittence

C1=0.05;

%identifiant de la random Seed
idSeed=7;

%valeur de la random seed
Seed=1060060143;%déterminé ; disp(rng('Shuffle')
%OU
%Seed=rng('shuffle');%tirage aléatoire sur l'horloge

filename=File_name(idSeed,nb_lat,nb_long,alpha,C1,H)
% disp(filename);

%Génération des matrices de nombres aléatoires à partir de la random seed
[phi,expp] = random_values(nb_lat,nb_long,Seed); % déterminé par la seed ou de l'horloge

%----generation du FIF----------------------------------------

%simulation
FIF=spher_sim(alpha, C1 ,H , phi, expp);%tableau de taille nb_lat*nb_long contenant la topo

minFIF = min(min(FIF));
maxFIF = max(max(FIF));

FIFnorm = (FIF - minFIF)/(maxFIF - minFIF); % rescale between 0 and 1
height=FIFnorm*2000-1000; %rescale between -1000 and +1000

%Save the result
fitsave(path, height, idSeed, Seed, alpha, H, C1); %in fits

%----plot--------------------------------------
%2D
figure; 
imagesc(height)% affichage de la matrice avec une colorbar

% %3D
% %rayon de la planète utilisé pour le plot 3D
% radius = 10000;
% %generate 3D points
% 
% lat=linspace(-90,90,nb_lat);
% long=linspace(-180,180,nb_long);
% [longgrid, latgrid] = meshgrid(long, lat);
% 
% latpoint = reshape(latgrid, 1, nb_lat*nb_long);
% longpoint = reshape(longgrid, 1, nb_lat*nb_long);
% heightpoint = reshape(height, 1, nb_lat*nb_long);
% 
% %points = [latpoint; longpoint; heightpoint+radius];
% X = (heightpoint+radius).*cos(latpoint*pi/180).*cos(longpoint*pi/180);
% Y = (heightpoint+radius).*cos(latpoint*pi/180).*sin(longpoint*pi/180);
% Z = (heightpoint+radius).*sin(latpoint*pi/180);
% %points_XYZ = [X;Y;Z];
%
% pcshow([X(:),Y(:),Z(:)],heightpoint)%nuage de pts avec une couleur en fctn de l'altitude
%
%
% %%%approximation en ensemble de triangles (chronophage) 
%
% % dt = delaunayTriangulation(points_XYZ');
% % tetramesh(dt, 'FaceColor', 'cyan');
% % 
% % triboundary = convexHull(dt)% enveloppe convexe de dt
% % trisurf(triboundary, X(:),Y(:),Z(:), 'FaceColor', 'cyan') %affiche des surfaces
