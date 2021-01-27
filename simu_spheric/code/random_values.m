%Générer les nombres aléatoires phi et expp
function [phi,expp] = random_values(a,b,Seed)
    rng(Seed); %random seed determinée
    %rng('shuffle');%tirage aléatoire sur l'horloge
    phi=rand(a,b);%matrice distrib aleatoire entre 0 1 uniforme; dim a,b
    expp=exprnd(1,a,b);%matrice champ aleatoire distibution exp; moyenne=1, a,b dim matrice  
    %info=rng() ;