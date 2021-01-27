%Defintion du bruit de levy
function levyed=levy(alpha, phi, expp)
	phi=-pi/2+pi*phi;
	phi0=(-pi/2)*(1-abs(1-alpha))/alpha;
	levyed=sign(alpha-1)*sin(alpha*(phi-phi0)).*(cos(phi).*abs(alpha-1)).^(-1/alpha).*(cos(phi-alpha*(phi-phi0))./expp).^((1-alpha)/alpha);

    
    

