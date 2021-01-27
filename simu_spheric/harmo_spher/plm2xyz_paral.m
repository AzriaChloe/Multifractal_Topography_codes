function varargout=plm2xyz_paral(Lrange,lmcosi,option,degres,c11cmn,lmax,latmax,Plm)
% [r,lon,lat,Plm]=PLM2XYZ(lmcosi,degres,c11cmn,lmax,latmax,Plm)
% [r,lon,lat,Plm]=PLM2XYZ(lmcosi,lat,lon,lmax,latmax,Plm)
%
% Inverse (4*pi-normalized real) spherical harmonic transform.
%
% Compute a spatial field from spherical harmonic coefficients given as 
% [l m Ccos Csin] (not necessarily starting from zero, but sorted), with
% degree resolution 'degres' [default: approximate Nyquist degree].
%
% INPUT:
%
% lmcosi     Matrix listing l,m,cosine and sine expansion coefficients
%            e.g. those coming out of ADDMON
% degres     Longitude/ latitude spacing, in degrees [default: Nyquist] OR
%            "lat": a column vector with latitudes [degrees]
% c11cmn     Corner nodes of lon/lat grid [default: 0 90 360 -90] OR
%            "lon": a column vector with longitudes [degrees]
% lmax       Maximum bandwidth expanded at a time [default: 720]
% latmax     Maximum linear size of the latitude grid allowed [default: Inf]
% Plm        The appropriate Legendre polynomials should you already have them
% 
% OUTPUT:
%
% r          The field (matrix for a grid, vector for scattered points)
% lon,lat    The grid (matrix) or evaluation points (vector), in degrees
% Plm        The set of appropriate Legendre polynomials should you want them
% 
% EXAMPLE:
%
% plm2xyz('demo1') % Illustrates forward and inverse transform
% plm2xyz('demo2',fra) % with 'fra' a data fraction
% plm2xyz('demo3') % Plots EGM96 versus EGM2008 free-air anomaly
% plm2xyz('demo4') % Plots GTM3AR versus EGM2008 topography
% plm2xyz('demo5') % Plots EGM96 geopotential
% plm2xyz('demo6') % EGM2008 topography somewhere
% plm2xyz('demo7') % Tests the expansion to scattered points
% plm2xyz('demo8') % EGM2008 gravity globally
%
% NOTE:
%
% The degree range is now split intelligently in blocks of degrees whose
% memory requirements never exceed the initial expansion from 0 to lmax.
% For very high bandwidth models, specify a small region, and this
% together with 'degres' will determine memory requirements. Built-in
% maxima using 'latmax' and 'lmax'.
%
% NOTE:
%
% Compare to YLM which has a different (-1)^m phase factor.
%
% See also XYZ2PLM, PLM2SPEC, TH2PL, PL2TH, YLM
%
% Special thanks to kwlewis-at-princeton.edu for spotting a bug.
% Last modified by fjsimons-at-alum.mit.edu, 07/13/2012

disp('Reconstruction...')
t0=clock;
% Lowest degree of the expansion
lmin=lmcosi(1);
% Highest degree (bandwidth of the expansion)
L=lmcosi(end,1);
% Never use Libbrecht algorithm... found out it wasn't that good
defval('libb',0)
% Default resolution is the Nyquist degree; return equal sampling in
% longitude and latitude; sqrt(L*(L+1)) is equivalent wavelength
degN=180/sqrt(L*(L+1));
defval('degres',degN);

% When do you get a task bar?
taskmax=100;
taskinf=0;

% Default grid is all of the planet
defval('c11cmn',[0 90 360 -90]);

% Build in maxima to save memory
defval('latmax',Inf); 
defval('lmax',720);

% But is it a grid or are they merely scattered points?
if length(degres)==length(c11cmn)
% It's a bunch of points!
nlat=length(degres);
nlon=length(c11cmn);
% Colatitude vector in radians
theta=[90-degres(:)']*pi/180;
% Longitude vector in radians
phi=c11cmn(:)'*pi/180;

% Initialize output vector
r=repmat(0,nlat,1);

% Now if this is too large reduce lmax, our only recourse to hardcode
ntb=256;
if round(sqrt(nlat)) >= ntb || round(sqrt(nlon)) >= ntb
lmax=round(ntb);
end
elseif length(degres)==1 && length(c11cmn)==4
% It's a grid
if degres>degN
disp('PLM2XYZ: You can do better! Ask for more spatial resolution')
disp(sprintf('Spatial sampling ALLOWED: %8.3f ; REQUESTED: %6.3f',...
degN,degres))
end
% The number of longitude and latitude grid points that will be computed
nlon=min(ceil([c11cmn(3)-c11cmn(1)]/degres+1),latmax);
nlat=min(ceil([c11cmn(2)-c11cmn(4)]/degres+1),2*latmax+1); %%%modif -1 FL

% Initialize output grid
r=repmat(0,nlat,nlon);

% Longitude grid vector in radians
phi=linspace(c11cmn(1)*pi/180,c11cmn(3)*pi/180,nlon);
% Colatitude grid vector in radians
theta=linspace([90-c11cmn(2)]*pi/180,[90-c11cmn(4)]*pi/180,nlat);

%    disp(sprintf('Creating %i by %i grid with resolution %8.3f',nlat,nlon,degres))
else
error('Make up your mind - is it a grid or a list of points?')
end

% Here we were going to build an option for a polar grid
% But abandon this for now
% [thetap,phip]=rottp(theta,phi,pi/2,pi/2,0);

% Piecemeal degree ranges
% Divide the degree range increments spaced such that the additional
% number of degrees does not exceed addmup(lmax)
% If this takes a long time, abort it
els=0; ind=0;
while els<L
ind=ind+1;
% Take positive root
els(ind+1)=min(floor(max(roots(...
[1 3 -els(ind)^2-3*els(ind)-2*addmup(lmax)]))),L);
if any(diff(els)==0)
error('Increase lmax as you are not making progress')
end
end
% Now els contains the breakpoints of the degrees
if ~all(diff(addmup(els))<=addmup(lmax))
error('The subdivision of the degree scale went awry')
end

% Here's the lspacings
if length(els)>2
els=pauli(els,2)+...
[0 0 ; ones(length(els)-2,1) zeros(length(els)-2,1)];
end

for ldeg=1:size(els,1)
    ldown=els(ldeg,1);
    lup=els(ldeg,2);
    % Construct the filename
    fnpl=sprintf('%s/LSSM_back-%i-%i-%i_',fullfile(getenv('IFILES'),'LEGENDRE'),ldown,lup,nlat);


    % ONLY COMPLETE LINEARLY SPACED SAMPLED VECTORS ARE TO BE SAVED!
    if [exist([myname(fnpl,lup), '.mat'],'file')==2 && length(c11cmn)==4 && all(c11cmn==[0 90 360 -90])]...
        & ~[size(els,1)==1 &&  exist('Plm','var')==1]
        disp('Precomputed Plm will be used...')
    elseif option==0
        % Get Legendre function values at linearly spaced intervals
        % disp(sprintf('Using preloaded %s',fnpl))
        %load(fnpl)
        % AND TYPICALLY ANYTHING ELSE WOULD BE PRECOMPUTED, BUT THE GLOBAL
        % ONES CAN TOO! The Matlabpool check doesn't seem to work inside 
        %elseif size(els,1)==1 &&  exist('Plm','var')==1 && matlabpool('size')==0
        % disp(sprintf('Using precomputed workspace Legendre functions'))
        %else
        % Evaluate Legendre polynomials at selected points
        disp('computing and saving legendre polynomials')

        for l=ldown:lup
            Plm=(legendre(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
            save(myname(fnpl,l),'Plm') 
        end

    elseif option==1
        disp('calcul a la vole')  
    end

    % Loop over the degrees
    more off
    %    disp(sprintf('PLM2XYZ Expansion from %i to %i',max(lmin,ldown),lup))
    for l=max(lmin,ldown):lup
    
        if length(find(Lrange==l))>0
		l
            % Compute Schmidt-normalized Legendre functions at 
            % the cosine of the colatitude (=sin(lat)) and 
            % renormalize them to the area of the unit sphere

            % Remember the Plm vector always starts from ldown
            b=addmup(l-1)+1-addmup(ldown-1);
            e=addmup(l)-addmup(ldown-1);

            if option==0
            load(myname(fnpl,l))
            else
            Plm=(legendreFL1(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';    
            end

            plm=Plm';

            m=0:l;
            mphi=m(:)*phi(:)';

            defval('tst',0)
            if tst
                f=(plm.*plm)'.*repmat(sin(theta(:)),1,l+1);
                c=(cos(mphi).*cos(mphi))';
                ntest=simpson(theta,f).*simpson(phi,c)/4/pi;
                %disp(sprintf('Mean Normalization Error l= %3.3i: %8.3e',l,(sum(abs(1-ntest)))/(l+1)))
                % For a decent test you would use "legendreprodint"
            end

            % Find the cosine and sine coefficients for this degree
            clm=shcos(lmcosi,l);
            slm=shsin(lmcosi,l);

            fac1=repmat(clm,1,nlon).*cos(mphi)+...
            repmat(slm,1,nlon).*sin(mphi);

            % Sum over all orders and (through loop) over all degrees
            if length(degres)==length(c11cmn)
                expa=sum(plm.*fac1,1)'; 
                % Or diag(plm'*fac1) if you will
            elseif length(degres)==1 & length(c11cmn)==4
                expa=plm'*fac1;
            end

            r=r+expa;
        end
    end
end

lon=phi*180/pi;
lat=90-theta*180/pi;

% Prepare output
vars={r,lon,lat,Plm};
varargout=vars(1:nargout);


function name=myname(fnpl,l)
name =[fnpl, 'L_',num2str(l)];














