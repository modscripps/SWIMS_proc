function [ThRef] = thalweg_refer(lons, lats, ThalFile, ThVars, minkm, maxkm)
% thalweg_refer.m - find nearest points on thalweg to given locations.
% ThRef = thalweg_refer(lons, lats, ThalFile, ThVars, minkm, maxkm);
%   lons, lats = specified positions to reference to thalweg (req'd)
%   ThalFile = name of thalweg file (req'd);
%   ThVars = {'lonname','latname','distname','anglename'} = variable/field
%       names from thalweg file to use for calculation (angle is optional),
%       Default = {'TW.lon_hires','TW.lat_hires','TW.dist_hires','TW.ang_hires'};
%       (Use full structure reference if appropriate)
%   minkm,maxkm = optional along-thalweg range.
% Outputs the following for each closest thalweg reference point:
%   ThRef.thal_lon, ThRef.thal_lat = position
%   ThRef.thal_dist = reference distance along thalweg
%   ThRef.thal_ang = thalweg angle at reference (from thalweg file)
%   ThRef.thal_sepkm = distance in km to thalweg
%   ThRef.ang_to_thal = angle from point to thalweg (0=E,90=N,180=W,-90=S)
% All fields are row vectors
% Dave Winkel, 25-Jan-2002

ThRef=[];
if nargin<3
	error(['Inputs lons,lats,ThalFile are required.']);
end
if nargin<4 | isempty(ThVars)
    ThVars = {'TW.lon_hires','TW.lat_hires','TW.dist_hires','TW.ang_hires'};
end
if length(ThVars)<3 | length(ThVars)>4
    error('ThVar requires variable names for thalweg lon,lat,distance,angle');
end
if nargin<5 | isempty(minkm)
   minkm = -inf;
end
if nargin<6 | isempty(maxkm)
   maxkm = inf;
end
if maxkm<=minkm
    error('minkm,maxkm out of order')
end

npt = size(lats);
x=size(lons);
if min(npt)~=1 | max(npt)<1 | ~isequal(npt,x)
   error('LONS and LATS must be identically sized vectors')
end
if npt(1)>1, lats=lats'; lons=lons'; end  % force row vectors
npt = max(npt);
inbad = find(isnan(lats)|isnan(lons));

% Load thalweg points
TH = load(ThalFile);
vnams = {'lon_th','lat_th','dis_th','ang_th'};
lon_th=[]; lat_th=[]; dis_th=[]; ang_th=[];
% Retrieve thalweg values from specified file (and variable/field names)
for iv=1:length(ThVars)
    Vst=[]; Vnm=ThVars{iv}; ip=findstr(Vnm,'.');
    if ip>1 & ip<length(Vnm) % structure field in TW file was specified
        Vst = ['.' Vnm(1:ip-1)]; Vnm = Vnm(ip+1:end);
    end
    ip = eval(['isfield(TH' Vst ', Vnm)']);
    if ip
        eval([vnams{iv} ' = TH' Vst '.' Vnm ';']);
    else
        error(['variable/field=' ThVars{iv} ' not in thalweg file'])
    end
end
if isempty(ang_th), ang_th = NaN*dis_th; end % angle is optional
% get subset of points in  distance range
idd = find(dis_th >= minkm & dis_th <= maxkm);
lat_th=lat_th(idd); lon_th=lon_th(idd);
dis_th=dis_th(idd); ang_th=ang_th(idd);
if size(lat_th,1)>1 % force row vectors
    lat_th=lat_th'; lon_th=lon_th';
    dis_th=dis_th'; ang_th=ang_th';
end
nth = length(idd);

clear TH

M_PER_NM=1.852e3; % meters per nautical mile
m_per_dlat=60*M_PER_NM; % meters per degree of latitude
% use average latitude to determine meters per degree of longitude
lon_per_lat = cos(mean(lat_th)*pi/180);
m_per_dlon = m_per_dlat * lon_per_lat;

% For every input point, compute distance (in km) to ALL thalweg points
clear dx dy
for ip = 1:npt
   dx(:,ip) = lon_th' - lons(ip);
   dy(:,ip) = lat_th' - lats(ip);
end
ds = sqrt( (dx*lon_per_lat).^2 + dy.^2 ) * m_per_dlat/1000;

% find the two closest thalweg points
[dn1,it1] = min(ds); % finds closest for each one
% look only at neighboring points for next closest
for ip = 1:npt
   il=it1(ip)-1; ig=it1(ip)+1;
   if il<1  % start of tw range
      it2(ip)=ig; dn2(ip)=ds(ig,ip);
   elseif ig>nth  % end of tw range
      it2(ip)=il; dn2(ip)=ds(il,ip);
   else  % choose the closer one 
      it2(ip)=il; dn2(ip)=ds(il,ip);
      if ds(ig,ip) < ds(il,ip)
         it2(ip)=ig; dn2(ip)=ds(ig,ip);
      end
   end   
end  % now have two closest distances, indices

% interpolate between points: use law of cosines, but allow for
% cases where thalweg changes direction
Dt = abs( dis_th(it1) - dis_th(it2) ); % separation along thalweg
delt = 0.5 * (Dt.^2 + dn1.^2 - dn2.^2) ./ Dt;
delt(find(delt<0)) = 0; % point near outside of corner (obtuse)
delt = delt .* sign(it2-it1); % move up or down the thalweg
% Save in output structure
ThRef.thal_dist = dis_th(it1) + delt;
rat = delt./Dt;
%keyboard
ThRef.thal_lat = interp1(dis_th, lat_th, ThRef.thal_dist);
ThRef.thal_lon = interp1(dis_th, lon_th, ThRef.thal_dist);
ThRef.thal_ang = ang_th(it1);  % use nearest refn point to avoid +/-PI problem
if ~isempty(inbad) % Make sure input (lat,lon)=NaN are NaN on output
   ThRef.thal_dist(inbad)=NaN;
   ThRef.thal_lat(inbad)=NaN;
   ThRef.thal_lon(inbad)=NaN;
   ThRef.thal_ang(inbad)=NaN;
end
% compute distance,angle(0=E,90=N,180=W,-90=S) TO thalweg FROM points
clear dx dy ds
dy = ThRef.thal_lat - lats;
dx = (ThRef.thal_lon - lons) * lon_per_lat;
ThRef.thal_sepkm = sqrt( dx.^2 + dy.^2 ) * m_per_dlat/1000;
ThRef.ang_to_thal = atan2(dy,dx) * 180/pi;

return