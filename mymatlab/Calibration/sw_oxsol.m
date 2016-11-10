function oxs = sw_oxsol(T,S)
% oxs = sw_oxsol(T,S); Oxygen saturation value (ml/l) for seawater at T,S and 1 atm
%	T = in situ temperature (deg C);
%	S = salinity (PSU);
%	oxs = oxygen saturation value (ml/l), figured at 1 atmosphere of pressure
%
%	Test values: sw_oxsol([4 20], [0.010 0.035]) = [8.5754 5.1749];
%	Coded 25-Jan-2012 by David Winkel, per REVISED SeaBird Appl. Note 64
%		from Garcia and Gordon, 1992 (L&O)
%
% Earlier version coded in sw_oxsat.m (also in sw_satO2.m) was from
% 2001 version of Note 64, using Weiss, 1970 (DSR)

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('sw_oxsol.m: Must pass 2 parameters')
end %if

oxs = NaN*T; % initialize

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   warning('check_stp: S & T must have same dimensions')
end %if

% Formula coeffecients:
A0 = 2.00907;
A1 = 3.22014;
A2 = 4.0501;
A3 = 4.94457;
A4 = 0.256847;
A5 = 3.88767;
B0 = -0.00624523;
B1 = -0.00737614;
B2 = -0.010341;
B3 = -0.00817083;
C0 = -4.88682e-7;

%% Compute:

Ts = log( (298.15-T)./(273.15+T) );
oxs = exp( (A0+Ts.*(A1+Ts.*(A2+Ts.*(A3+Ts.*(A4+Ts.*A5))))) + ...
    S.*(B0+Ts.*(B1+Ts.*(B2+Ts.*B3))) + C0.*S.^2);
