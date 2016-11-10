function [CycD] = SWIMS_Cycle_Times(Ship_speed, Line_rates, Depths, Line_Max, Flow)
% Usage: [CycD] = SWIMS_Cycle_Times(Ship_speed,Line_rates, Depths);
% With inputs: (DEF = defaults)
%   Ship_speed = Vessel speed in Knots while SWIMS in being cycled. DEF=4;
%   Line_Rates = [Out, In] in m/min for cycling SWIMS. DEF=[100,70];
%   Depths = column vector of depths/m to be reached by SWIMS. DEF=[0:100:700]';
%   Line_max = amount of seacable on drum. DEF=700 m;
%   Flow = current speed in Knots with(+) or against(-) ship. DEF=0; 
% Returned fields in structure = CycD:
%   Speed,Rates,Line_max,Flow = same as Ship_speed,Line_rates,Line_max,Flow;
%   Depths = same as input=Depths up to depth that Line_max is reached,
%      then maximum depth range followed by NaNs;
%   Line_out = vector with amount of line [m] that SWIMS takes to
%      reach each Depths(:) at given cycling rate and tow speed;
%   Cycle_time = vector of elasped time [s] to cycle SWIMS down
%      and back up to each Depths(:);
%   Down_time = same as Cycle_time, but only for downward portion;
%   Cycle_distance = distances [m] that vessel travels for one cycle.
% Dave W, 6/21/06

clear CycD
%% Defaults: %%
Deps = [0:100:700]'; % depths/m that SWIMS attains (table columns, also)
CycD.Speed = 4; % ship speed [knots]
CycD.Rates = [100, 70]; % line rate in,out [m/min]
CycD.Depths = Deps;
CycD.Line_max = 700; % maximum line out on drum [m]
CycD.Flow = 0; % absolute water flow speed [knots], with(+) or against(-) ship
%% Check for input arguments: %%
if nargin>0 & ~isempty(Ship_speed)
    if Ship_speed>=0 & Ship_speed<=6
        CycD.Speed = Ship_speed;
    else
        warning('Ship_speed out of range (0-6), using default 4 kts')
    end
end
rmxo=140;rmxi=100;
if nargin>1 & ~isempty(Line_rates)
    if Line_rates(1)>0
        CycD.Rates(1) = min(Line_rates(1), rmxo);
        if Line_rates(1)>rmxo
            warning(['Max rate out is ' int2str(rmxo) '  m/min, using that'])
        end
    end
    if length(Line_rates)>1
        if Line_rates(2)>0
            CycD.Rates(2) = min(Line_rates(2), rmxi);
            if Line_rates(2)>rmxi
                warning(['Max rate in is ' int2str(rmxi) '  m/min, using that'])
            end
        end
    end
end % of line rate checking
if nargin>2 & ~isempty(Depths)
    if size(Depths,2)>1 & size(Depths,1)==1
        Depths = Depths'; % coulmn vector
    end
    x = Depths; Depths = sort(x);
    CycD.Depths = Depths;
end
if nargin>3 & ~isempty(Line_max)
    if Line_max>100 & Line_max<1000
        CycD.Line_max = Line_max;
    end
end
if nargin>4 & ~isempty(Flow)
    if abs(Flow)<6
        CycD.Flow = Flow;
    end
end
% Finished with input

%%%  Lookup Tables %%%%%%%%%%%%%%%%%%%%%%
% Deps: specified above with 'Defaults:'
Spds = [0, 2, 3, 4, 5, 6]; % ship speed, in knots (relative flow, actually)
% Arrays of LineOut/m by Speed,Depth;  for 100 m/min, 120 m/min line rates
Lmx = 2000; % maximum for extrapolating into unknown/invalid realm
Lrates = [100, 120];
LD_SbyD100 = ...
     [0, 0, 0, 0, 0, 0;
     100, 100, 100, 100, 110, 120;
     200, 200, 230, 250, 300, 350;
     300, 300, 360, 450, 530, 640;
     400, 400, 490, 660, 760, 930;
     500, 500, 620, 870, Lmx, Lmx;
     600, 600, 750, Lmx, Lmx, Lmx;
     700, 700, Lmx, Lmx, Lmx, Lmx];
 LD_SbyD120 = ...
     [0, 0, 0, 0, 0, 0;
     100, 100, 100, 100, 100, 100;
     200, 200, 220, 240, 270, 310;
     300, 300, 330, 400, 490, 600;
     400, 400, 440, 580, 710, 880;
     500, 500, 560, 770, 950, Lmx;
     600, 600, 680, Lmx, Lmx, Lmx;
     700, 700, 800, Lmx, Lmx, Lmx];
%%%%%%%%%%%%%%%%%%%
% Interpolate for specified Speed and Line rate out (allow for extrapolation)
Vrel = CycD.Speed - CycD.Flow; % SWIMS_relative flow
Lout = Deps * NaN;
for i=1:length(Deps)
    a = interp1(Spds, LD_SbyD100(i,:), Vrel, 'linear','extrap');
    b = interp1(Spds, LD_SbyD120(i,:), Vrel, 'linear','extrap');
    Lout(i) = interp1(Lrates, [a b], CycD.Rates(1), 'linear','extrap');
end
%keyboard
% used specified depths, check for depth where cable runs out
x = interp1(Deps, Lout, CycD.Depths, 'linear','extrap');
ix = find(diff(x)<=0);
for i=1:length(ix)
    x(ix(i)+1) = x(ix(i)+1)+i; % fudge to allow interp at max cable
end
if max(x)>CycD.Line_max % last non-NaN is max depth for available cable
    y = interp1(x, CycD.Depths, CycD.Line_max);
    ix = find(x >= CycD.Line_max);
    x(ix(1)) = CycD.Line_max;
    CycD.Depths(ix(1)) = y;
    if length(ix)>1
        x(ix(2:end)) = NaN;
        CycD.Depths(ix(2:end)) = NaN;
    end
end
CycD.Line_out = x;
% figure cycle times using line_out and line rates:
CycD.Down_time = 60 * (CycD.Line_out/CycD.Rates(1));
CycD.Cycle_time = CycD.Down_time + 60*(CycD.Line_out/CycD.Rates(2));
% m/s = knots/2 for distance
CycD.Cycle_distance = CycD.Cycle_time * (CycD.Speed/2);

CycD.Comment = [];
if Vrel<0 | Vrel>6.5
    CycD.Comment = ['Relative flow=' num2str(Vrel) ' yields uncertain results; '];
end
if CycD.Rates(1)<90 | CycD.Rates(1)>130
    CycD.Comment = [CycD.Comment ...
            'LineRate out=' num2str(CycD.Rates(1)) ' yields uncertain results; '];
end

%%%%%%%%%%%%%%% That's All %%%%%%%%%%%%%%%%%%%%
    
    