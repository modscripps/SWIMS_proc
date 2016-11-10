% Table_SWIMS_Cycles.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot for variety of Speed(kts), Line Rates(m/min);
Lr = [100, 70];
%Lr = [120, 70];
Spds = [2, 3, 4, 5, 6];
Deps = [0:100:700]';
clear TB ByT ByX
for i=1:length(Spds)
    sp = Spds(i);
    XX = SWIMS_Cycle_Times(sp, Lr, Deps);
    TB(i) = XX;
    ix = find(~isnan(XX.Depths)); ix = ix(end);
    ByT(1:length(XX.Depths),i) = XX.Cycle_time/60;
    ByX(1:length(XX.Depths),i) = XX.Cycle_distance/1000;
    % save max depth in place of 0-meter stats
    ByT(1,i) = XX.Depths(ix); ByX(1,i) = XX.Depths(ix);
%     plot(XX.Depths,XX.Cycle_time/60,['-' cc 'o']);
%     plot(XX.Depths,XX.Cycle_distance,['-' cc 'o']);
    clear XX
end

clear Tout

Tout{1} = ['SWIMS cycle times/minute, line rates=(' ...
    num2str(Lr(1)) ',' num2str(Lr(2)) ') m/min'];
Tout{end+1} = ['depth  ' sprintf('      %dkt', Spds)];
Tout{end+1} = '----------------------------------------------------------';
Tout{end+1} = [' max/m:' sprintf('%9.1f',ByT(1,:))];
Tout{end+1} = '----------------------------------------------------------';
for i=length(Deps):-1:2
    Tout{end+1} = [sprintf(' %3d   ',Deps(i)) ...
        sprintf('%9.1f',ByT(i,:))];
end

Tout{end+1} = ' ';
Tout{end+1} = '==========================================================';
Tout{end+1} = ' ';
Tout{end+1} = ['SWIMS cycle distances/km, line rates=(' ...
    num2str(Lr(1)) ',' num2str(Lr(2)) ') m/min'];
Tout{end+1} = ['depth  ' sprintf('      %dkt', Spds)];
Tout{end+1} = '----------------------------------------------------------';
Tout{end+1} = [' max/m:' sprintf('%9.1f',ByT(1,:))];
Tout{end+1} = '----------------------------------------------------------';
for i=length(Deps):-1:2
    Tout{end+1} = [sprintf(' %3d   ',Deps(i)) ...
        sprintf('%9.2f',ByX(i,:))];
end
Tout{end+1} = ' ';

% Okay, print it out
for i=1:length(Tout)
    fprintf('%s\n',Tout{i})
end