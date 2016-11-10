% Plot_SWIMS_Cycles.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
agumisccolumn(5,7)
wysiwyg

ht = subplot(2,1,1);
hd = subplot(2,1,2);
% plot for variety of Speed(kts), Line Rates(m/min);
%Lr = [100, 70];
Lr = [120, 70];
Spds = [2, 3, 4, 5, 6];
ccs = 'rckgbkkkkkkkk';
for i=1:length(Spds)
    cc=ccs(i);
    sp = Spds(i);
    XX = SWIMS_Cycle_Times(sp, Lr);
    axes(ht);
    plot(XX.Depths,XX.Cycle_time/60,['-' cc 'o']);
    hold on
    axes(hd);
    plot(XX.Depths,XX.Cycle_distance,['-' cc 'o']);
    hold on
    clear XX
end

axes(ht);
lgdc = cellstr(num2str(Spds'));
grid on
xlabel('SWIMS depth / m');
ylabel('Total cycle time / min');
title({'Time,Dist to cycle SWIMS to various depths at given';
    ['ship speeds(kts); line rates= ' int2str(Lr(1)) ',' int2str(Lr(2)) ' (m/min)']})
legend(lgdc, 'location','SouthEast');
axes(hd);
grid on
xlabel('SWIMS depth / m');
ylabel('Total cycle distance / m');
legend(lgdc, 'location','NorthEast');

