% AD_SWIMS_Proc_Fix_wanliwapr2013.m - Assign pitch,roll values to ADUP
%  using ADDN data (ADUP pitch jumps +4 deg sometimes)

DNidx = ['E:\swims\' cruise '\indexes\ADDN_' cruise '_matfiles.mat'];
DNfld = ['E:\swims\' cruise '\data_mat\ADDN'];

EnsDNtoUP = [113.9,0];
EnsYDgap = [113.90,114.25,114.32,inf]; % ens no's reset after these ydays
ADFix = get_ADupdn_data(yd_b-ex_yday,yd_e+ex_yday, DNidx, DNfld, 1);

% remove unneeded array fields, offset ADDN ens_no's to align w/ADUP's
ADFix = rmfield(ADFix,{'v1_bm','ec1_bm','cor1_bm', 'v2_bm','ec2_bm','cor2_bm', ...
    'v3_bm','ec3_bm','cor3_bm', 'v4_bm','ec4_bm','cor4_bm'});
ix = find( EnsDNtoUP(:,2) ); % non-zero offsets
ydl = [ EnsDNtoUP(:,1); inf ]; % yearday limits
% Do NOT alter ADUP ens_no's at this processing stage
for i=1:length(ix)
    ie = find( ADFix.yday>=ydl(ix(i)) & ADFix.yday<ydl(ix(i)+1) );
    ADFix.ens_no(ie) = ADFix.ens_no(ie) + EnsDNtoUP(ix(i),2);
end

warning off MATLAB:interp1:NaNinY
for iGn=1:length(EnsYDgap)-1
    irU = find(ADRaw.yday>=EnsYDgap(iGn) & ADRaw.yday<EnsYDgap(iGn+1));
    irD = find(ADFix.yday>=EnsYDgap(iGn) & ADFix.yday<EnsYDgap(iGn+1));
    if length(irU)<2 || length(irD)<2
        continue % need two to interpolate
    end
    % copy values from ADDN, unless ens_no's out-of-range
    xx = interp1(ADFix.ens_no(irD),ADFix.pitch(irD), ADRaw.ens_no(irU));
    ix = find(~isnan(xx));
    if ~isempty(ix)
        ADRaw.roll(irU(ix)) = - xx(ix); % ADUP.roll = - ADDN.pitch
    end
    xx = interp1(ADFix.ens_no(irD),ADFix.roll(irD), ADRaw.ens_no(irU));
    ix = find(~isnan(xx));
    if ~isempty(ix)
        ADRaw.pitch(irU(ix)) = xx(ix); % ADUP.pitch = ADDN.roll
    end
end

clear ADFix xx ix irU irD ie