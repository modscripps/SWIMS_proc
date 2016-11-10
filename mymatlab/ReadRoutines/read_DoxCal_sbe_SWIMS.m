function [DoxCoefs]=read_DoxCal_sbe_SWIMS(sn,year,yday)%function [DoxCoefs]=read_DoxCal_sbe_SWIMS(sn,year,yday)%Return the dissolved oxygen cal coefficients for a given S/N, year and yearday.%The function searches the global directory stored in SWIMS_cal/dox for%a file whose name is the desired S/N, which contains a list of dates%and coefficients for each calibration.  It reads the last ones before %the desired date and returns them in structure DoxCoefs.%  apr,2002 - dpw (for SBE 43 model)% jan 2012 (dpw) - modified for additional SB coeffs by using cell arrayglobal SWIMS_calDoxCoefs.Soc=NaN;dtnum_for = datenum(year,1,1) + yday;%open this file.filename=fullfile(SWIMS_cal, 'dox' , sn);fid=fopen(filename,'r');if fid == -1	disp(['Error: no cal file exists for DOX S/N=' sn])	return;end%Read the first three header lines.line=fgetl(fid); line=fgetl(fid); line=fgetl(fid);done=0; cfs = [];count=1;line=fgetl(fid);while isstr(line)    cfs{count} = [];    dstr = [line(6:7) '/' line(9:10) '/' line(1:4) ' ' ...            line(12:13) ':' line(15:16) ':' line(18:19)]; % mm/dd/yyyy hh:mm:ss    dtnum_cal(count) = datenum(dstr);    [dstr,x] = strtok(line);    % Save coeff sets in rows of cfs(cal#,:)     for i=1:100        [cstr,x] = strtok(x);        cfs{count} = [cfs{count} str2num(cstr)];        % all entries in file must have same number of coeffs        if isempty(deblank(x))            break        end    end	count=count+1;		line=fgetl(fid);endfclose(fid);%Find the last calibration before the desired date.ind = find(dtnum_cal<=dtnum_for);if isempty(ind)    disp(['No calib for S/N=' sn ' in effect for ' datestr(dtnum_for)])    returnend[x, im] = max(dtnum_cal(ind)); % may not be monotonic (esp. for testing)ind=ind(im); cfs = cfs{ind};if length(cfs)==5 % older SB formula    DoxCoefs.Soc=cfs(1); DoxCoefs.Boc=cfs(2);    DoxCoefs.Voffset=cfs(3);     DoxCoefs.TCor=cfs(4); DoxCoefs.PCor=cfs(5);elseif length(cfs)==8 % newer one, with time constant Tau(P,T)    DoxCoefs.Soc=cfs(1); DoxCoefs.Tau20=cfs(2);    DoxCoefs.Voffset=cfs(3);     DoxCoefs.A=cfs(4); DoxCoefs.B=cfs(5); DoxCoefs.C=cfs(6);    DoxCoefs.D1=cfs(7); DoxCoefs.D2=cfs(8);else    disp(['Unknown calib format, Dox S/N=' sn ' for ' datestr(dtnum_for)])end