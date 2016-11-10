function [SWIMSctd]=Read_raw2_SWIMS(filename, Version)%Read all data from a SWIMS data file.  Variables returned are...%%SWIMS 2 CTD files have:% first 4 bytes are integer clock value - in seconds% next 8 bytes are double elapsed time value - % next 30 (was:24) bytes are a CTD data scan%if nargin<2 | isempty(Version)    Version = 0; % for SWIMS 2end%Open the file as little-endian since it was written on a PC.[fid,message]=fopen(filename,'r','ieee-le');% read the clock values, skip 32 bytes after each[SWIMSctd.SWIMStime, count] = fread(fid, inf, 'uint32', 38);% read elapsed time values, skip 28 bytes after eachstatus = fseek(fid, 4, 'bof');[SWIMSctd.SWIMSfasttime, count] = fread(fid, inf, 'double', 34);%These both produce the time in seconds.%plot(SWIMStime-SWIMStime(1),SWIMSfasttime-SWIMSfasttime(1),SWIMStime-SWIMStime(1),SWIMStime-SWIMStime(1),'k--')frewind(fid);[d, count] = fread(fid, [42, inf], 'uchar');%Get temp - first three bytes after the time stampsSWIMSctd.tfreq=d(12+1,:)*256+d(12+2,:)+d(12+3,:)/256;%Get cond - next threeSWIMSctd.cfreq=d(12+4,:)*256+d(12+5,:)+d(12+6,:)/256;%Get pres - next threeSWIMSctd.pfreq=d(12+7,:)*256+d(12+8,:)+d(12+9,:)/256;%Get temp2 - first three bytes after the time stampsSWIMSctd.tfreq2=d(12+10,:)*256+d(12+11,:)+d(12+12,:)/256;%Get cond2 - next threeSWIMSctd.cfreq2=d(12+13,:)*256+d(12+14,:)+d(12+15,:)/256;%Get pressure-temperature comp byte%td=(d(12+22,:)/5)-10;%A/D Data.%Now rewind and load in data as 12-bit words.%The 42-byte scan (4 bytes + 8 bytes + 30 bytes of CTD data) %is then 28 12-bit words.  The first 8 are the time stamps.%Then the next 6 12-bit words are the first 9 bytes of the CTD scan%(the P, T and C data).  Then 4 are the six bytes of T2 and C2.  %Then the next 8 (19-26) are the a/d channels.  The last %three bytes make two more 12-bit words, for a total of 28.%Note it is important to close and reopen the file as bigendian so the %bit orders are correct.fclose(fid);[fid,message]=fopen(filename,'r','ieee-be');frewind(fid);%[a,count]=fread(fid,inf,'uchar');[d2, count] = fread(fid, [28, inf], 'ubit12');%close the filestatus = fclose(fid);addata=d2(19:26,:);%Convert it to volts (0 to 5 for SWIMs2)%addata=10*(1 - 2*addata/4096);SWIMSctd.addata=5*(1 - addata/4095);td=d2(27,:); % next 12-bits after A-to-DAD590M = .01286587327; AD590B = -10.05189772; % from SBE9+ config sheetSWIMSctd.td = AD590M*td + AD590B; % pressure-temp comp (degC?)SWIMSctd.status=d2(28,:)/256; % first 4 bitsSWIMSctd.modCT = mod(d2(28,:),256); % last 8 bits%