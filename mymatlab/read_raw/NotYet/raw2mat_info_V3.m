function raw2mat_info_V3(transname,matbasename,inrule,outrule,inrange);
%
% RAW2MAT_INFO_V3 translates raw transect files and outputs a series of matlab 
% 
% USAGE: raw2mat_info_V3(transname,matbasename,inrule,outrule,inrange);
% 
%  This program translates a set of transect files.  It only requires 
%  the prefix from the raw files to run; it incrementally will
%  read all the files it finds in the same directory.  i.e. if you
%  have a deployment MC97001.000, MC97001.001, MC97001.002,...
%  then a call like raw2mat('MC97001',MC97001_') will create
%  a series of files MC97001_1.mat MC97001_2.mat, ....  The numbering
%  of the output files does not correspond to the numbering of the
%  TRANSECT files.  A new matlab file is written every 900 ensembles.
%  The 'info' version EXCLUDES profile data (vels, echo, corr, status, per good)
%

%  Copyright: Jody Klymak, Nov 19. 1998: jklymak@apl.washington.edu
%  Conditions:  This routine is supplied with no warranty.  It is freely
%  distributable under the condition that this copyright notice is kept
%  intact.  V3 by Dave Winkel checks for valid ensembles, skips deviants.


%% Check optional arguments (new: formats(in,out) for multi-filenames, range of seqnos)
if nargin<3
   inrule = {3}; % default 3-digit zero padded file counter;
   % use [] or {0} for single file conversion
end
if nargin<4
   outrule = {1; '.mat'}; % default non-padded out-file counter
end
if nargin<5 | isempty(inrange)
   inrange = [0 inf]; % default, all in-file indices
end
if length(inrange)==1
   inrange = [1 1]*inrange;
end

% Build format arguments for input file names
infmt = ['%s']; inflds = ['transname']; inmult=0;
if iscell(inrule)
   nfi = length(inrule); 
   for i=1:nfi
      x = inrule{i};
      if isnumeric(x) & x>0 % zero-pad for counter
         infmt = [infmt '%0' num2str(x) 'd'];
         inflds = [inflds ',file_no'];
         inmult = inmult+1;
      elseif ischar(x)
         infmt = [infmt '%s'];
         inflds = [inflds ',''' x ''''];
      end
   end
end
if inmult>1
   disp('WARNING: raw2mat_V2, multiple fields for in-file counter');
end
%keyboard

%% Check to see how many input files we have:
nfiles=0; in_files=[];
if inmult<1
   inrange=[1 1]; % single in-file
end
for file_no=inrange(1):inrange(2)
   eval(['inname = sprintf(infmt,' inflds ');']);
   fin=fopen(inname,'r');
   if fin<0
      break;
   end
   fclose(fin);
   nfiles = nfiles+1;
   in_files{nfiles} = inname;
end

if nfiles~=1
   fprintf(1,'Number of in-files found: %d\n',nfiles);
end
if nfiles<1
   return
end

% Build format arguments for output file names
outfmt = ['%s']; outflds = ['matbasename']; outmult=0;
if isempty(outrule)
   outrule={'.mat'};
end
nfo = length(outrule); 
for i=1:nfo
   x = outrule{i};
   if isnumeric(x)
      outfmt = [outfmt '%0' num2str(x) 'd'];
      outflds = [outflds ',matfile'];
      outmult = outmult+1;
   else
      outfmt = [outfmt '%s'];
      outflds = [outflds ',''' x ''''];
   end
end
if outmult>1
   disp('WARNING: raw2mat_V2, multiple fields for out-file counter');
end
%keyboard

% set some parameters now 
MAXENS =900;
if outmult<1
   MAXENS = 5400; % 3000; % hope they all fit (FIX later)
end

% set structures with variable names, their type, and their offset
% within their data block (in bytes)
% I do this so that we have some flexibility should RDI change
% their formats.  It also saves retyping all the names over and over
% during the read and then the write phases of the translation.
name ={'firmwareversion', 'char', 3,
    'firmwarerevison', 'char', 4,
    'sysconfig', 'short', 5,
    'numberbeams','char',9,
    'nbins','char',10,
    'npings','ushort',11,
    'binlen','short',13,
    'blanklen','short',15,
    'watermode','char',17,
    'lowcorr','char',18,
    'codereps','char',19,
    'goodthresh','char',20,
    'errthresh','short',21,
    'tpmin','char',23,
    'tpsec','char',24,
    'tphun','char',25,
    'coords','char',26,
    'headoffset','short',27,
    'headbias','short',29,
    'sensors','char',31,
    'sensorson','char',32,
    'dis1','short',33,
    'pulselen','short',35,
    'fishthresh','char',39,
    'pulselag','short',41
};
 offset=name(:,3);type=name(:,2);name=name(:,1);
 fixleader=struct('name',name,'offset',offset,'type',type);
 
 name={
 'ensemble_number', 'ushort', 3,
    'year','char',5,
    'month','char',6,
    'day','char',7,
    'hour','char',8,
    'minute','char',9,
    'second','char',10,
    'hundreths','char',11,
    'ensMSB','char',12,
    'soundspeed','short',15,
    'xducerdepth','short',17,
    'heading','ushort',19,
    'pitch', 'short' ,21,
    'roll', 'short',23,
    'mptmin','char',29,
    'mptsec','char',30,
    'mpthun','char',31,
    'stdhed','char',32,
    'stdpitch','char',33,
    'stdroll','char',34
    };

 offset=name(:,3);type=name(:,2);name=name(:,1);
 varleader=struct('name',name,'offset',offset,'type',type);
 
 name={
 'btnpings','short',3,
  'btdelay','short',5,
  'btcorrthresh','char',7,
  'btminamp','char',8,
  'btpgmin','char',9,
  'btmode','char',10,
  'btmaxerr','short',11,
  'btrange1','ushort',17,
  'btrange2','ushort',19,
  'btrange3','ushort',21,
  'btrange4','ushort',23,
  'btvel1','short',25,
  'btvel2','short',27,
  'btvel3','short',29,
  'btvel4','short',31,
  'btcor1','uchar',33,
  'btcor2','uchar',34,
  'btcor3','uchar',35,
  'btcor4','uchar',36,
  'btevalamp1','char',37,
  'btevalamp2','char',38,
  'btevalamp3','char',39,
  'btevalamp4','char',40,
  'btpg1','char',41,
  'btpg2','char',42,
  'btpg3','char',43,
  'btpg4','char',44,
  'btmaxdep','short',71,
  'rcvstr1','char',73,
  'rcvstr2','char',74,
  'rcvstr3','char',75,
  'rcvstr4','char',76
};
 offset=name(:,3);type=name(:,2);name=name(:,1);
 btleader=struct('name',name,'offset',offset,'type',type);
 
 name={
 'navutcday','uchar',3,
  'navutcmonth','uchar',4,
  'navutcyear','ushort',5,
  'navutctime1','uint',7,
  'navpctimeoff','int',11,
  'navlat1','int',15,
  'navlon1','int',19,
  'navutctime2','uint',23,
  'navlat2','int',27,
  'navlon2','int',31,
  'navavgspd','short',35,
  'navavgtrackT','short',37,
  'navavgtrackM','short',39,
  'navSMG','short',41,
  'navDMG','ushort',43,
  'navflags','ushort',47,
  'navEnsno','uint',51,
  'navadcpyear','ushort',55,
  'navadcpday','uchar',57,
  'navadcpmonth','uchar',58,
  'navadcptime1','uint',59
};  %more later
 offset=name(:,3);type=name(:,2);name=name(:,1);
 navdata=struct('name',name,'offset',offset,'type',type);
 
 % set a bunch of global variables we will need
 trimstr=[];savestr=[];
 for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);  % make all the fixleader variables global
    % make a string that will let us save these variables whenever we write a 
    % matlab file.
    savestr = [savestr ' ' fixleader(i).name];
 end;
 for i=1:length(varleader)
    eval(['global ',varleader(i).name]);
    % make a string that lets us make these variables the right size if
    % we write before the buffer is full...
    trimstr = [trimstr varleader(i).name '=' varleader(i).name '(1:ensembles);'];
    savestr = [savestr ' ' varleader(i).name ];
 end;
 % Initialize optional Bottom track data and VmDas-(.STA,.LTA) Navg data
 for i=1:length(btleader)
    eval(['global ' btleader(i).name]);
    eval([btleader(i).name ' = [];']);
    savestr = [savestr ' ' btleader(i).name];
    trimstr = [trimstr btleader(i).name '=' btleader(i).name ...
          '(1:min(ensembles,length(' btleader(i).name ' )));'];
 end
 for i=1:length(navdata)
    eval(['global ' navdata(i).name]);
    eval([navdata(i).name ' = [];']);
    savestr = [savestr ' ' navdata(i).name];
    trimstr = [trimstr navdata(i).name '=' navdata(i).name ...
          '(1:min(ensembles,length(' navdata(i).name ' )));'];
 end
 % Just count profiles detected (no data saved)
 savestr = [savestr ' velread corread echoread percread statusread'];
 
%% Ready to convert ...
matfile=1;
ensembles=1;
eval(['matname = sprintf(outfmt,' outflds ');']); % first matlab file-name
% now loop through all the files
for file_no=1:nfiles
   % open next in-file
   %fprintf(1,'Translating %s\n',in_files{file_no});
   fin=fopen(in_files{file_no},'r','ieee-le'); % open as a little endian
   % Revised to find Ensemble starts=7F7F, where to stop for multi-Ens files
   A = fread(fin, inf, 'uchar');
   %keyboard
   iEpos = find(A(1:end-1)==127 & A(2:end)==127); % potential Ens starts
   if isempty(iEpos), continue; end % no valid ensembles in this file
   if iEpos(end)+10 > length(A), iEpos(end)=[]; end
   Enbyt = A(iEpos+3)*256+A(iEpos+2); % bytes in ensemble (w/o checksum)
   ChkAdd=NaN*Enbyt; ChkSum=ChkAdd;
   for iE=1:length(iEpos)
       lb=iEpos(iE)+Enbyt(iE)-1;
       if lb>iEpos(iE) & lb+2<=length(A)
           ChkAdd(iE) = mod(sum(A(iEpos(iE):lb)), 65536); % computed
           ChkSum(iE) = A(lb+2)*256+A(lb+1); % recorded
       end
   end
   iZ=find(ChkAdd-ChkSum==0); % exclude if computed~=recorded
   keyboard
   if isempty(iZ), continue; end % no valid ensembles in this file
   iEpos = iEpos(iZ); Enbyt = Enbyt(iZ);
   % Final check, ALL ensembles should be same length
   EnLEN = median(Enbyt);
   ix = find(Enbyt ~= EnLEN);
   if ~isempty(ix)
       iEpos(ix) = []; Enbyt(ix) = [];
       disp(['   Excluding ' num2str(length(ix)) ' bogus ensemble(s).'])
   end
   iEpos = iEpos-1; % positions just before good ensembles
   %keyboard
   Hdrbytes=[]; % Save for diagnosis
   velread = 0;
   corread = 0;
   echoread = 0;
   percread = 0;
   statusread = 0;

   for iEns=1:length(iEpos)
      file_pos = iEpos(iEns);
      fseek(fin,file_pos,'bof');
      A=local_hdread(fin);
      ndata=A(1);nbytes=A(2);offsets=A(3:2+ndata);
      Hdrbytes=[Hdrbytes ndata nbytes offsets];
      %if iEns>16 & iEns<22, keyboard; end
      % the raw data is variable length, with variable pieces
      % for storing the data so we need to be flexible in 
      % translating the files.  ndata tells us how many data
      % blocks we have.  offsets tells us where each of these
      % data blocks is in the file stream...
      for data_type=1:ndata
         fseek(fin,file_pos+offsets(data_type),'bof');
         data_id=fread(fin,1,'short');
         switch data_id
         case 0  % hex bytes = 00 00
            local_fixread(fixleader,fin,file_pos+offsets(data_type));
         case 128 % 80 00
            local_varread(varleader,fin,file_pos+offsets(data_type),...
                ensembles);
        case 256 % 00 01
            velread = velread + 1;
        case 512 % 00 02
            corread= corread + 1;
         case 768 % 00 03
            echoread = echoread + 1;
         case 1024 % 00 04
            percread = percread + 1;
         case 1280 % 00 05
            statusread = statusread + 1;
         case 1536 % 00 06
            local_btread(btleader,fin,file_pos+offsets(data_type),...
               ensembles);
         case 8192 % 00 20
            local_navread(navdata,fin,file_pos+offsets(data_type),...
               ensembles);
         end; % switch for conditional data reads....
      end; % going through the data types
      
      % Check to see if ensembles is high enough to write a
      % matlab file yet...
      if ensembles==MAXENS
         eval(trimstr);
         eval(['save ' matname savestr]);
         matfile=matfile+1;
         eval(['matname = sprintf(outfmt,' outflds ');']); % next matlab file-name
         ensembles=0;
      end;
      ensembles=ensembles+1;
      % done with this ensemble, on to next
end; % while going through headers in the file
	%MHA 1/1/99 change: close the file.  If not fid gets updated and it can't do
	%more than 512.
	fclose(fin);
    %keyboard
end; % for file_no=nfiles
% write the last mat file...
if ensembles>1 % was 0, DPW 3-2001
   ensembles=ensembles-1;
   eval(trimstr);
   eval(['save ' matname savestr]);
   matfile=matfile+1;
   ensembles=0;
else
    % No ensembles found, just save empty 'ensemble_number'
    ensemble_number=[];
    save(matname, 'ensemble_number');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%% LOCALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ndata]=local_hdread(fin);

offset=ones(1,12)*NaN;
% reads the header information

id = fread(fin,1,'uchar');
if (id~=127)
   ndata=[];
   return;
end;

src= fread(fin,1,'uchar');
nbytes=fread(fin,1,'ushort');
junk=fread(fin,1,'uchar');
ndata=fread(fin,1,'uchar');
for i=1:ndata
   offset(i)=fread(fin,1,'ushort');
end;
ndata = [ndata nbytes offset(1:ndata)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_readfix %%%%%%%%%%%%%%%%%%%%%%
function local_fixread(fixleader,fin,data_pos);
% fixleader is set above.  It allows us to define a bunch of variables
% (stored in the name portion of the structure)
for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);
 end;
 for i=1:length(fixleader)
    fseek(fin,data_pos+fixleader(i).offset-1,'bof');
    
    eval([fixleader(i).name '=fread(fin,1,''' fixleader(i).type ''');']);
    % i.e. firmwareversion=fread(fin,1,'char');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_readvar %%%%%%%%%%%%%%%%%%%%%%
function local_varread(fixleader,fin,data_pos,ensemble);
for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);
end;
for i=1:length(fixleader)
   fseek(fin,data_pos+fixleader(i).offset-1,'bof');
   eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
   % i.e. ensemble_number(ensemble)=fread(fin,1,'short');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_btread %%%%%%%%%%%%%%%%%%%%%%
function local_btread(fixleader,fin,data_pos,ensemble);
% same tricks are used here as for local_varread...
for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);
end;
for i=1:length(fixleader)
   fseek(fin,data_pos+fixleader(i).offset-1,'bof');
   eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_navread %%%%%%%%%%%%%%%%%%%%%%
function local_navread(fixleader,fin,data_pos,ensemble);
% same tricks are used here as for local_varread...
for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);
end;
for i=1:length(fixleader)
   fseek(fin,data_pos+fixleader(i).offset-1,'bof');
   eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
end;


