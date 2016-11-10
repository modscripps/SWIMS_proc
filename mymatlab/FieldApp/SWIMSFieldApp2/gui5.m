function varargout = GUI5(varargin)
% GUI5 Application M-file for GUI5.fig
%    FIG = GUI5 launch GUI5 GUI.
%    GUI5('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 01-May-2002 14:51:12

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    
    
    %Initialize everything here
    disp 'Initializing...'
    %SetParams;
    handles.SW=SetInitialParamsGUI;
    set(handles.edit1,'String',handles.SW.printfilename)
    set(handles.edit2,'String',num2str(handles.SW.beg_time))
    set(handles.edit3,'String',num2str(handles.SW.end_time))
    set(handles.edit4,'String',num2str(handles.SW.zmin))
    set(handles.edit5,'String',num2str(handles.SW.zmax))
    set(handles.edit6,'String',num2str(handles.SW.dmin1))
    set(handles.edit7,'String',num2str(handles.SW.dmax1))
    set(handles.edit8,'String',num2str(handles.SW.dmin2))
    set(handles.edit9,'String',num2str(handles.SW.dmax2))
    set(handles.edit10,'String',num2str(handles.SW.dmin3))
    set(handles.edit11,'String',num2str(handles.SW.dmax3))
    set(handles.edit12,'String',num2str(handles.SW.lonmin))
    set(handles.edit13,'String',num2str(handles.SW.lonmax))
    set(handles.edit14,'String',num2str(handles.SW.latmin))
    set(handles.edit15,'String',num2str(handles.SW.latmax))
    set(handles.edit16,'String',handles.SW.remotebasepath)
    set(handles.edit19,'String',handles.SW.localbasepath)
    set(handles.edit20,'String',handles.SW.prpath)
    set(handles.edit21,'String',num2str(handles.SW.printfilenum))
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
disp 'Copying files...'
err=GetNewFiles(handles.SW);
disp 'Done.'


% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
%Load in the data
disp 'loading...'
handles.SW.SWIMSgrid = get_swims_data(handles.SW.beg_timereq, handles.SW.end_timereq, fullfile(handles.SW.localindexdir,handles.SW.index_fileSWIMS), handles.SW.data_pathSWIMS,[],1);
%disp('adcp.')
handles.SW.ADCP = get_adcp_data(handles.SW.beg_timereq, handles.SW.end_timereq, fullfile(handles.SW.localindexdir,handles.SW.index_fileADCP), handles.SW.data_pathADCP,[],1);

%Set beg time
if handles.SW.beg_timereq < min(handles.SW.SWIMSgrid.yday)
    handles.SW.beg_time=min(handles.SW.SWIMSgrid.yday);
end

%Now set end_time
if handles.SW.end_timereq > max(handles.SW.SWIMSgrid.yday)
    handles.SW.end_time=max(handles.SW.SWIMSgrid.yday);
end

%Now update the GUIs
set(handles.edit2,'String',num2str(handles.SW.beg_time))
set(handles.edit3,'String',num2str(handles.SW.end_time))

disp 'done.'
guidata(h,handles);



% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)
%disp 'Plot; OK'
disp 'Plotting...'
PlotFigureGUI(handles.SW)
disp 'Done.'
% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
%disp 'Print; OK'
handles.SW=PrintFigureGUI(handles.SW);
set(handles.edit1,'String',handles.SW.printfilename)
set(handles.edit21,'String',handles.SW.printfilenum)
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
%Print filename edit box.
%If it is changed, reset the counter to 1.
tmp=get(h,'String');
if strcmp(tmp,handles.SW.printfilename)~=1
handles.SW.printfilename=tmp;
handles.SW.printfilenum=1;
set(handles.edit21,'String',num2str(handles.SW.printfilenum))
end
guidata(h,handles)

% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
%tmin
str=get(h,'String');
handles.SW.beg_time=str2num(str);
handles.SW.beg_timereq=handles.SW.beg_time;
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
%tmax
str=get(h,'String');
handles.SW.end_time=str2num(str);
%handles.SW.end_timereq=handles.SW.end_time;
guidata(h,handles);




% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
%zmin
str=get(h,'String');
handles.SW.zmin=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)
%zmax
str=get(h,'String');
handles.SW.zmax=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)
%
str=get(h,'String');
handles.SW.dmin1=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit7_Callback(h, eventdata, handles, varargin)

str=get(h,'String');
handles.SW.dmax1=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit8_Callback(h, eventdata, handles, varargin)

str=get(h,'String');
handles.SW.dmin2=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit9_Callback(h, eventdata, handles, varargin)

str=get(h,'String');
handles.SW.dmax2=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit10_Callback(h, eventdata, handles, varargin)

str=get(h,'String');
handles.SW.dmin3=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = edit11_Callback(h, eventdata, handles, varargin)

str=get(h,'String');
handles.SW.dmax3=str2num(str);

guidata(h,handles);



% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)
val = get(h,'Value');
string_list = get(h,'String');
handles.SW.varstr1=string_list{val};
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = popupmenu2_Callback(h, eventdata, handles, varargin)

val = get(h,'Value');
string_list = get(h,'String');
handles.SW.varstr2=string_list{val};
guidata(h,handles);



% --------------------------------------------------------------------
function varargout = popupmenu3_Callback(h, eventdata, handles, varargin)

val = get(h,'Value');
string_list = get(h,'String');
handles.SW.varstr3=string_list{val};

guidata(h,handles);


% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)
handles.SW.sm=get(h,'Value');

guidata(h,handles);


% --------------------------------------------------------------------
function varargout = edit12_Callback(h, eventdata, handles, varargin)
%lonmin
str=get(h,'String');
handles.SW.lonmin=str2num(str);

guidata(h,handles);




% --------------------------------------------------------------------
function varargout = edit13_Callback(h, eventdata, handles, varargin)
%lonmax
str=get(h,'String');
handles.SW.lonmax=str2num(str);

guidata(h,handles);




% --------------------------------------------------------------------
function varargout = edit14_Callback(h, eventdata, handles, varargin)
%latmin 
str=get(h,'String');
handles.SW.latmin=str2num(str);

guidata(h,handles);




% --------------------------------------------------------------------
function varargout = edit15_Callback(h, eventdata, handles, varargin)
%latmax
str=get(h,'String');
handles.SW.latmax=str2num(str);

guidata(h,handles);




% --------------------------------------------------------------------
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
%Return SW in a global variable that the workspace can access.
global SWglob;

SWglob=handles.SW;

disp 'Data are now available as SWglob.SWIMSgrid and SWglob.ADCP.  Type ''global SWglob'' to access these.'


% --------------------------------------------------------------------
function varargout = edit16_Callback(h, eventdata, handles, varargin)
%remote path
str=get(h,'String');
handles.SW.remotebasepath=str;
%remote paths
handles.SW.remotedata_pathSWIMS=[str 'griddata\'];
handles.SW.remotedata_pathADCP=[str 'data_mat\ADCP\'];
handles.SW.remoteindexdir=[str 'indexes\'];


guidata(h,handles);




% --------------------------------------------------------------------
function varargout = edit19_Callback(h, eventdata, handles, varargin)
str=get(h,'String');
handles.SW.localbasepath=str;
%local paths
handles.SW.data_pathSWIMS=[str 'griddata'];
handles.SW.data_pathADCP=[str 'data_mat\ADCP\'];
handles.SW.localindexdir=[str 'indexes\'];

guidata(h,handles);





% --------------------------------------------------------------------
function varargout = edit20_Callback(h, eventdata, handles, varargin)

str=get(h,'String');
handles.SW.prpath=str;

guidata(h,handles);




% --------------------------------------------------------------------
function varargout = edit21_Callback(h, eventdata, handles, varargin)
%The number of the print file
str=get(h,'String');
handles.SW.printfilenum=str2num(str);

guidata(h,handles);
