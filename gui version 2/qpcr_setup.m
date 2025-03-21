function varargout = qpcr_setup(varargin)
% QPCR_SETUP M-file for qpcr_setup.fig
%      QPCR_SETUP, by itself, creates a new QPCR_SETUP or raises the existing
%      singleton*.
%
%      H = QPCR_SETUP returns the handle to a new QPCR_SETUP or the handle to
%      the existing singleton*.
%
%      QPCR_SETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QPCR_SETUP.M with the given input arguments.
%
%      QPCR_SETUP('Property','Value',...) creates a new QPCR_SETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before qpcr_setup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to qpcr_setup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qpcr_setup

% Last Modified by GUIDE v2.5 25-Nov-2015 21:08:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qpcr_setup_OpeningFcn, ...
                   'gui_OutputFcn',  @qpcr_setup_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before qpcr_setup is made visible.
function qpcr_setup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qpcr_setup (see VARARGIN)

% Choose default command line output for qpcr_setup
handles.output = hObject;

% custom stuff
handles.alluserdat.rawdata=[];

handles.alluserdat.num_grps=0;
handles.alluserdat.grp_names={};
handles.alluserdat.grp_samples={};

handles.alluserdat.housekeeping_idx=1;
handles.alluserdat.secondarynorm_idx=0;

param.log=false;
param.geo=false;
param.legend=true;
param.grps2plot=[];
param.tars2plot=[];
handles.alluserdat.graph_param=param;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes qpcr_setup wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = qpcr_setup_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in check_log.
function check_log_Callback(hObject, eventdata, handles)
% hObject    handle to check_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_log
val=get(handles.check_log,'Value');
handles.alluserdat.graph_param.log=boolean(val);
guidata(hObject,handles);

% --- Executes on button press in push_loadbttn.
function push_loadbttn_Callback(hObject, eventdata, handles)
% hObject    handle to push_loadbttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newses=isempty(handles.alluserdat.rawdata);
if ~newses
    ansq=questdlg('Do you want to replace what you have currently?','Warning','Yes','No','Yes');
else
    ansq='Yes';
end
switch ansq
    case 'Yes'
        temp=dir;
        ind_xls_files=[];
        for n=1:length(temp)
            tempname=temp(n).name;
            len_tn=length(tempname);
            if len_tn>4
                if strcmpi(tempname((len_tn-3):(len_tn)),'.xls') ||  strcmpi(tempname((len_tn-4):(len_tn)),'.xlsx')
                    ind_xls_files=[ind_xls_files,n];
                end
            end
        end
        
        
        if length(ind_xls_files)==1
            dirn=pwd;
            dirn=[pwd,'\'];
            fn=temp(ind_xls_files).name;
        else
            [fn,dirn] = uigetfile('*.xls;*.xlsx','All Excel files');
        end
        if strcmpi(fn((end-3):(end)),'.xls')
            ext='.xls';
        else
            ext='.xlsx';
        end
        filename=[dirn,fn];
        if length(ext)==4
            filename2save=[filename(1:(length(filename)-4))];
        else
            filename2save=[filename(1:(length(filename)-5))];
        end
        
        handles.alluserdat.filename=filename2save;
        handles.alluserdat.fileext=ext;
        
        handles.alluserdat.rawdata=gui_read_QuantStudio6_qpcrdata(filename,filename2save);
        cd(dirn)
        if ~newses
            handles.alluserdat.num_grps=0;
            handles.alluserdat.grp_names={};
            handles.alluserdat.grp_samples={};
            handles.alluserdat.housekeeping_idx=1;
            handles.alluserdat.secondarynorm_idx=0;
            handles.alluserdat.graph_param.grps2plot=[];
            handles.alluserdat.graph_param.tars2plot=[];
            set(handles.popup_sel_housekeeping,'Value',handles.alluserdat.housekeeping_idx)
            set(handles.list_secondary_norm,'Value',handles.alluserdat.secondarynorm_idx+1)
        end
        
        guidata(hObject,handles);
        update_early_guis(hObject,handles);
        if ~newses
            handles=update_numgrp_guis(hObject,handles);
            guidata(hObject,handles);
            update_curgrp_guis(hObject,handles);
            
            update_grp_list(hObject,handles);
            update_graph_params(hObject,handles);
        end
        
    case 'No'
end

% --- Executes on button press in push_gengraph.
function push_gengraph_Callback(hObject, eventdata, handles)
% hObject    handle to push_gengraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
raw=handles.alluserdat.rawdata;
tempidx=get(handles.popup_sel_housekeeping,'Value');
templist=get(handles.popup_sel_housekeeping,'String');
ctrlTar{1}=templist{tempidx};
tempidx=get(handles.popup_selgrp_ctrl,'Value');
grpsamps=handles.alluserdat.grp_samples;
ctrlSamp=tempidx;%grpsamps{tempidx};
grpnames=get(handles.popup_selgrp,'String');


gui_proc_QuantStudio6_qpcrdata(raw,ctrlTar,handles.alluserdat.graph_param.rem_tarlist,ctrlSamp,grpnames,grpsamps,handles.alluserdat.graph_param);


% --- Executes on selection change in list_chosen_samp.
function list_chosen_samp_Callback(hObject, eventdata, handles)
% hObject    handle to list_chosen_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_chosen_samp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_chosen_samp


% --- Executes during object creation, after setting all properties.
function list_chosen_samp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_chosen_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_avail_samp.
function list_avail_samp_Callback(hObject, eventdata, handles)
% hObject    handle to list_avail_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_avail_samp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_avail_samp


% --- Executes during object creation, after setting all properties.
function list_avail_samp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_avail_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_choose.
function push_choose_Callback(hObject, eventdata, handles)
% hObject    handle to push_choose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
grpsamps=handles.alluserdat.grp_samples;
if ~isempty(grpsamps)
sel=get(handles.list_avail_samp,'Value');
availsamps=get(handles.list_avail_samp,'String');

grpsel=get(handles.popup_selgrp,'Value');

oldval=grpsamps{grpsel};
grpsamps{grpsel}=unique([oldval;availsamps(sel)]);
handles.alluserdat.grp_samples=grpsamps;

guidata(hObject,handles)
update_curgrp_guis(hObject,handles)
end

% --- Executes on button press in push_unchoose.
function push_unchoose_Callback(hObject, eventdata, handles)
% hObject    handle to push_unchoose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

grpsamps=handles.alluserdat.grp_samples;
if ~isempty(grpsamps)
grpsel=get(handles.popup_selgrp,'Value');

if ~isempty(grpsamps{grpsel})
    sel=get(handles.list_chosen_samp,'Value');
    availsamps=get(handles.list_chosen_samp,'String');
    
    newlist={};
    ct=1;
    for m=1:length(availsamps)
        if ~any(m==sel)
            newlist{ct,1}=availsamps{m};            
            ct=ct+1;
        end        
    end
    grpsamps{grpsel}=newlist;
    if get(handles.list_chosen_samp,'Value')>length(newlist)
        set(handles.list_chosen_samp,'Value',length(newlist))
    end
    
    handles.alluserdat.grp_samples=grpsamps;
    guidata(hObject,handles)
    update_curgrp_guis(hObject,handles)
end
end

% --- Executes on button press in check_legend.
function check_legend_Callback(hObject, eventdata, handles)
% hObject    handle to check_legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_legend
val=get(handles.check_legend,'Value');
handles.alluserdat.graph_param.legend=boolean(val);
guidata(hObject,handles);

function edit_numgrps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numgrps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numgrps as text
%        str2double(get(hObject,'String')) returns contents of edit_numgrps as a double
numgrps=str2double(get(handles.edit_numgrps,'String'));
if isnan(numgrps)
    badinput=true;
elseif (round(numgrps)-numgrps)~=0
    badinput=true;
else
    badinput=false;
end

if badinput
    numgrps=handles.alluserdat.num_grps;
    set(handles.edit_numgrps, 'String',num2str(numgrps));
end


% --- Executes during object creation, after setting all properties.
function edit_numgrps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numgrps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_updatenumgrps.
function push_updatenumgrps_Callback(hObject, eventdata, handles)
% hObject    handle to push_updatenumgrps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
numgrps=str2double(get(handles.edit_numgrps,'String'));
handles.alluserdat.num_grps=numgrps;

guidata(hObject,handles)
handles=update_numgrp_guis(hObject,handles);
guidata(hObject,handles)
update_curgrp_guis(hObject,handles)

% --- Executes on selection change in popup_selgrp.
function popup_selgrp_Callback(hObject, eventdata, handles)
% hObject    handle to popup_selgrp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_selgrp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_selgrp
set(handles.list_chosen_samp,'Value',1)
update_curgrp_guis(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popup_selgrp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_selgrp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_grpname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_grpname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_grpname as text
%        str2double(get(hObject,'String')) returns contents of edit_grpname as a double
grpnames=handles.alluserdat.grp_names;
grpsel=get(handles.popup_selgrp,'Value');

newname=get(handles.edit_grpname,'String');
grpnames{grpsel}=newname;
handles.alluserdat.grp_names=grpnames;

guidata(hObject,handles)
set(handles.popup_selgrp,'String',grpnames)
set(handles.popup_selgrp_ctrl,'String',grpnames);
update_grp_list(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_grpname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_grpname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popup_sel_housekeeping.
function popup_sel_housekeeping_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sel_housekeeping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_sel_housekeeping contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sel_housekeeping
handles=update_early_guis(hObject,handles);

handles.alluserdat.housekeeping_idx=get(handles.popup_sel_housekeeping,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popup_sel_housekeeping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sel_housekeeping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_secondary_norm.
function list_secondary_norm_Callback(hObject, eventdata, handles)
% hObject    handle to list_secondary_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_secondary_norm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_secondary_norm
handles.alluserdat.secondarynorm_idx=get(handles.list_secondary_norm,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function list_secondary_norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_secondary_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_selgrps2plot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selgrps2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selgrps2plot as text
%        str2double(get(hObject,'String')) returns contents of edit_selgrps2plot as a double
old_grps2plot=handles.alluserdat.graph_param.grps2plot;
grpnames=handles.alluserdat.grp_names;

numgrps=length(grpnames);
new_grps2plot=str2num(get(handles.edit_selgrps2plot,'String'));
if numgrps>0
    if any(isnan(new_grps2plot))
        newstr=num2str(old_grps2plot);
    elseif max(new_grps2plot)>numgrps || min(new_grps2plot)<1
        newstr=num2str(old_grps2plot);
    elseif any((round(new_grps2plot)-new_grps2plot)~=0)
        newstr=num2str(old_grps2plot);
    else
        newstr=num2str(new_grps2plot);
        handles.alluserdat.graph_param.grps2plot=new_grps2plot;
    end
else
    newstr='';
end

set(handles.edit_selgrps2plot,'String',newstr)
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_selgrps2plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selgrps2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_seltar2plot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_seltar2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_seltar2plot as text
%        str2double(get(hObject,'String')) returns contents of edit_seltar2plot as a double
old_tars2plot=handles.alluserdat.graph_param.tars2plot;
tarnames=handles.alluserdat.rawdata.TarNames;

numtars=length(tarnames);
new_tars2plot=str2num(get(handles.edit_seltar2plot,'String'));
if numtars>0
    if any(isnan(new_tars2plot))
        newstr=num2str(old_tars2plot);
    elseif max(new_tars2plot)>numtars || min(new_tars2plot)<1
        newstr=num2str(old_tars2plot);
    elseif any((round(new_tars2plot)-new_tars2plot)~=0)
        newstr=num2str(old_tars2plot);
    else
        newstr=num2str(new_tars2plot);
        handles.alluserdat.graph_param.tars2plot=new_tars2plot;
    end
else
    newstr='';
end

set(handles.edit_seltar2plot,'String',newstr)
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_seltar2plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_seltar2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in check_geometric.
function check_geometric_Callback(hObject, eventdata, handles)
% hObject    handle to check_geometric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_geometric
val=get(handles.check_geometric,'Value');
handles.alluserdat.graph_param.geo=boolean(val);
guidata(hObject,handles);


% --- Executes on selection change in popup_selgrp_ctrl.
function popup_selgrp_ctrl_Callback(hObject, eventdata, handles)
% hObject    handle to popup_selgrp_ctrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_selgrp_ctrl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_selgrp_ctrl


% --- Executes during object creation, after setting all properties.
function popup_selgrp_ctrl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_selgrp_ctrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------- Menus ------------------------

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_save_ses_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_ses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,path]=uiputfile('*.mat', 'Save Session');
if ~isequal(filename,0)
    if ~isempty(handles.alluserdat.rawdata)
        alluserdat=handles.alluserdat;
        save([path,filename],'alluserdat')
    end
end

% --------------------------------------------------------------------
function menu_load_ses_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_ses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,path]=uigetfile('*.mat', 'Save Session');
if ~isequal(filename,0)
    load([path,filename])
    if exist('alluserdat','var')
        handles.alluserdat=alluserdat;
        guidata(hObject,handles)
        
        
        update_early_guis(hObject,handles);
        update_numgrp_guis(hObject,handles);
        update_curgrp_guis(hObject,handles);
        update_grp_list(hObject,handles);
        update_graph_params(hObject,handles);   
        
        set(handles.popup_sel_housekeeping,'Value',handles.alluserdat.housekeeping_idx)
        set(handles.list_secondary_norm,'Value',handles.alluserdat.secondarynorm_idx+1)
        set(handles.edit_numgrps,'String',num2str(length(handles.alluserdat.grp_names)))
        update_early_guis(hObject,handles);
    else
        errordlg('File is not appropriate')
    end
end




% --------------------------------------------------------------------
function menu_copy_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_copyall_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copyall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
raw=handles.alluserdat.rawdata;
tempidx=get(handles.popup_sel_housekeeping,'Value');
templist=get(handles.popup_sel_housekeeping,'String');
ctrlTar{1}=templist{tempidx};
tempidx=get(handles.popup_selgrp_ctrl,'Value');
grpsamps=handles.alluserdat.grp_samples;
ctrlSamp=tempidx;%grpsamps{tempidx};
grpnames=get(handles.popup_selgrp,'String');

spec_gp=handles.alluserdat.graph_param;
spec_gp.no_graph=true;

package=gui_proc_QuantStudio6_qpcrdata(raw,ctrlTar,handles.alluserdat.graph_param.rem_tarlist,ctrlSamp,grpnames,grpsamps,spec_gp);

text_mean=dat2D_to_str(package.mean,package.SampleNames);
if length(size(package.sd))==2
    text_sd=dat2D_to_str(package.sd,package.SampleNames);
elseif length(size(package.sd))==3
    text_sd_low=dat2D_to_str(package.sd(:,:,1),package.SampleNames);
    text_sd_high=dat2D_to_str(package.sd(:,:,2),package.SampleNames);
    text_sd=[text_sd_low,sprintf('\t<--to this point = lower bounds\r'),text_sd_high,sprintf('\t<-- to this point = higher bounds')];
end
text_tar='';
for n=1:length(package.TargetNames)
    if n<length(package.TargetNames)
        text_tar=[text_tar,package.TargetNames{n},sprintf('\t')];
    else
        text_tar=[text_tar,package.TargetNames{n}];
    end
end

alltext=[sprintf('MEAN\t'),text_tar,sprintf('\r'),text_mean,sprintf('\r\r'),sprintf('STD-DEV\t'),text_tar,sprintf('\r'),text_sd];
clipboard('copy',alltext)


% --------------------------------------------------------------------
function menu_copygene_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copygene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
raw=handles.alluserdat.rawdata;
tempidx=get(handles.popup_sel_housekeeping,'Value');
templist=get(handles.popup_sel_housekeeping,'String');
ctrlTar{1}=templist{tempidx};
tempidx=get(handles.popup_selgrp_ctrl,'Value');
grpsamps=handles.alluserdat.grp_samples;
ctrlSamp=tempidx;%grpsamps{tempidx};
grpnames=get(handles.popup_selgrp,'String');

spec_gp=handles.alluserdat.graph_param;
spec_gp.no_graph=true;

package=gui_proc_QuantStudio6_qpcrdata(raw,ctrlTar,handles.alluserdat.graph_param.rem_tarlist,ctrlSamp,grpnames,grpsamps,spec_gp);

text_tar='';
for n=1:length(package.TargetNames)
    if n<length(package.TargetNames)
        text_tar=[text_tar,package.TargetNames{n},sprintf('\t')];
    else
        text_tar=[text_tar,package.TargetNames{n}];
    end
end
clipboard('copy',text_tar)

% --------------------------------------------------------------------
function menu_copygrp_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copygrp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
raw=handles.alluserdat.rawdata;
tempidx=get(handles.popup_sel_housekeeping,'Value');
templist=get(handles.popup_sel_housekeeping,'String');
ctrlTar{1}=templist{tempidx};
tempidx=get(handles.popup_selgrp_ctrl,'Value');
grpsamps=handles.alluserdat.grp_samples;
ctrlSamp=tempidx;%grpsamps{tempidx};
grpnames=get(handles.popup_selgrp,'String');

spec_gp=handles.alluserdat.graph_param;
spec_gp.no_graph=true;

package=gui_proc_QuantStudio6_qpcrdata(raw,ctrlTar,handles.alluserdat.graph_param.rem_tarlist,ctrlSamp,grpnames,grpsamps,spec_gp);

text_grp='';
for n=1:length(package.SampleNames)
    text_grp=[text_grp,package.SampleNames{n},sprintf('\r')];
end
clipboard('copy',text_grp)

% --------------------------------------------------------------------
function menu_copymean_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copymean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
raw=handles.alluserdat.rawdata;
tempidx=get(handles.popup_sel_housekeeping,'Value');
templist=get(handles.popup_sel_housekeeping,'String');
ctrlTar{1}=templist{tempidx};
tempidx=get(handles.popup_selgrp_ctrl,'Value');
grpsamps=handles.alluserdat.grp_samples;
ctrlSamp=tempidx;%grpsamps{tempidx};
grpnames=get(handles.popup_selgrp,'String');

spec_gp=handles.alluserdat.graph_param;
spec_gp.no_graph=true;

package=gui_proc_QuantStudio6_qpcrdata(raw,ctrlTar,handles.alluserdat.graph_param.rem_tarlist,ctrlSamp,grpnames,grpsamps,spec_gp);
text_mean=dat2D_to_str(package.mean,package.SampleNames);

text_tar='';
for n=1:length(package.TargetNames)
    if n<length(package.TargetNames)
        text_tar=[text_tar,package.TargetNames{n},sprintf('\t')];
    else
        text_tar=[text_tar,package.TargetNames{n}];
    end
end
alltext=[sprintf('MEAN\t'),text_tar,sprintf('\r'),text_mean];
clipboard('copy',alltext)

% --------------------------------------------------------------------
function menu_copysd_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copysd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
raw=handles.alluserdat.rawdata;
tempidx=get(handles.popup_sel_housekeeping,'Value');
templist=get(handles.popup_sel_housekeeping,'String');
ctrlTar{1}=templist{tempidx};
tempidx=get(handles.popup_selgrp_ctrl,'Value');
grpsamps=handles.alluserdat.grp_samples;
ctrlSamp=tempidx;%grpsamps{tempidx};
grpnames=get(handles.popup_selgrp,'String');

spec_gp=handles.alluserdat.graph_param;
spec_gp.no_graph=true;

package=gui_proc_QuantStudio6_qpcrdata(raw,ctrlTar,handles.alluserdat.graph_param.rem_tarlist,ctrlSamp,grpnames,grpsamps,spec_gp);

if length(size(package.sd))==2
    text_sd=dat2D_to_str(package.sd,package.SampleNames);
elseif length(size(package.sd))==3
    text_sd_low=dat2D_to_str(package.sd(:,:,1),package.SampleNames);
    text_sd_high=dat2D_to_str(package.sd(:,:,2),package.SampleNames);
    text_sd=[text_sd_low,sprintf('<--to this point = lower bounds\r'),text_sd_high,sprintf('<-- to this point = higher bounds')];
end

text_tar='';
for n=1:length(package.TargetNames)
    if n<length(package.TargetNames)
        text_tar=[text_tar,package.TargetNames{n},sprintf('\t')];
    else
        text_tar=[text_tar,package.TargetNames{n}];
    end
end
alltext=[sprintf('STD-DEV\t'),text_tar,sprintf('\r'),text_sd];
clipboard('copy',alltext)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%my functions 

function handles = update_early_guis(hObject,handles)

set(handles.popup_sel_housekeeping,'String',handles.alluserdat.rawdata.TarNames)

%building secondary normalization list:
prim_idx=get(handles.popup_sel_housekeeping,'Value');
if prim_idx==1
    rem_list=handles.alluserdat.rawdata.TarNames(2:end);
elseif prim_idx==length(handles.alluserdat.rawdata.TarNames)
    rem_list=handles.alluserdat.rawdata.TarNames(1:(end-1));
else
    rem_list=[handles.alluserdat.rawdata.TarNames(1:(prim_idx-1));handles.alluserdat.rawdata.TarNames((prim_idx+1):end)];
end
set(handles.list_secondary_norm,'String',['none';rem_list])

set(handles.list_avail_samp,'String',handles.alluserdat.rawdata.SampleNames)

tarnames=rem_list;
text=sprintf('    Gene List:\n');
for n=1:length(tarnames)
    addtext=sprintf([num2str(n),'. ',tarnames{n},'\n']);
    text=[text,addtext];    
end
set(handles.text_tarlist,'String',text)
handles.alluserdat.graph_param.rem_tarlist=rem_list;
guidata(hObject,handles)

function handles=update_numgrp_guis(hObject,handles)
grpnames=handles.alluserdat.grp_names;
grpsamps=handles.alluserdat.grp_samples;
if isempty(grpnames)
    %initialize grps
    for n=1:handles.alluserdat.num_grps
        grpnames{n}=num2str(n);
        grpsamps{n}={};
    end
elseif length(grpnames)<handles.alluserdat.num_grps
    for n=(length(grpnames)+1):handles.alluserdat.num_grps
        grpnames{n}=num2str(n);
        grpsamps{n}={};
    end
elseif length(grpnames)>handles.alluserdat.num_grps
    %be sure to create warning/confirmation dlg
    grpnames=grpnames(1:handles.alluserdat.num_grps);
    grpsamps=grpsamps(1:handles.alluserdat.num_grps);
    set(handles.popup_selgrp,'Value',1);
    set(handles.popup_selgrp_ctrl,'Value',1);
end
handles.alluserdat.grp_names=grpnames;
handles.alluserdat.grp_samples=grpsamps;

if isempty(grpnames)
    curgrpidx=0;
    set(handles.edit_numgrps,'String','0');
    set(handles.edit_grpname,'String','');
    set(handles.popup_selgrp,'String',' ');
    set(handles.popup_selgrp_ctrl,'String',' ');
    set(handles.popup_selgrp,'Value',1);
    set(handles.popup_selgrp_ctrl,'Value',1);
    
else
    curgrpidx=get(handles.popup_selgrp,'Value');
    set(handles.edit_grpname,'String',grpnames{curgrpidx});
    set(handles.popup_selgrp,'String',grpnames);
    set(handles.popup_selgrp_ctrl,'String',grpnames);
end


update_grp_list(hObject,handles)
guidata(hObject,handles)

function update_curgrp_guis(hObject,handles)
grpnames=handles.alluserdat.grp_names;
grpsamps=handles.alluserdat.grp_samples;

if isempty(grpnames)
    curgrpidx=0;
    set(handles.edit_grpname,'String','');
    set(handles.list_chosen_samp,'String','')
    set(handles.popup_selgrp_ctrl,'String',' ');
else
    curgrpidx=get(handles.popup_selgrp,'Value');
    set(handles.edit_grpname,'String',grpnames{curgrpidx});
    set(handles.list_chosen_samp,'String',grpsamps{curgrpidx})
    set(handles.popup_selgrp_ctrl,'String',grpnames);
end
update_grp_list(hObject,handles)
guidata(hObject,handles)

function update_grp_list(hObject,handles)
grpnames=handles.alluserdat.grp_names;
text=sprintf('    Group List:\n');
for n=1:length(grpnames)
    addtext=sprintf([num2str(n),'. ',grpnames{n},'\n']);
    text=[text,addtext];    
end
set(handles.text_grplist,'String',text)

function update_graph_params(hObject,handles)
set(handles.check_geometric,'Value',handles.alluserdat.graph_param.geo)
set(handles.check_log,'Value',handles.alluserdat.graph_param.log)
set(handles.check_legend,'Value',handles.alluserdat.graph_param.legend)

set(handles.edit_selgrps2plot,'String',num2str(handles.alluserdat.graph_param.grps2plot))
set(handles.edit_seltar2plot,'String',num2str(handles.alluserdat.graph_param.tars2plot))


function outtext=dat2D_to_str(indat,col)
outtext='';
sz=size(indat);
if nargin==1
    for m=1:sz(1)
        for n=1:sz(2)
            outtext=[outtext,sprintf([num2str(indat(m,n)),'\t'])];
        end
        if m<sz(1)
            outtext=[outtext,sprintf('\r')];
        end
    end
else
    for m=1:sz(1)
        outtext=[outtext,sprintf([col{m},'\t'])];
        for n=1:sz(2)
            if n<sz(2)
                outtext=[outtext,sprintf([num2str(indat(m,n)),'\t'])];
            else
                outtext=[outtext,num2str(indat(m,n))];
            end
        end
        if m<sz(1)
            outtext=[outtext,sprintf('\r')];
        end
    end
end