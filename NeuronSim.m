function varargout = NeuronSim(varargin)
% NEURONSIM MATLAB code for NeuronSim.fig
%      NEURONSIM, by itself, creates a new NEURONSIM or raises the existing
%      singleton*.
%
%      H = NEURONSIM returns the handle to a new NEURONSIM or the handle to
%      the existing singleton*.
%
%      NEURONSIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURONSIM.M with the given input arguments.
%
%      NEURONSIM('Property','Value',...) creates a new NEURONSIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuronSim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuronSim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuronSim

% Last Modified by GUIDE v2.5 04-Feb-2016 02:31:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuronSim_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuronSim_OutputFcn, ...
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


% --- Executes just before NeuronSim is made visible.
function NeuronSim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuronSim (see VARARGIN)

% Choose default command line output for NeuronSim
handles.output = hObject;

handles.noise_60Hz=0;
handles.noise_TCS=0;
handles.input_time=5;
handles.sampling_freq=5;
handles.neurons=1;
handles.isi=1;
handles.snr=1;
handles.feedback=0;



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NeuronSim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuronSim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in OPTION_feedback.
function OPTION_feedback_Callback(hObject, eventdata, handles)
% hObject    handle to OPTION_feedback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OPTION_feedback
handles.feedback = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in OPTION_60hz_noise.
function OPTION_60hz_noise_Callback(hObject, eventdata, handles)
% hObject    handle to OPTION_60hz_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OPTION_60hz_noise
handles.noise_60Hz = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in OPTION_TCS_noise.
function OPTION_TCS_noise_Callback(hObject, eventdata, handles)
% hObject    handle to OPTION_TCS_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OPTION_TCS_noise
handles.noise_TCS = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);


function INPUT_time_Callback(hObject, eventdata, handles)
% hObject    handle to INPUT_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of INPUT_time as text
%        str2double(get(hObject,'String')) returns contents of INPUT_time as a double
handles.input_time = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function INPUT_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to INPUT_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function INPUT_Fs_Callback(hObject, eventdata, handles)
% hObject    handle to INPUT_Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of INPUT_Fs as text
%        str2double(get(hObject,'String')) returns contents of INPUT_Fs as a double
handles.sampling_freq=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function INPUT_Fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to INPUT_Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function INPUT_neurons_Callback(hObject, eventdata, handles)
% hObject    handle to INPUT_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of INPUT_neurons as text
%        str2double(get(hObject,'String')) returns contents of INPUT_neurons as a double
handles.neurons=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function INPUT_neurons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to INPUT_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function INPUT_isi_Callback(hObject, eventdata, handles)
% hObject    handle to INPUT_isi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of INPUT_isi as text
%        str2double(get(hObject,'String')) returns contents of INPUT_isi as a double
array_input = get(hObject, 'String');
array_converted = str2num(char(array_input));
handles.isi=array_converted;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function INPUT_isi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to INPUT_isi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function INPUT_snr_Callback(hObject, eventdata, handles)
% hObject    handle to INPUT_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of INPUT_snr as text
%        str2double(get(hObject,'String')) returns contents of INPUT_snr as a double
array_input = get(hObject, 'String');
array_converted = str2num(char(array_input));
handles.snr=array_converted;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function INPUT_snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to INPUT_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in INPUT_start.
function INPUT_start_Callback(hObject, eventdata, handles)
% hObject    handle to INPUT_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Call Options and Inputs
make_60hz = handles.noise_60Hz;
make_tcs = handles.noise_TCS;
signal_time=handles.input_time;
fs=handles.sampling_freq;
neurons=handles.neurons;
model_fs=100000;
isi=handles.isi;
snr=handles.snr;
play_sound=handles.feedback;

%Generate Model Neuron Template 1ms @ 100 kHz 
neuron_model_0=nan(1,100);
neuron_model_1=[0,0,0,0,0,0,0,0,0,.02,.05,.08,.1,.104,.110,.115,.123,.205,.308,.405,.605,.706,.89,.98,1,.98,.89,.706,.505,.305,.105,-.155,-.310,-.525,-.712,-.87,-.93,-.97,-.99,-1,-.99,-.97,-.93,-.87,-.712,-.525,-.310,-.155,.08,.11,.14,.15,.155,.158,.16,.161,.16,.158,.155,.15,.14,.132,.127,.122,.115,.114,.111,.108,.106,.100,.095,.09,.085,.08,.075,.06,.05,.04,.03,.02,.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

%Generate Time Vector
model_time_vector = 0:1/model_fs:signal_time;
total_datapoints = length(model_time_vector);

%Pre-Allocate Final Noise Matrix
Final_Noise=zeros(1,total_datapoints);

%Pre-Allocate and Generate Noise Matrix
if make_60hz == 1 && make_tcs == 1
    noise_chans = zeros(4,total_datapoints);
    noise_chans(1,:) = sin(2*pi*model_time_vector*60);
    noise_chans(2,:) = wgn(1,total_datapoints,0);
    noise_chans(3,:) = wgn(1,total_datapoints,0);
    noise_chans(4,:) = wgn(1,total_datapoints,0);
    Final_Noise = noise_chans(1,:)+noise_chans(2,:)+noise_chans(3,:)+noise_chans(4,:);
elseif make_60hz == 1
    Final_Noise = sin(2*pi*model_time_vector*60);
elseif make_tcs == 1
    noise_chans = zeros(3,total_datapoints);
    noise_chans(1,:) = wgn(1,total_datapoints,0);
    noise_chans(2,:) = wgn(1,total_datapoints,0);
    noise_chans(3,:) = wgn(1,total_datapoints,0);
    Final_Noise = noise_chans(1,:)+noise_chans(2,:)+noise_chans(3,:);
end

%Pre-Allocate Neuron Firing Matrix
spike_data=zeros(neurons,(100000*signal_time));

%Generate Placement Raster Data for Neurons @ 100 kHz
time_step = 1/100000;
for i = 1:neurons
    spike_data(i,:) = rand(1, (100000*signal_time)) < isi(1,i)*time_step;
end

%Pre-Allocate Final Neuron Data Matrix
Neuron_Matrix=nan(neurons,total_datapoints);
write_tracker=zeros(neurons,2);
%Reconstruct neuron vectors based on firing patterns generated
for i = 1:(100000*signal_time)

    for j = 1:neurons
        
        if write_tracker(j,1)==0
            %check to write
            if spike_data(j,i)==1
                write_tracker(j,1)=1;
                write_tracker(j,2)=1;
                
                Neuron_Matrix(j,i)=snr(1,j)*neuron_model_1(1,write_tracker(j,2));
                write_tracker(j,2)=write_tracker(j,2)+1;
                
            else
                    Neuron_Matrix(j,i)=NaN;
            end
            
        else
                if write_tracker(j,2) > 99  %currently all neuron model sizes are 100 datapoints so we're writing from 0 to 99.
                    write_tracker(j,1)=0;
                    Neuron_Matrix(j,i)=NaN;
                else
                    Neuron_Matrix(j,i)=snr(1,j)*neuron_model_1(1,write_tracker(j,2));
                    write_tracker(j,2)=write_tracker(j,2)+1;
                end
        end
    end

end

%Combine Final Neuron Vector
[x,~]=size(Neuron_Matrix);
if x == 1
    Final_Neurons=Neuron_Matrix;
else
    Final_Neurons=nansum(Neuron_Matrix);
end


%Generate Plot Data Vector @ 100 kHz
Oversampled_Data=nan(2,total_datapoints);
Oversampled_Data(1,:) = Final_Noise;
Oversampled_Data(2,:) = Final_Neurons;

Oversampled_Plot_Data=nansum(Oversampled_Data);

%Downsample to specified recording frequency
down_sample_value=floor(model_fs/fs);
Final_Plot_Data = downsample(Oversampled_Plot_Data,down_sample_value);

%Plot Data
figure;plot(Final_Plot_Data);

if play_sound==true
    soundsc(Final_Plot_Data,fs);
end
