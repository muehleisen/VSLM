  function varargout = vslm(varargin)
% VSLM_M-file for vslm.fig
%      VSLM, by itself, creates a new VSLM or raises the existing
%      singleton*.
%
%      H = VSLM returns the handle to a new VSLM or the handle to
%      the existing singleton*.
%
%      VSLM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VSLM.M with the given input arguments.
%
%      VSLM('Property','Value',...) creates a new VSLM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vslm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vslm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE'Slow Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to xproducthelp vslm

% V0.4.2 Update 01-Oct-2016 move code from sourceforge to github and start using git revision control.  Remove explicit references to v0_4_1 and now use revision
% numbering with release and in these comments info but not within code filenameesitself%
% V0.4.2 Update 16-Aug-2016 by Elsa Piollet, <elsa.piollet@polymtl.ca: changed wavread to audioread for compatibility with Matlab 2015 and later
% V0.4.1 Initial public release on sourceforge in 2011

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vslm_OpeningFcn, ...
                   'gui_OutputFcn',  @vslm_OutputFcn, ...
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


% --- Executes just before vslm is made visible.
function vslm_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vslm (see VARARGIN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%
if 1==1
    disp('GUIstarted'); %tells the Splash screen that everything is ok and it may close.
else
    disp ('GUIerror'); %tells the Splash screen that something went wrong and it should close.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%


% Choose default command line output for vslm
handles.output = hObject;

% set max array size to load in one hunk.
% With sampling rates of 22055, 44100, 48k and 96k, the smallest time slice
% with an integer number of samples for each sampling rate is 20 ms
% for 20 ms there is 441 samples at 22k and 1920 at 96k.  Make the max
% array size a multiple 441*1920
handles.meas.maxsize=2540160;  % this number is 3x 441*1920
handles.meas.mindt=0.02; % set the minimum time slice for further analysis


% set default settings for meas file, cal file, measurement length, sample
% rate and calibration factor
handles.cal.file='None Loaded';
handles.cal.maxlen=30;  % set the max cal file length to 30s 
handles.cal.scale=1.0;

handles.meas.Nseg=1;  % set the default number of buffers
handles.meas.file='None Loaded';
handles.meas.len=0;

handles.meas.fs=0;
handles.separatedplot=0;  % default with no separate display figure

% update the text info string 
handles.infotxt(1)={handles.meas.file};
handles.infotxt(2)={[num2str(handles.meas.len) ' sec']};
handles.infotxt(3)={[num2str(handles.meas.fs) ' Hz']};
handles.infotxt(4)={handles.cal.file};
handles.infotxt(5)={num2str(handles.cal.scale)};
set(handles.text4,'string',handles.infotxt);


% set up default analysis settings
handles.weighting='A';  % default to A weighting
set(handles.A,'value',1); % set the button to on
handles.speed='Slow';  % default to slow speed
set(handles.Slow,'value',1); % set the button to on 
handles.analysis='lpbtn';  % default to sound pressure level
set(handles.xlpbtn,'value',1);


% set up lpplot defaults
handles.lp.dt=0.100;  % set default lp plot time spacing to 100 ms
handles.lp.dth=handles.xlpdt100ms;
handles.lp.plotscale=0; % default lp plot scale to autoscale


% set up xleqbtn defaults
handles.leq.mindt=0.1;
handles.leq.dt=1; % default user LEQ  time segment to 1 s
handles.leq.dth=handles.xleq1s; % set handle 
handles.leq.userpct=5; % default LEQ user percentile to 5
handles.leq.plotscale=0;  % default to autoscale
handles.leq.Dtype='NIOSH'; % set noise dose to NIOSH by default
handles.leq.Dtypeh=handles.xniosh;
handles.leq.Dex=3; % set the user noise doise exchange rate to NIOSH value of 3 dB by default
handles.leq.DLc=85; % set the user noise dose criterion to NIOSH value of 85 by default
handles.leq.DLt=80; % set the user noise dose threshold to NIOSH value of 80 by default

% set up band defaults
handles.band.type='oct';
handles.band.typeh=handles.xoct;
handles.band.method='FFT';
handles.band.methodh=handles.xbandfft;

% set up xpsdbtn defaults
handles.psd.Nfft=4096; % set default PSD FFT size to 4096
handles.psd.Nffth=handles.xpsd4k; % set handle to point to 4k fft
handles.psd.Nwindow=handles.psd.Nfft; % default window size to FFT Size
handles.psd.window='Hanning'; % set the default PSD window to be a hanning window
handles.psd.windowh=handles.xpsdhann;
handles.psd.overlap=0.50; % set default PDF overlap to 50%
handles.psd.overh=handles.xpsdover50;

% set up spectrogram defaults
handles.spec.Nfft=512; % set default fft for spectrogram to 512 points
handles.spec.Nffth=handles.xspec512; % set menu to point to 512 pt fft
handles.spec.dt=1.0;  % set default spectrogram time slice length to 1 sec
handles.spec.dth=handles.xspecdt1s; % set menu to point to 1 s
handles.spec.dnr=0; % set default spectrogram range to auto
handles.spec.dnrh=handles.xspecdnrauto; % set menu to point to auto
handles.spec.view='2d';  % set the default view to 2D
handles.spec.viewh=handles.xspecview2d; % set menu to point to auto

% clear the plot display when the program starts
subplot(1,1,1,'Parent',handles.plots)
cla reset
axis off


handles.last.analysis='none';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vslm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vslm_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in quit.
function quit_Callback(~, ~, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% delete the main gui figure handle to quit the gui and program
delete(handles.figure1)

% --- Executes on button press in playsound.
function playsound_Callback(hObject, ~, handles)
% hObject    handle to playsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispc
    if ~strcmp(handles.meas.file,'None Loaded')       
        m=figure;
        x=get(m,'position');
        set(m,'toolbar','none','menubar','none','position',[x(1) x(2) 300 48] )
        h=actxcontrol('WMPlayer.OCX.7', [0 0 300 48], m);
        h.URL=handles.meas.fullfile;
        h.controls.play;     
    else
        dstring={'Cannot Play Measurement File','No Measurement File Loaded'};
        plotdisplay(handles,dstring);
    end
else
    
    dstring={'Playing Measurement Files only supported on PC'};
    plotdisplay(handles,dstring);
end


pause(1)
guidata(hObject, handles);


% --- Executes on button press in analyzebtn.
function analyzebtn_Callback(hObject, ~, handles)
% hObject    handle to analyzebtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~strcmp(handles.meas.file,'None Loaded')
    switch handles.analysis
        case 'xspecbtn'
            %do the specgtrogram calculation
            analyze_spec(hObject,handles);
        case 'xleqbtn'
            % do the LEQ computation
            analyze_leq(hObject,handles);
        case 'xbandbtn'
            % do octave or 1/3 octave band analysis
            if handles.band.method(1)=='F'
                analyze_bandfft(hObject,handles);
            else
                analyze_bandansi(hObject,handles);
            end
        case 'xncrcbtn'
            % do NC/RC analysis
            analyze_ncrc(hObject,handles);
        case 'xpsdbtn'
            % do the power spectral density calculation
            analyze_psd(hObject,handles);           
        otherwise
            % if we get here we must be asking to do an Lp computation
            analyze_lp(hObject,handles);
    end
else %if we get here we haven't loaded a file into memory   
    dstring={'Cannot Analyze Data','No Measurement File Loaded'};
    plotdisplay(handles,dstring);
end

function outdata=apply_weighting(indata,fs,weighting)
allowed_fs=[22050,44100,48000,96000];
I=find(fs==allowed_fs);
switch weighting
    case {'A'}
        % scale by cal factor and filter with A-weighted filter
        
        % created with adesign4e22k and N=1399, adesign4e44k and N=66,
        % adesign4e48k and N=1331, adesign4e96k and N=???
        a_coeff=[1.00000000000000,-0.820486371642917,-0.954967112698490,0.819769137019238,-0.0308133150879904;...
            1.00000000000000,-3.16548213635874,3.58370028028808,-1.66952594339746,0.251312323839929;...
            1.00000000000000,-3.21098152161547,3.70934390304575,-1.78459093140491,0.286231936413823;...
            1.00000000000000,-3.52292032961222,4.59276248551295,-2.61658164419795,0.546739755260711];
        
        b_coeff=[0.411130459349308,0.489163899078409,-1.81511690197873,0.505982131319440,0.408895901518515;...
            0.218436284008832,-0.0271197371604133,-1.23075568024016,1.66912752925371,-0.629688397944900;...
            0.180459241142510,0.0860957271342434,-1.34233754417886,1.70455253682352,-0.628769962323045;...
            -0.0930424680674304,0.885061570433706,-2.09736158642518,1.91170857911392,-0.606366095212699];
        a=a_coeff(I,:);
        b=b_coeff(I,:);
        outdata=filter(b,a,indata);
    case 'C'
        %  filter with the C weighting filter     
        %created with cdesign4e22k and N=1470, cdesign4e44k and N=1186,
        %cdesign4e48k and N=1102, and cdesign4e96k and N=1201
        a_coeff=[1.00000000000000,0.487110438920731,-1.19411745882897,-0.473099276797724,0.208201319236022; ...
            1.00000000000000,0.720524473037467,-1.27460941890720,-0.711991070511022,0.283142933498142; ...
            1.00000000000000,0.678016119469390,-1.28615510047168,-0.687899298538665,0.312313266344913;...
            1.00000000000000,-2.93939917852114,3.10634442611301,-1.39372064867588,0.226775918992396];
        b_coeff=[0.403673160421101,0.867939940009175,-0.387869031048492,-0.871337058723055,-0.0191991534215092; ...
            0.205377794209128,0.731122249663226,0.318817527196078,-0.733113968701055,-0.526187064231238;...
            0.210715717951749,0.718855215372106,0.245507756808660,-0.716973939314165,-0.462053117024141;...
            0.0484465929621835,0.151813149050837,-0.459174213032651,0.269122623377235,-0.0102081523621007];
        a=a_coeff(I,:);
        b=b_coeff(I,:);
        outdata=filter(b,a,indata);
        
    otherwise
        % if we get here we have flat weighting and we do no filtering
        outdata=indata;
end % switch weighting

% --- Executes on button press in set_cal.
function set_cal_Callback(hObject, ~, handles)
% hObject    handle to set_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string=inputdlg('Input Calibration Factor','Set Cal Value');
if ~isempty(string)
    % if the string is not empty, convert
    handles.cal.scale=str2double(string{1});
    handles.cal.file='None Loaded';
    handles.infotxt(4)={handles.cal.file};
    handles.infotxt(5)={num2str(handles.cal.scale)};
    set(handles.text4,'string',handles.infotxt); 
    dstring=sprintf('Calibration Factor Set to %G',handles.cal.scale);
    plotdisplay(handles,dstring);
end

guidata(hObject, handles);

% --- Executes on button press in loadcalwav.
function loadcalwav_Callback(hObject, ~, handles)
% hObject    handle to loadcalwav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile('.wav','Open Calibration File');
switch filename
    case {0}
        % User cancelled out
    otherwise
        fullfilename=[pathname,filename];
        % get the sample rate and filesize
        [~,fs]=audioread(fullfilename,1);
        fsize=size(audioread(fullfilename));
        
        if all(fs~=[22050,44100,48000,96000])           
            % if we get here, fs is not one of our four allowed fs values
            h=warndlg('Only files with fs=22.05k, 44.1k, 48k and 96k sampling rate allowed');
            waitfor(h);s
            return %if we got here we cancelled out so end routine
        end % if all(fs~=[22050,44100,48000,96000])  
             
        string=inputdlg('Input Calibrator RMS Output in DB','Set Cal Value');
        if ~isempty(string{1});
            % convert input string to number
            caldB=str2double(string{1});
            % generate a linear cal factor from the xlpbtn level
            calpa = 10.0.^((caldB-94)/20);
            
            if fsize(2)>1;
                h=warndlg('Warning, Multiple Channels Detected.  Only Channel 1 is read');
                waitfor(h);
            end %if fsize(2)>1;
            
            % check to see if the file is longer than the max allowed
            % length and set the max size to read in if it is
            Nsamp=handles.cal.maxlen*fs;
            if fsize(1)>handles.cal.maxlen*fs;
                dstring=sprintf('Calibration file truncated to %d sec',handles.cal.maxlen);
                h=warndlg(dstring,'Cal File Truncated');
                waitfor(h);
            else
                % if filesize is less than the max
                Nsamp=fsize(1);
            end %if fsize(1)>handles.cal.maxlen*fs;            
            
            data=audioread(fullfilename,[1,Nsamp]);
            caldata=data(:,1);
            calrms=sqrt(mean(abs(caldata).^2));
            calfactor=calpa/calrms;
            handles.infotxt(4)={filename};
            handles.infotxt(5)={num2str(calfactor)};
            set(handles.text4,'string',handles.infotxt);
            
            dstring={['Cal File Name:',filename],...
                ['Cal File Length: ',num2str(round(10*length(caldata)/fs)/10),' sec'], ...
                ['Cal Level: ',num2str(caldB),' dB'], ...
                ['Cal File RMS: ',num2str(calrms)], ...
                ['Cal Factor: ',num2str(calfactor)]};
            plotdisplay(handles,dstring);
            
            handles.cal.file=fullfilename;
            handles.cal.scale=calfactor;
            
        end % if ~isempty(string{1});  
end % switch filename
    
% Update handles structure
guidata(hObject, handles);
    
% --------------------------------------------------------------------

% --- Executes on button press in loadmeaswav.
function loadmeaswav_Callback(hObject, ~, handles)
% hObject    handle to loadmeaswav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile('.wav','Open Measurement File');
switch filename
    case {0}
        % If we get here the user cancelled out - do nothing
    otherwise
        fullfilename=[pathname,filename];  % create the full filename for subsequent processing
        [~,fs]=audioread(fullfilename,[1,1]); % read the first sample to get the sample rate
        if all(fs~=[22050,44100,48000,96000])
            
            % if we get here, fs is not one of our four allowed fs values
            h=warndlg('Only files with fs=22.05k, 44.1k, 48k and 96k sampling rate allowed');
            waitfor(h);
 %           filename=0;  % reset filename to 0 as if we cancelled out
            return
        end
        % 20 ms segments are the smallest ones that give an integer number
        % of samples for all four sample rates.  Break everything into
        % blocks of 20 ms
        bsize=fs*0.02;  % set the basic block size for a 20 ms segment
        
        % find the file size of file truncated to integer number bsize segments
        filesize=size(audioread(fullfilename));  % get file size in samples
        if filesize(2)>1
            h=warndlg('Multiple Channels Detected.  Only Channel 1 will be Analyzed');
            waitfor(h)
        end
        fsize=floor(filesize(1)/bsize)*bsize;  % truncate to a multiple of bsize
        measlen=fsize/fs; % compute and store the length in s
        % update the info display string with the meas filename
        
        % The following code figures out how many segments we'll break the
        % file into for later processing
        maxsize=handles.meas.maxsize;   % read in the maximum block size
        
        if fsize>maxsize
            % break into segments of maxsize samples and analyze the file in parts
            Nseg=floor(fsize/maxsize);  % compute the number of segments
            N=1:maxsize:fsize; % generate sample indexes into file for each segment
            if (fsize/maxsize)~=floor(fsize/maxsize)
                % if the filesize is not a multiple of maxsize, add the last seg
                N=[N,fsize];
                Nseg=Nseg+1;
            end
        else
            Nseg=1;
            N=[1,fsize];
        end
        
        % update the info text and the main display
        handles.infotxt(1)={filename};
        handles.infotxt(2)={[num2str(measlen) ' sec']};
        handles.infotxt(3)={[num2str(fs) ' Hz']};
        set(handles.text4,'string',handles.infotxt);
        
        dstring={['Meas File Name:',filename];...
            ['Meas File Length: ',num2str(measlen),' sec']};
        plotdisplay(handles,dstring);
        
        handles.meas.fs=fs;
        handles.meas.fullfile=fullfilename;
        handles.meas.file=filename;
        handles.meas.path=path;
        handles.meas.len=measlen;
        handles.meas.Nseg=Nseg;
        handles.meas.N=N;
        lp.dthsize=fsize;
        handles.meas.bsize=bsize;
end % switch

% Update handles structure
guidata(hObject, handles);

function xfilemenu_Callback(~, ~, ~)
% hObject    handle to xfilemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function File_CreateFcn(~, ~, ~)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in wtgbtns.
function wtgbtns_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in wtgbtns 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles.weighting=get(eventdata.NewValue,'Tag');
guidata(hObject, handles);


% --- Executes when selected object is changed in speedbtns.
function speedbtns_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in speedbtns 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles.speed=get(eventdata.NewValue,'Tag');
guidata(hObject, handles);

% --- Executes when selected object is changed in analysis.
function analysis_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in analysis 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles.analysis=get(eventdata.NewValue,'Tag');
guidata(hObject, handles);

function analyze_lp(hObject, handles)
filename=handles.meas.fullfile;
Nseg=handles.meas.Nseg;
N=handles.meas.N;
fs=handles.meas.fs;
udt=handles.lp.dt; % load up the user defined dt for display
%weighting=handles.weighting;
cal=handles.cal.scale;

hw=waitbar(0,'Please Wait ...');
dI=floor(fs.*udt);  % compute the sample spacing for uder display
d2last=0; % initialize a variable to hold the last data point in a given segment

% read the file in large block sizes and compute Lp using our ARMA filter
% and save out all the points every udt seconds.
p2=zeros(floor(N(end)/dI),1);
Klast=0;  % starting index for the square pressure 
ppeak=0;
tpeak=0;
for I=1:Nseg
    data=cal.*audioread(filename,N(I:I+1));
    data2=(apply_weighting(data(:,1),fs,handles.weighting)).^2;
 
    % because we are reading the file in separated segments we cannot simply use
    % filter.  Also, the impulse meter has a response that changes depending upon if
    % the signal is rising and falling.  So, let'slow just implement the
    % filters by hand with for loops
    switch handles.speed
        case 'Slow'
            wbstring='Please Wait ... Simulating Slow Meter';
            % set the exponential time constant to 1 sec and compute the
            % digital filter coefficents
            tau=1.0;
            b1=1.0-exp(-1/(fs*tau)); a2=b1-1;
            % implement the ARMA filter.  The first output uses the
            % last datapoint from the last segment so do that separately
            data2(1)=b1.*data2(1)-a2.*d2last;
            for J=2:length(data2)
                data2(J)=b1*data2(J)-a2*data2(J-1);
            end % for J=2:length(data2)
            d2last=data2(end);
            % square the data then put through exponential filter
        case 'Fast'
            wbstring='Please Wait ... Simulating Fast Meter';
            tau=0.125;  % set fast meter time constant to 125 ms.
            b1=1.0-exp(-1/(fs*tau));  % compute filter coefficients
            a2=b1-1;
            % do the 1st datapoint separately
            data2(1)=b1.*data2(1)-a2.*d2last;
            for J=2:length(data2)
                data2(J)=b1*data2(J)-a2*data2(J-1);
            end % for J=2:length(data2)
            d2last=data2(end);
        otherwise
            wbstring='Please Wait ... Simulating Impulse Meter';
            taur=0.035; % set rising pressure time constant is 35 ms
            tauf=1.5; % set falling pressure time constant to 1.5 s
            br=1.0-exp(-1/(fs*taur));  % get rising pressure b and a
            ar=br-1.0;
            bf=1.0-exp(-1/(fs*tauf));  % get falling pressure b and a
            af=bf-1.0;
            
            % do the IIR filter using the fast time constant on the rising part
            % of the signal and the slow time constant on the falling part
            
            % do the 1st datapoint separately using
            if data2(1)>d2last
                b1=br; a2=ar;
            else
                b1=bf; a2=af;
            end
            data2(1)=b1.*data2(1)-a2.*d2last;
            for J=2:length(data2)
                if data2(J)>=data2(J-1)  % if data2(J)> data2(J-1) the data is rising
                    b1=br; a2=ar;  % set the filter coefficients to the rising ones
                else % otherwise it is falling
                    b1=bf; a2=af;  % set the filter coefficients to the falling ones
                end
                data2(J)=b1*data2(J)-a2*data2(J-1);  % filter the data using the selected coefficients
            end % for J=2:length(data2)
            d2last=data2(end);  % save the last data point for filtering with the next segment
    end % switch handles.speed
    
    % this code looks for the peak pressure
    [pmax,J]=max(data2);
    if pmax>ppeak
        ppeak=pmax;
        tpeak=(N(I)+J)/fs;
    end
 
    % now we have filtered data for our segment.  We now select points with
    % dt spacing for plotting
    J=1:dI:length(data2);
    
    p2((Klast+1):(Klast+length(J)))=data2(J);  % extract the points spaced dI apart
    Klast=Klast+length(J);
    waitbar(I/Nseg,hw,wbstring);
end % for I=1:Nseg
delete(hw)
Lpdt=10*log10(p2)+94;       %convert the square pressure to Lp
Lppeak=10*log10(ppeak)+94;  % convert the peak square pressue to an Lp

t=(0:length(Lpdt)-1).*udt;  % generate a time vector with spacing of user defined dt

h=subplot(1,1,1,'Parent',handles.plots);
plot(t,Lpdt);
xlabel('t (sec)')
ylabel('Lp (dB)')
title('Sound Pressure Level');
        
if ~isscalar(handles.lp.plotscale)
    axis(handles.lp.plotscale);% if handles.leq.plotscale is not scalar, it must be an axis vector
end

set(h,'position',[0.13,0.2,0.775,0.71])
text(0.0,-0.13,sprintf('Weighting: %s        Meter Speed: %s',handles.weighting,handles.speed),'units','normalized');
wstring=handles.weighting(1);
if wstring=='F'
    wstring='';
end

pstring=sprintf('L_{%s%smax}=%2.1f dB at t=%1.1f sec.',wstring,handles.speed(1),Lppeak,tpeak);
text(0.0,-0.2,pstring,'units','normalized');

handles.lp.t=t;
handles.lp.Lpdt=Lpdt;
set(handles.analysis_save,'string','Save Lp Data');
handles.last.analysis='lp';
% save out changed handles
guidata(hObject, handles);

function analyze_leq(hObject,handles)
% function to compute the LEQ of the measurement.  
% This routine computes a 100 ms LEQ for statistic purpose and
% then combines them to get the desired LEQ size for plotting and saving
%
% you can change the 100ms to another number by changing handles.leq.mindt

filename=handles.meas.fullfile;
filesize=size(audioread(filename));
fsize=filesize(1);
fs=handles.meas.fs;

dt=handles.leq.mindt;
dtI=floor(fs*dt);  % compute the number of dt length samples in each dt length segment

if fsize<dtI
    h=warndlg('File is too short for analysis - LEQ analysis skipped');
    waitfor(h)
else
    % round off to an even number of dt sized segments
    Nseg=floor(fsize/dtI);
    if fsize/dtI==Nseg  % if we are exactly an even number of segments we need to subtract one from Nseg
        Nseg=Nseg-1; 
    end
    N=(1:dtI:fsize); 
end

msq=zeros(Nseg,1);  % create the array for the saving mean square computations
h=waitbar(0,'Please Wait ... Computing 0.1s LEQs');
for I=1:Nseg
    % read in the dt length segment, calibrate, apply weighting, and square
    data=handles.cal.scale.*audioread(filename,N(I:I+1));
    data2=(apply_weighting(data(:,1),fs,handles.weighting)).^2;  % make sure to take only channel 1
    
    % for each dt time period, find the msq
    msq(I)=mean(data2);
    waitbar(I/Nseg,h); 
end % I=1:Nseg
delete(h);

leqdata=10*log10(msq)+94;  % compute the LEQ in each 100 ms segment
overallmsq=mean(msq); % compute the overall mean square pressure for all segments
handles.leq.overall_leq=10*log10(overallmsq)+94;  % compute overall LEQ


%compute and save the statistics
handles.leq.Lmax=max(leqdata);  % find max LEQ
handles.leq.Lmin=min(leqdata); % find min LEQ
for I=1:9
    handles.leq.Ln(I)=percentile(leqdata,100-10*I);  % find the Ln values
end
handles.leq.Luser=percentile(leqdata,100-handles.leq.userpct); % find the user selected Ln value

% compute and save the noise dose
q=handles.leq.Dex/log10(2);  % normalized exchange rate
Tn=(480*60);   % reference time in seconds
I= leqdata>handles.leq.DLt;       % find the indices for when the data is above the criterion
D=sum(dt.*10.^((leqdata(I)-handles.leq.DLc)./q))/Tn;  % add up all the energy that counts
Ltwa=max(10*log10(D)+handles.leq.DLc,0);  % if Ltwa<0 dB , truncate at 0 dB        

%
% we've computed the  dt length LEQ  now combine these into the LEQ lengths
% requested by the user.
Udt=handles.leq.dt;  % load up the user selected dt

dI=Udt/dt;
NJ=floor(Nseg/dI);
msqJ=zeros(NJ,1);
for J=1:NJ
    msqJ(J)=mean(msq(((J-1)*dI+1):(J*dI)));
end
leqJ=10*log10(msqJ)+94;

handles.leq.leqdata=leqJ;  % save this LEQ data into our structure

t=(0:(NJ-1))*Udt;  % create the time vector for plotting as from 0 to end with a Udt spacing

subplot(2,1,1,'Parent',handles.plots);
cla reset
lineplot(t,leqJ); % plot the lines for each LEQ
axis tight
xlabel('t (sec)')
ylabel('LEQ (dB)')
title(sprintf('%2.1d s LEQ',Udt));

if ~isscalar(handles.leq.plotscale)
    % if handles.leq.plotscale is a vector then it must be an axis vector
    axis(handles.leq.plotscale)
end


subplot(2,1,2,'Parent',handles.plots);
cla reset
axis off

lc=-0.1;
mc=0.3;
hlc=0.0;
hrc=0.5;
rc=0.7;

wstring=handles.weighting(1);
if wstring=='F'
    wstring='';
end;
text(lc,1.0,'LEQ Statistics computed from 100 ms LEQ');
text(hlc,0.85,sprintf('Overall LEQ: %1.3G dB',handles.leq.overall_leq),'fontweight','bold');
text(hrc,0.85,sprintf('%s Weighting Applied',handles.weighting));

text(lc,.70,sprintf('L_{%smax}: %1.3G dB',wstring,handles.leq.Lmax));
text(mc,.70,sprintf('L_{%smin}: %1.3G dB',wstring,handles.leq.Lmin));
text(rc,.70,sprintf('L_{%s%s}: %1.3G dB',wstring,num2str(handles.leq.userpct),handles.leq.Luser));

text(lc,0.50,sprintf('L_{%s10}: %1.3G dB',wstring,handles.leq.Ln(1)));
text(mc,0.50,sprintf('L_{%s20}: %1.3G dB',wstring,handles.leq.Ln(2)));
text(rc,0.50,sprintf('L_{%s30}: %1.3G dB',wstring,handles.leq.Ln(3)));

text(lc,0.40,sprintf('L_{%s40}: %1.3G dB',wstring,handles.leq.Ln(4)));
text(mc,0.40,sprintf('L_{%s50}: %1.3G dB',wstring,handles.leq.Ln(5)));
text(rc,0.40,sprintf('L_{%s60}: %1.3G dB',wstring,handles.leq.Ln(6)));

text(lc,0.30,sprintf('L_{%s70}: %1.3G dB',wstring,handles.leq.Ln(7)));
text(mc,0.30,sprintf('L_{%s80}: %1.3G dB',wstring,handles.leq.Ln(8)));
text(rc,0.30,sprintf('L_{%s90}: %1.3G dB',wstring,handles.leq.Ln(9)));

if handles.weighting(1)=='A'
    text(lc,0.1,sprintf('Noise Dose Computations Using %s Criterion',handles.leq.Dtype));
    text(lc,-0.02,sprintf('Exchange Rate: %d dB',handles.leq.Dex));
    text(mc,-0.02,sprintf('Criterion Level: %d dB',handles.leq.DLc));
    text(rc,-0.02,sprintf('Threshold Level: %d dB',handles.leq.DLt));
    text(hlc,-0.18,sprintf('L_{TWA} = %1.3G dB',Ltwa),'fontweight','bold');
    text(hrc,-0.18,sprintf('Noise Dose D = %1.3f',D),'fontweight','bold');
    handles.leq.Ltwa=Ltwa;
    handles.leq.D=D;
else
    text(lc,0,'Cannot compute noise dose, A weighted LEQ is Required');
end


handles.leq.t=t;
handles.leq.leqdata=leqJ;
handles.last.analysis='leq';
set(handles.analysis_save,'string','Save LEQ Data');
% save out changed handles
guidata(hObject,handles);

function p=percentile(x,y)
% computes the yth percentile of the input vector x

N=length(x);
sx=sort(x);
idx=y*N/100;  % find the fractional index into the vector sx

n1=floor(idx);
if n1==idx
    p=sx(idx);
else
    n2=ceil(idx);
    x1=sx(n1);
    x2=sx(n2);
    p=interp1([n1,n2],[sx(n1),sx(n2)],idx);
end

% 
function lineplot(x,y)
% function to plot a set of connected lines for each data point.  Kind of
% like bar except we just print the top line rather than the full bar
%
% Copyright 2011 Ralph T Muehleisen
for I=1:(length(x)-1);
    line([x(I);x(I+1)],[y(I);y(I)]);  % draw the line across
    line([x(I+1);x(I+1)],[y(I);y(I+1)]);  % draw the line up/down
end


function analyze_psd(hObject,handles)
% function to compute the power spectral density using the Welch modified
% method.
%
%
filename=handles.meas.fullfile;
Nseg=handles.meas.Nseg;
N=handles.meas.N;
fs=handles.meas.fs;

Nfft=handles.psd.Nfft;
Noverlap=round(handles.psd.Nfft*handles.psd.overlap);
NPxx=Nfft/2+1;
Pxx=zeros(Nseg,NPxx);

dstring=sprintf('Please Wait ... Computing %d pt PSDs',Nfft);
h=waitbar(0,dstring);

switch handles.psd.window
    case {'Hanning'}
        win=window(@hann,Nfft);
    case 'Hamming'
        win=window(@hamming,Nfft);
    case 'Flattop'
        win=window(@flattopwin,Nfft);
end
for J=1:Nseg
    data=handles.cal.scale.*audioread(filename,N(J:J+1));
    [Pxx(J,:),fxx]=pwelch(data(:,1),win,Noverlap,Nfft,fs);
    H=acfilter(fxx',handles.weighting);
    Pxx(J,:)=Pxx(J,:).*abs(H).^2;

    waitbar(J/Nseg,h);
end
delete(h);

[Nr,~]=size(Pxx);
if Nr==1
    mPxx=Pxx;
else
    mPxx=mean(Pxx);  %compute the average of all the 
end;

Lpxx=10*log10(mPxx)+94;  % convert mean square pressure to a dB
Lpmax=max(Lpxx);
Lpd=max(Lpxx,Lpmax-60);


% clear display and get plot ready
subplot(1,1,1,'Parent',handles.plots)
cla reset
axis off
semilogx(fxx,Lpd);
axis tight
xlabel('Frequency (Hz)');
ylabel('Pxx (dB/Hz)');
tstring=sprintf('Power Spectral Density using %d pt FFT and %s window',Nfft,handles.psd.window);
title(tstring);

handles.psd.f=fxx;
handles.psd.LPxx=Lpxx;
handles.last.analysis='psd';
set(handles.analysis_save,'string','Save PSD Data');
guidata(hObject,handles);


function analyze_spec(hObject,handles)
% this routine implements the spectrogram analysis

% break time signal into slices and compute the PSD at each time slice
% then present in a color plot
dt=handles.spec.dt;
fs=handles.meas.fs;
dtI=floor(fs*dt);  % compute the number of samples in each segment
filename=handles.meas.fullfile;
filesize=size(audioread(filename));
fsize=filesize(1);

if fsize>dtI
    Nseg=floor(fsize/dtI);  % find the number of dt length segments
    N=1:dtI:fsize; % find the indices in the file for each segment
    N=[N,fsize]; % add in the last segment
else
    Nseg=1;
    N=[1,fsize];
end

%data=zeros(floor(fs*dt),1);
Nfft=handles.spec.Nfft;
Noverlap=round(Nfft*0.5);
NPxx=Nfft/2+1;
Pxx=zeros(NPxx,Nseg);

dstring=sprintf('Please Wait ... Computing %1.1f s X %d pt PSDs',dt,Nfft);
h=waitbar(0,dstring);
for I=1:Nseg
    % read in the dt length segment, calibrate, apply weighting, and square
    data=handles.cal.scale.*audioread(filename,N(I:I+1)); % read in the segment of the wave file
    [Pxx(:,I),fxx]=pwelch(data(:,1),Nfft,Noverlap,Nfft,fs);  % compute the PSD
    H=acfilter(fxx,handles.weighting); % compute the Weighting filter based on fxx
    Pxx(:,I)=Pxx(:,I).*abs(H).^2; % apply the weighting filter in the freq domain
    waitbar(I/Nseg,h);  % update waitbar
end % I=1:Nseg
delete(h);  % delete the waitbar

SdB=10*log10(Pxx)+94;
t=(0:(Nseg-1))*dt;

% if we have a value for the dnr, find the max SdB value and truncate
% everything below that
if handles.spec.dnr~=0
    Smax=max(max(SdB));
    SdB=max(SdB,Smax-handles.spec.dnr);
end


subplot(1,1,1,'parent',handles.plots);
cla reset
surf(t,fxx,SdB,'edgecolor','none'); 

if strcmp(handles.spec.view,'2d')
    axis tight; 
    view(0,90); % give a 2d view of spectrogram
else
    axis tight
    view(-30,46);  % give a 3D view of spectrogram
end
xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');
tstring=sprintf('Spectrogram with %d pt FFT and %1.1f s time slices',Nfft,dt);
title(tstring)
colorbar
handles.last.analysis='spec';
set(handles.analysis_save,'string','Save Spectrogram');
guidata(hObject,handles);



function analyze_bandfft(hObject,handles)
% function to find octave and 1/3 octave band LEQ using high resolution FFT
% This routine uses a rectangular window for FFTs which means there is some
% FFT sidelobe leakage which can increase the sound level of adjacent bands for
% a loud tone in a particular band
% This routine has been checked against a software routine that implements
% true buttworth digital filters that follow ANSI standards

% get file info from the measurement routine
filename=handles.meas.fullfile;
fs=handles.meas.fs;

Nfft=2^16;

% initialize the waitbar
h=waitbar(0,'Please Wait ... Computing Band Levels');

% break the wave file into segments that are multiples of NFFT, truncating
% the very end of the file. 
filesize=size(audioread(filename));  % get file size in samples
Nseg=floor(filesize(1)/Nfft); % break into Nfft size segments
N=1:Nfft:Nseg*Nfft;

SB=0;
fxx=linspace(0,fs/2,Nfft/2+1);
for J=1:Nseg-1
    data=handles.cal.scale.*audioread(filename,[N(J),N(J+1)-1]);  % read in the data for the segment and scale
    dataFFT=fft(data)/Nfft*sqrt(2);  % Take FFT and scale so the spectral amplitude is correct
    SP2=dataFFT.*conj(dataFFT);  % get square magnitude of FFT to get power spectrum
    SP=SP2(1:1+Nfft/2);  %extract the data from 0 to fs/2
    SB=SB+SP;
    waitbar(J/Nseg,h);
end

% Divide spectrum by Nseg to get an average and apply the weighting filter
H=acfilter(fxx,handles.weighting);
Pxx=SB(:)./Nseg.*abs(H(:)).^2;

delete(h)    
    
% generate the 1/3 octave band center frequencies
fcb=1000*2.^((-19:1:13)/3);
% select the fc that are no more than fs/2/sqrt(2) to ensure that our
% octave band edge will not exceed fs/2
J= fcb<(fs/2/sqrt(2));
fct=fcb(J);
% compute the upper and lower 1/3 octave band frequencies
flt=fct*2^(-1/6);
fut=fct*2^(1/6);

msqt=zeros(length(fct),1);
% 
for I=1:length(fct)
    J=(fxx>=flt(I))&(fxx<fut(I));
    msqt(I)=sum(Pxx(J));
end
Lpt=max(0,10*log10(msqt)+94);  % convert from a mean square to a dB with a min of 0 cB

Nfco=floor(length(fct)/3);

fco=fct(3*((1:Nfco)-1)+2);  %generate the octave band frequency vector from the 1/3 OB
msqo=0*fco;
for I=1:length(fco)
    msqo(I)=sum(msqt(3*(I-1)+1:1:3*(I-1)+3),1);  % add the 1/3 octave band mean square pressures
end
Lpo=max(0,10*log10(msqo)+94);  % convert from a mean square to a dB with a min of 0 cB


fctxt={'16','31.5','63','125','250','500','1k','2k','4k','8k','16k'};
if strcmp(handles.band.type,'oct')
    Lp=Lpo;
    fc=fco;
    tickmarks=1:11;
    tstring=sprintf('%s Weighted Octave Band LEQ computed using %d pt FFT',handles.weighting,Nfft);
else
    Lp=Lpt;
    fc=fct;
    tickmarks=2:3:32;
    tstring=sprintf('%s Weighted 1/3 Octave Band LEQ computed using %d pt FFT',handles.weighting,Nfft);
end

subplot(1,1,1,'Parent',handles.plots);
bar(1:length(fc),Lp);
set(gca,'xtick',tickmarks,'xticklabel',fctxt);
axis tight
xlabel('Frequency (Hz)')
ylabel('Lp (dB)');
title(tstring);
grid

% select and save the data we want to  keep
handles.band.Lp=Lp;
handles.band.fc=fc;
handles.band.NFFT=Nfft;
handles.band.Nseg=Nseg;  
handles.band.method='FFT';
handles.last.analysis='band';
set(handles.analysis_save,'string','Save Band Data')
guidata(hObject,handles);  % remember to include this so data is saved

function analyze_bandansi(hObject, handles)

filename=handles.meas.fullfile;
%fsize=handles.meas.fsize;
fs=handles.meas.fs;
filesize=size(audioread(filename));
fsize=filesize(1);

dt=0.1; % set system to read 100 ms buffers
dtI=floor(fs*dt);  % compute the number of samples in each dt length segment

Nseg=floor(fsize/dtI);  % break file into dt size segments and toss out the rest
if Nseg==0
    h=warndlg('File is too short for ANSI analysis - LEQ analysis skipped');
    waitfor(h)
    return
else
    N=(0:dtI:fsize)+1;
end

switch handles.band.type
    case 'oct'
        fcb=1000*2.^((-6:1:4));
        tstring=sprintf('%s Weighted Octave Band LEQ computed using ANSI Filters',handles.weighting);
        dstring='Please Wait Analyzing Data Using ANSI Octave Bands';
        tickmarks=1:length(fcb);
    otherwise
        fcb=1000*2.^((-19:1:13)/3);
        tstring=sprintf('%s Weighted 1/3 Octave Band LEQ computed using ANSI Filters',handles.weighting);
        dstring='Please Wait Analyzing Data Using ANSI 1/3 Octave Bands';
        tickmarks=2:3:length(fcb);
end

J= fcb<(fs/2/sqrt(2));
fc=fcb(J);

Hd=zeros(length(fc));
for I=1:length(fc)
   switch handles.band.type
        case 'oct'
            Hd(I)=nthoctdesign(fc(I),fs,1,5);  
        otherwise
            Hd(I)=nthoctdesign(fc(I),fs,3,5);     
   end % switch handles.band.type
   
   % set the filter to save the internal coefficients in between invocation so we
   % have a continuous filter outout even though we are submitting separate
   % segments
   Hd(I).PersistentMemory=1;  
end % for I=1:length(fc)



ms=zeros(Nseg,length(fc));
h=waitbar(0,dstring);
for J=1:Nseg
    data1=handles.cal.scale.*audioread(filename,[N(J),N(J+1)-1]);
    data1=apply_weighting(data1(:,1),handles.meas.fs,handles.weighting); % apply weighting filter and overwrite data1
    
    % do all the filters on that segment before reading the next segment
    for I=1:length(fc);
        data2=filter(Hd(I),data1);
        
        ms(J,I)=mean(abs(data2).^2);  % filter data and save the filter status so we 
    end
    waitbar(J/Nseg,h);
end
    delete(h) 
if ~isvector(ms)
    msq=mean(ms);
else
    msq=ms;
end


% convert mean square to abtn sound level
Lp=10*log10(msq)+94;

% set the minimum lpbtn=0 - ibtn.e. don't plot anything lower
I= Lp<0;
Lp(I)=0;



% plot the octave band levels
fctxt={'16','32','63','125','250','500','1k','2k','4k','8k','16k'};
bar(1:length(fc),Lp);
set(gca,'xtick',tickmarks,'xticklabel',fctxt);
axis tight
xlabel('Frequency (Hz)')
ylabel('Lp (dB)');
title(tstring);
grid



% select and save the data we want to  keep
handles.band.Lp=Lp;
handles.band.fc=fc;
handles.band.Nseg=Nseg;  
handles.last.analysis='band';
handles.band.method='ANSI';
set(handles.analysis_save,'string','Save Band Data')
% save out changed handles
guidata(hObject,handles);



function [Hd] = nthoctdesign(fc,fs,nth,N) 
% OCT3DSGN  Design of abtn one-nth-octave filter.
%    [hd] = nthoctdesign(fc,fs,nth,N)  designs an order N 1/nth-octave 
%    band digital filter with center frequency fc for sampling frequency fs. 
%    This returns a filter object to be used with filter (not the a,b
%    coefficients)
%
%    Requires the Signal Processing Toolbox. 

% Written by Ralph T Muehleisen
% The basic idea was inspired by a matlab routine octdesign by Christophe Couvreur,
% but this routine is general nth order and to uses z,p,k methods for better
% accuracy at low sample rates
% Last modification: 15-Mar-2011

f1 = fc/(2^(1/(2*nth))); 
f2 = fc*(2^(1/(2*nth))); 
d=fdesign.bandpass('N,F3dB1,F3dB2',2*N,f1,f2,fs);
Hd=design(d,'butter');






function analyze_ncrc(hObject,handles)
% function to find unweighted octave band LEQ using high resolution FFT and
% then compute NC and RC values
% the LEQ is a streamlined version of the one used for band analysis
filename=handles.meas.fullfile;
fs=handles.meas.fs;

Nfft=2^16;

% break the wave file into segments that are multiples of NFFT, truncating
% the very end of the file. 

filesize=size(audioread(filename));  % get file size in samples
Nseg=floor(filesize(1)/Nfft); % break into Nfft size segments
N=1:Nfft:Nseg*Nfft;

SB=0;
PD=1;
fxx=linspace(0,fs/2,Nfft/2+1);
h=waitbar(0,'Please Wait ... Computing Unweighted Octave Band Levels');
for J=1:Nseg-1
    
    data=handles.cal.scale.*audioread(filename,[N(J),N(J+1)-1]);
    dataFFT=fft(data)/Nfft*sqrt(2);  % scale FFT data properly so that we
    SP2=dataFFT.*conj(dataFFT);  % get square magnitude of FFT
    SP=SP2(1:1+Nfft/2);  % grab only half the
    SB=SB+SP;
    PD=PD+Nfft;
    waitbar(J/Nseg,h);
end
delete(h)

% Divide spectrum by Nseg to convert from sum to average over all segments
% note:  apply no weighting for NC/RC computations
Pxx=SB./Nseg;
    
% generate the octave band center frequencies from 16-8khz (or 4 kHz for
% a sampling rate of 22050 Hz 
if fs==22050
    fc=1000*2.^((-6:1:2));
    fctxt={'16','32','63','125','250','500','1k','2k','4k'};
    
else
    fc=1000*2.^((-6:1:3));
    fctxt={'16','32','63','125','250','500','1k','2k','4k','8k'};
end
% compute the upper and lower 1/3 octave band frequencies
fl=fc*2^(-1/2);
fu=fc*2^(1/2);

msq=zeros(length(fc),1);
for I=1:length(fc)
    J=(fxx>=fl(I))&(fxx<fu(I));
    msq(I)=sum(Pxx(J));
end
Lp=max(0,10*log10(msq)+94);  % convert from a mean square to a dB with a min of 0 cB



% define the NC curves in a matrix for plotting on top of octave band LEQ
nccurves=[     78	61	47	36	28	22	18	14	12	11; ...
                79	63	50	40	33	26	22	20	17	16; ...
                80	65	54	44	37	31	27	24	22	22; ...
                81	68	57	48	41	35	32	29	28	27; ...
                82	71	60	52	45	40	36	34	33	32; ...
                84	74	64	56	50	44	41	39	38	37; ...
                85	76	67	60	54	49	46	44	43	42; ...
                87	79	71	64	58	54	51	49	48	47; ...
                89	82	74	67	62	58	56	54	53	52; ...
                90	85	77	71	66	63	60	59	58	57; ...
                90.01	88	80	75	71	68	65	64	63	62; ...
                90.02	90	84	79	75	72	71	70	68	68];
ncval=[15;20;25;30;35;40;45;50;55;60;65;70];
[~,~,~]=NCRC(Lp);

tstring=sprintf('Unweighted Octave Band LEQ computed using %d pt FFT',Nfft);
tickmarks=1:length(fctxt);
h1=subplot(2,1,1,'Parent',handles.plots);

bar(1:length(fc),Lp,0.6);
set(gca,'xtick',tickmarks,'xticklabel',fctxt);
xlabel('Frequency (Hz)')
ylabel('Lp (dB)');
title(tstring);
grid;
hold on


for I=1:12
    text(length(fc)+0.5,nccurves(I,length(fc)),num2str(ncval(I)));
    line([length(fc);length(fc)+0.4],[nccurves(I,length(fc)),nccurves(I,length(fc))],'color','r');
    plot(1:length(fc),nccurves(I,1:length(fc)),'r');
end

hold off



h2=subplot(2,1,2,'Parent',handles.plots);
set(h1,'position',[0.13,0.31,0.775,0.615]);
set(h2,'position',[0.13,0.11,0.775,0.1]);
axis off
[ncstring,rcstring,~]=NCRC(Lp);

text(0,0.7,sprintf('NC Rating: %s',ncstring));
text(0,0.0,sprintf('RC Mark II Rating: %s',rcstring));

handles.ncrc.nc=ncstring;
handles.ncrc.rc=rcstring;
handles.ncrc.fc=fc;
handles.ncrc.Lp=Lp;
handles.last.analysis='ncrc';
set(handles.analysis_save,'string','Save NC/RC')
guidata(hObject,handles);  % remember to include this so data is saved

function [NCstring,RCstring,NCcurve]=NCRC(Lp)
% function to compute the NC and RNC level accorind to ANSI S12.2-2008
% The RC Mark II Level according to ASHRAE HVAC Apps 2009
%


% Copyright 2011, Ralph Muehleisen


nccurves=[     78	61	47	36	28	22	18	14	12	11; ...
                79	63	50	40	33	26	22	20	17	16; ...
                80	65	54	44	37	31	27	24	22	22; ...
                81	68	57	48	41	35	32	29	28	27; ...
                82	71	60	52	45	40	36	34	33	32; ...
                84	74	64	56	50	44	41	39	38	37; ...
                85	76	67	60	54	49	46	44	43	42; ...
                87	79	71	64	58	54	51	49	48	47; ...
                89	82	74	67	62	58	56	54	53	52; ...
                90	85	77	71	66	63	60	59	58	57; ...
                90.01	88	80	75	71	68	65	64	63	62; ...
                90.02	90	84	79	75	72	71	70	68	68];

f=[16 32 63 125 250 500 1000 2000 4000 8000];
ncval=[15;20;25;30;35;40;45;50;55;60;65;70];
NCf=0;


% first check to see if spectra is below NC curve everywhere.  If so we
% cannot rate it
if all(Lp(:)'<nccurves(1,1:length(Lp)))  % if we are under all NC curves at all freq 
     NCstring='NC is Below NC 15'; % set NC string to denote we are too low
     NCcurve=nccurves(1,1:length(Lp)); % set NC curve to NC 15
elseif any(Lp(:)'>nccurves(12,1:length(Lp))) % if any values is over the NC 70 curve
      NCstring='NC is Above NC 70'; % set NC string to denote we are too high
      NCcurve=nccurves(12,1:length(Lp));  % set NC curve to NC 70 
else
    % get SIL value for the current spectra, round to nearest dB
    sil=round(mean(Lp(6:9)));
    NCsil(1:length(Lp))=round(interp1(ncval,nccurves(:,1:length(Lp)),sil));
       
    % look for any Lp exceeding the NC_sil curve.
    if any(Lp(:)>NCsil(:))
        
        % if we are above NC_sil at any freq, find the tangent NC at each frequency
        NCtan=0*Lp;  %preallocate NCtan
        for I=1:length(Lp)
            if Lp(I)<nccurves(1,I);  % look for Lp < min NC value, and set NC=0 if it is
                NCtan(I)=0;
            elseif Lp(I)>nccurves(12,I) % look for Lp > max NC value, and set NC=100 if it is
                NCtan(I)=100;
            else
                NCtan(I)=round(interp1(nccurves(:,I),ncval(:),Lp(I))); % for all others, find the tangent NC through table interpolation
            end
        end
        NC=max(NCtan);  % set the NC value to the max tangent NC curve % 
        NCf=f(NC==NCtan);  % get the frequency of the tangency too
        NCcurve(1:length(Lp))=round(interp1(ncval,nccurves(:,1:length(Lp)),NC));  % get the NC curve through interpolation of the tables
        
    else  % if no Lp>NCsil assign the sil value as the NC value and the NC curve is the NC_sil curve
        NC=sil;
        NCcurve=NCsil;
    end
    if NCf==0
        NCstring=sprintf('NC %d',NC);
    else
        NCstring=sprintf('NC %d (%d)',NC,NCf);
        
    end
end

% now find the RC Mark II
psil=round(mean(Lp(6:8)));

RCcurve=psil+25-5*(1:length(Lp));
RCcurve(1)=RCcurve(2);
dev=Lp(:)-RCcurve(:);
LF=10*log10(mean(10.^(0.1*dev(1:3))));
MF=10*log10(mean(10.^(0.1*dev(4:6))));
HF=10*log10(mean(10.^(0.1*dev(7:9))));
D=[LF,MF,HF];
QAI=max(D)-min(D);

if Lp(1)>75 || Lp(2) > 75 || Lp(3)>80
    RCmod='HVa';
elseif Lp(1)>65 || Lp(2) > 65 || Lp(3)>75
    RCmod='HVb';
elseif QAI<5
    RCmod='N';
else
    if max(D)==LF
        RCmod='LF';
    elseif max(D)==MF
        RCmod='MF';
    elseif max(D)==HF
        RCmod='HF';
    end
end
RCstring=sprintf('RC %d (%s) with QAI: %1.1f',psil,RCmod,QAI);


function plotdisplay(handles,dstring,fontsize,fontstyle)
% get string to display.
%
% use default font 
if nargin<4
    fontstyle='default';
end

% use default size of 16
if nargin<3
    fontsize=16;
end
% get rid of the colorbar if there is one
colorbar('delete')
% get back to one axes if we had two and clear it
subplot(1,1,1,'Parent',handles.plots)
cla reset
axis off
text(-0.1,0.5,dstring,'fontsize',fontsize,'units','normalized', ...
    'horizontalalignment','left','fontname',fontstyle,'interpreter','none');
pause(0.25)


% --------------------------------------------------------------------
function xlpmenu_Callback(~, ~, ~)
% hObject    handle to xlpmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xlpdt_Callback(~, ~, ~)
% hObject    handle to xlpdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xlpscalemnu_Callback(~, ~, ~)
% hObject    handle to xlpscalemnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [H]=acfilter(f,weight)
% H = acfilter(fast,weight) designs a digital A- or C-Weighted filter
% evaluated at the frequencies fast.
%
% inputs are fast = a vector of frequencies at which to evaluate
% and weight, the weighting to apply, 'A' or 'C' are valid entries
% if no weight is given, a flat weighting is assumed and H=1;
%
% Reference ANSI S1.42

%    sampling frequency Fs. Usage: Y = FILTER(B,A,X). 
% A weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole at 108
% Hz and 738 Hz according to S1.42.
% C weighting has only the poles at 20 Hz and 12.2 kHz

% Copyright 2011 by Ralph T. Muehleisen

if nargin<2
    weight='A';
end

f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;


% generate the numerator and denominator polynomials for A and C 
switch upper(weight)
    case 'C'
        
        C1000=  0.0619;
        NUMs = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
        DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
        H=freqs(NUMs,DENs,2*pi*f);

        
    case 'A'
        % if not C, assume A
        A1000 = 1.9997;
        NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
        DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
        DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]);
        H=freqs(NUMs,DENs,2*pi*f);
        
    otherwise
        H=ones(size(f));
end

function [NUMs,DENs]=acfilter2(weight)
% H = acfilter(fast,weight) designs a digital A- or C-Weighted filter
% evaluated at the frequencies fast.
%
% inputs are fast = a vector of frequencies at which to evaluate
% and weight, the weighting to apply, 'A' or 'C' are valid entries
% if no weight is given, a flat weighting is assumed and H=1;
%
% Reference ANSI S1.42

%    sampling frequency Fs. Usage: Y = FILTER(B,A,X). 
% A weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole at 108
% Hz and 738 Hz according to S1.42.
% C weighting has only the poles at 20 Hz and 12.2 kHz

% Copyright 2011 by Ralph T. Muehleisen

if nargin<1
    weight='A';
end

f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;


% generate the numerator and denominator polynomials for A and C 
switch upper(weight)
    case 'C'
        
        C1000=  0.0619;
        NUMs = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
        DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
        %H=freqs(NUMs,DENs,2*pi*f);

        
    case 'A'
        % if not C, assume A
        A1000 = 1.9997;
        NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
        DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
        DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]);
        %H=freqs(NUMs,DENs,2*pi*f);
        
    otherwise
        %H=ones(size(f));
        NUMs=1;
        DENs=1;
end


% --- Executes on button press in saveplots.
function saveplots_Callback(~, ~, handles)
% hObject    handle to saveplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
h=figure;
copyobj(get(handles.plots,'Children'),h);
set(h,'menubar','figure');
 exportsetupdlg(h);



% --------------------------------------------------------------------
function xlpdt100ms_Callback(hObject, ~, handles)
% hObject    handle to xlpdt100ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.lp.dt=0.1;
set(handles.lp.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.lp.dth=hObject;  % update the selection handle
% save out changed handles
guidata(hObject,handles);

% --------------------------------------------------------------------
function xlpdt1s_Callback(hObject, ~, handles)
% hObject    handle to xlpdt1s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.lp.dt=1;
set(handles.lp.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.lp.dth=hObject;  % update the selection handle
guidata(hObject,handles);
% --------------------------------------------------------------------
function xlpdt10s_Callback(hObject, ~, handles)
% hObject    handle to xlpdt10s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.lp.dt=10;
set(handles.lp.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.lp.dth=hObject;  % update the selection handle
% save out changed handles
guidata(hObject,handles);

% --------------------------------------------------------------------
function xlpdt1min_Callback(hObject, ~, handles)
% hObject    handle to xlpdt1min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.lp.dt=60000;
set(handles.lp.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.lp.dth=hObject;  % update the selection handle
% save out changed handles
guidata(hObject,handles);


% --------------------------------------------------------------------
function xleqmenu_Callback(~, ~, ~)
% hObject    handle to xleqmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xleqdt_Callback(~, ~, ~)
% hObject    handle to xleqdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xleqscales_Callback(~, ~, ~)
% hObject    handle to xleqscales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xlpplotauto_Callback(hObject, ~, handles)
% hObject    handle to xlpplotauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.lp.plotscale=1;
set(handles.xlpplotmanual,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
guidata(hObject,handles);
% --------------------------------------------------------------------
function xlpplotmanual_Callback(hObject, ~, handles)
% hObject    handle to xlpplotmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isscalar(handles.lp.plotscale)
    alims=axis;
else
    alims=handles.lp.plotscale;
end
prompt={'x_{min}','x_{max}','y_{min}','y_{max}'};
dlgoptions.Resize='on';
dlgoptions.WindowStyle='normal';
dlgoptions.Interpreter='tex';
default=cellstr(num2str(alims(:)));
answer=inputdlg(prompt,'Input Plot Limits',1,default,dlgoptions);
if ~isempty(answer)
    set(handles.xlpplotauto,'Checked','off');  % turn off the check on the current selection
    set(hObject,'checked','on');  % turn on the current checkmark
    handles.lp.plotscale=[str2num(answer{1}),str2num(answer{2}),str2num(answer{3}),str2num(answer{4})];
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function xleqauto_Callback(hObject, ~, handles)
% hObject    handle to xleqauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.plotscale=1;
set(handles.xleqmanual,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
guidata(hObject,handles);

% --------------------------------------------------------------------
function xleqmanual_Callback(hObject, ~, handles)
% hObject    handle to xleqmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isscalar(handles.leq.plotscale)
    alims=axis;
else
    alims=handles.leq.plotscale;
end
prompt={'x_{min}','x_{max}','y_{min}','y_{max}'};
default=cellstr(num2str(alims(:)));
dlgoptions.Resize='on';
dlgoptions.WindowStyle='normal';
dlgoptions.Interpreter='tex';
answer=inputdlg(prompt,'Input Plot Limits',1,default,dlgoptions);
if ~isempty(answer)
    set(handles.xleqauto,'Checked','off');  % turn off the check automatic
    set(hObject,'checked','on');  % turn on the check on manual
    handles.leq.plotscale=[str2num(answer{1}),str2num(answer{2}),str2num(answer{3}),str2num(answer{4})];
end

guidata(hObject,handles);

% --------------------------------------------------------------------
function xleq1s_Callback(hObject, ~, handles)
% hObject    handle to xleq1s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.dt=1;
set(handles.leq.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xleq10s_Callback(hObject, ~, handles)
% hObject    handle to xleq10s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.dt=10;
set(handles.leq.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xleq1min_Callback(hObject, ~, handles)
% hObject    handle to xleq1min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.dt=60;
set(handles.leq.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xleq15min_Callback(hObject, ~, handles)
% hObject    handle to xleq15min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.dt=900;
set(handles.leq.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xleq1hr_Callback(hObject, ~, handles)
% hObject    handle to xleq1hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.dt=3600;
set(handles.leq.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.dth=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xleq100ms_Callback(hObject, ~, handles)
% hObject    handle to xleq100ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.dt=0.1;
set(handles.leq.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.dth=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xlequser_Callback(hObject, ~, handles)
% hObject    handle to xlequser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer=inputdlg('Input Percentile for LEQ','Input Percentile',1,{num2str(handles.leq.userpct)});
[val,status]=str2num(answer{1}); 
if status % if the conversion is valid, reset the userpct value
  handles.leq.userpct=val;  
end;  % if the if fails we cancelled so do not change anything
guidata(hObject,handles);


% --------------------------------------------------------------------
function xpsdmenu_Callback(~, ~, ~)
% hObject    handle to xpsdmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xpsdNfftmnu_Callback(~, ~, ~)
% hObject    handle to xpsdNfftmnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xpsd256_Callback(hObject, ~, handles)
% hObject    handle to xpsd256 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.Nfft=256;
set(handles.psd.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsd512_Callback(hObject, ~, handles)
% hObject    handle to xpsd512 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.Nfft=512;
set(handles.psd.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsd1k_Callback(hObject, ~, handles)
% hObject    handle to xpsd1k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.Nfft=1024;
set(handles.psd.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsd2k_Callback(hObject, ~, handles)
% hObject    handle to xpsd2k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.Nfft=2048;
set(handles.psd.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsd4k_Callback(hObject, ~, handles)
% hObject    handle to xpsd4k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.Nfft=4096;
set(handles.psd.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsd8k_Callback(hObject, ~, handles)
% hObject    handle to xpsd8k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.Nfft=8192;
set(handles.psd.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xpsdoverlapmnu_Callback(~, ~, ~)
% hObject    handle to xpsdoverlapmnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xpsdwindowmnu_Callback(~, ~, ~)
% hObject    handle to xpsdwindowmnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xpsdhann_Callback(hObject, ~, handles)
% hObject    handle to xpsdhann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.psd.window='Hanning';
set(handles.psd.windowh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.windowh=hObject;  % update the selection handle
guidata(hObject,handles);
% --------------------------------------------------------------------
function xpsdhamm_Callback(hObject, ~, handles)
% hObject    handle to xpsdhamm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.window='Hamming';
set(handles.psd.windowh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.windowh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsdflat_Callback(hObject, ~, handles)
% hObject    handle to xpsdflat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.window='Flattop';
set(handles.psd.windowh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.windowh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsdover0_Callback(hObject, ~, handles)
% hObject    handle to xpsdover0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.overlap=0;  % set overlap to 0%
set(handles.psd.overh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.overh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsdover25_Callback(hObject, ~, handles)
% hObject    handle to xpsdover25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.overlap=0.25;  % set overlap to 25%
set(handles.psd.overh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.overh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsdover50_Callback(hObject, ~, handles)
% hObject    handle to xpsdover50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.overlap=0.55;  % set overlap to 50%
set(handles.psd.overh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.overh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsdover75_Callback(hObject, ~, handles)
% hObject    handle to xpsdover75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.overlap=0.75;  % set overlap to 75%
set(handles.psd.overh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.overh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xpsdover90_Callback(hObject, ~, handles)
% hObject    handle to xpsdover90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.psd.overlap=0.9;  % set overlap to 90%
set(handles.psd.overh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.psd.overh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xload_settings_Callback(~, ~, ~)
% hObject    handle to xload_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%this is the gui that will be closed once we load the new settings

theCurrentGUI = gcf; 
h1=waitbar(0,'Reloading VSLM Settings');
hgload('vslmsave')
close(h1);
close(theCurrentGUI);


% --------------------------------------------------------------------
function specgrammenu_Callback(~, ~, ~)
% hObject    handle to specgrammenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xfftsizemnu_Callback(~, ~, ~)
% hObject    handle to xfftsizemnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xdtmnu_Callback(~, ~, ~)
% hObject    handle to xdtmnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xspecdt100m_Callback(hObject, ~, handles)
% hObject    handle to xspecdt100m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.dt=0.1;
set(handles.spec.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspecdt1s_Callback(hObject, ~, handles)
% hObject    handle to xspecdt1s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.dt=1.0;
set(handles.spec.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspecdt10s_Callback(hObject, ~, handles)
% hObject    handle to xspecdt10s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.dt=10.0;
set(handles.spec.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.dth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspec256_Callback(hObject, ~, handles)
% hObject    handle to xspec256 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.Nfft=256;
set(handles.spec.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspec512_Callback(hObject, ~, handles)
% hObject    handle to xspec512 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.Nfft=512;
set(handles.spec.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspec1024_Callback(hObject, ~, handles)
% hObject    handle to xspec1024 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.Nfft=1024;
set(handles.spec.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspec2048_Callback(hObject, ~, handles)
% hObject    handle to xspec2048 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.Nfft=2048;
set(handles.spec.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspec128_Callback(hObject, ~, handles)
% hObject    handle to xspec128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.Nfft=128;
set(handles.spec.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspec4096_Callback(hObject, ~, handles)
% hObject    handle to xspec4096 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.Nfft=4096;
set(handles.spec.Nffth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.Nffth=hObject;  % update the selection handle
guidata(hObject,handles);



% --------------------------------------------------------------------
function xspecdnr30_Callback(hObject, ~, handles)
% hObject    handle to xspecdnr30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.dnr=30;
set(handles.spec.dnrh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.dnrh=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xspecdnr50_Callback(hObject, ~, handles)
% hObject    handle to xspecdnr50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.dnr=50;
set(handles.spec.dnrh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.dnrh=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xspecdnrauto_Callback(hObject, ~, handles)
% hObject    handle to xspecdnrauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.dnr=0.;
set(handles.spec.dnrh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.spec.dnrh=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xspecplotscale_Callback(~, ~, ~)
% hObject    handle to xspecplotscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xlpdt20ms_Callback(hObject, ~, handles)
% hObject    handle to xlpdt20ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.lp.dt=0.02;
set(handles.lp.dth,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.lp.dth=hObject;  % update the selection handle
% save out changed handles
guidata(hObject,handles);


% --------------------------------------------------------------------
function xbandmenu_Callback(~, ~, ~)
% hObject    handle to xbandmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xoct_Callback(hObject, ~, handles)
% hObject    handle to xoct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.band.type='oct';
set(handles.band.typeh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.band.typeh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xthird_Callback(hObject, ~, handles)
% hObject    handle to xthird (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.band.type='third';
set(handles.band.typeh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.band.typeh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xnoisedosemenu_Callback(~, ~, ~)
% hObject    handle to xnoisedosemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xniosh_Callback(hObject, ~, handles)
% hObject    handle to xniosh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.Dtype='NIOSH';
handles.leq.Dex=3; % set the noise doise exchange rate to OSHA value of 5 dBA
handles.leq.DLc=85; % set the noise dose criterion to OSHA value of 90 dBA
handles.leq.DLt=80; % set the noise dose threshold to OSHA value of 80 dBA
set(handles.leq.Dtypeh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.Dtypeh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xosha_Callback(hObject, ~, handles)
% hObject    handle to xosha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leq.Dtype='OSHA';
handles.leq.Dex=5; % set the noise doise exchange rate to OSHA value of 5 dBA
handles.leq.DLc=90; % set the noise dose criterion to OSHA value of 90 dBA
handles.leq.DLt=80; % set the noise dose threshold to OSHA value of 80 dBA
set(handles.leq.Dtypeh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.leq.Dtypeh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xnoisedoseuser_Callback(hObject, ~, handles)
% hObject    handle to xnoisedoseuser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


defaults=[handles.leq.Dex;handles.leq.DLc;handles.leq.DLt];
prompt={'Exchange Rate (dB):','Criterion (dB)','Threshold (dB)'};
defaultstr=cellstr(num2str(defaults));
answer=inputdlg(prompt,'Input Noise Dose Criterion',1,defaultstr);
if ~isempty(answer)
    handles.leq.Dtype='USER';
    set(handles.leq.Dtypeh,'Checked','off');  % turn off the check on the current selection
    set(hObject,'checked','on');  % turn on the current checkmark
    handles.leq.Dtypeh=hObject;  % update the selection handle
    handles.leq.Dex=str2double(answer{1});
    handles.leq.DLc=str2double(answer{2});
    handles.leq.DLt=str2double(answer{3});
end

guidata(hObject,handles);

% --- Executes on button press in xpsdbtn.
function xpsdbtn_Callback(~, ~, ~)
% hObject    handle to xpsdbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in xncrcbtn.
function xncrcbtn_Callback(~, ~, ~)
% hObject    handle to xncrcbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in xlpbtn.
function xlpbtn_Callback(~, ~, ~)
% hObject    handle to xlpbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in xleqbtn.
function xleqbtn_Callback(~, ~, ~)
% hObject    handle to xleqbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in xbandbtn.
function xbandbtn_Callback(~, ~, ~)
% hObject    handle to xbandbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xspecview_Callback(~, ~, ~)
% hObject    handle to xspecview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xspecview2d_Callback(hObject, ~, handles)
% hObject    handle to xspecview2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.view='2d';
set(handles.spec.viewh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.band.typeh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xspecview3d_Callback(hObject, ~, handles)
% hObject    handle to xspecview3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spec.view='3d';
set(handles.spec.viewh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.band.typeh=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xhelp_Callback(~, ~, ~)
% hObject    handle to xhelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xproducthelp_Callback(~, ~, ~)
% hObject    handle to xproducthelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ispc
    dos('vslm.chm');
else
    web vslm.htm;
end

% --------------------------------------------------------------------
function xabout_Callback(~, ~, ~)
% hObject    handle to xabout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message={'Virtual Sound Level Meter','Version 0.4.2','Apr 5, 2017', ...
    'Copyright 2017 Ralph T Muehleisen','vslm.info@gmail.com'};
msgbox(message,'About VSLM')


% --- Executes on button press in analysis_save.
function analysis_save_Callback(~, ~, handles)
% hObject    handle to analysis_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.last.analysis;
if strcmp(handles.last.analysis,'none')
       dstring={'Cannot Save Last Analysis','No Analysis Performed'};
        plotdisplay(handles,dstring);
else
    
    switch handles.last.analysis
        case 'spec'
            %cannot save the spectrogram in a text file so separate the
            %plot and allow user to save that
            h=figure;
            copyobj(get(handles.plots,'Children'),h);
            plotedit(h,'showtoolsmenu');
            exportsetupdlg(h);
        case 'leq'
            [filename,pathname]=uiputfile('vslmleq.txt','Save NC/RC Data');
             switch filename
                 case {0}
                     % User cancelled out
                 otherwise
                     fullfilename=[pathname,filename];
                     fid = fopen(fullfilename, 'w');
                     fprintf(fid,'File: %s\n',handles.meas.fullfile);
                     fprintf(fid,'Weighting: %s\n',handles.weighting);
                     fprintf(fid,'Overall LEQ: %1.3G dB\n',handles.leq.overall_leq);
                     fprintf(fid,'Lmax: %1.3G dB\n',handles.leq.Lmax);
                     fprintf(fid,'Lmin: %1.3G dB\n',handles.leq.Lmin);
                     fprintf(fid,'L%1.0f: %1.3G dB\n',handles.leq.userpct,handles.leq.Luser);
                     fprintf(fid,'L10: %1.3G dB\n',handles.leq.Ln(1));
                     fprintf(fid,'L20: %1.3G dB\n',handles.leq.Ln(2));
                     fprintf(fid,'L30: %1.3G dB\n',handles.leq.Ln(3));
                     fprintf(fid,'L40: %1.3G dB\n',handles.leq.Ln(4));
                     fprintf(fid,'L50: %1.3G dB\n',handles.leq.Ln(5));
                     fprintf(fid,'L60: %1.3G dB\n',handles.leq.Ln(6));
                     fprintf(fid,'L70: %1.3G dB\n',handles.leq.Ln(7));
                     fprintf(fid,'L80: %1.3G dB\n',handles.leq.Ln(8));
                     fprintf(fid,'L90: %1.3G dB\n',handles.leq.Ln(9));

                     if strcmp(handles.weighting,'A')
                         fprintf('Noise Dose Computations Using %s Criterion\n',handles.leq.Dtype);
                         fprintf('Exchange Rate: %d dB\n',handles.leq.Dex);
                         fprintf('Criterion Level: %d dB\n',handles.leq.DLc);
                         fprintf('Threshold Level: %d dB\n',handles.leq.DLt);
                         fprintf('LTWA: %3.2G dB\n',handles.leq.Ltwa);
                         fprintf('Noise Dose D: %3.3f',handles.leq.D);
                     else
                         fprintf('Noise Dose Not Computed, A weighted LEQ is Required\n');
                     end
                     fprintf('Integration Time: %1.3G dB\n',handles.leq.dt);
                     fprintf(fid,'t (s) \tLeq (dB)\n');
                     for I=1:length(handles.leq.t)
                         fprintf(fid,'%2.1f \t%2.1f\n',handles.leq.t(I),handles.leq.leqdata(I));
                     end
                     fclose(fid);
             end
        case 'band'
           %
            [filename,pathname]=uiputfile('vslmband.txt','Save Band Analysis to .txt');
            switch filename
                case {0}
                    % User cancelled out
                otherwise
                    fullfilename=[pathname,filename];
                    fid = fopen(fullfilename, 'w');
                    fprintf(fid,'File: %s\n',handles.meas.fullfile);
                    fprintf(fid,'Weighting: %s\n',handles.weighting);
                    fprintf(fid,'Type: %s\n',handles.band.method);
                    if handles.band.type(1)=='t'
                        fprintf(fid,'Resolution: 1/3 Octave\n');
                    else 
                        fprintf(fid,'Resolution: Octave\n');
                    end
                    fprintf(fid,'fc (Hz) \tLp (dB) \n');
                    for I=1:length(handles.band.fc)
                        fprintf(fid,'%2.1f \t%2.1f\n',handles.band.fc(I),handles.band.Lp(I));
                    end
                    fclose(fid);
            end 
        case 'ncrc'
            % save NC/RC data
             [filename,pathname]=uiputfile('vslmncrc.txt','Save NC/RC Data');
             switch filename
                 case {0}
                     % User cancelled out
                 otherwise
                     fullfilename=[pathname,filename];
                     fid = fopen(fullfilename, 'w');
                     fprintf(fid,'File: %s\n',handles.meas.fullfile);
                     fprintf(fid,'%s \n',handles.ncrc.nc);
                     fprintf(fid,'%s \n\n',handles.ncrc.rc);
                     fprintf(fid,'fc (Hz) \tLp (dB) \n');
                     for I=1:length(handles.ncrc.fc)
                         fprintf(fid,'%2.1f \t%2.1f \n',handles.ncrc.fc(I),handles.ncrc.Lp(I));
                     end
                     fclose(fid);
             end
 
        case 'psd'
            [filename,pathname]=uiputfile('vslmpsd.txt','Save PSD to .txt');
            switch filename
                case {0}
                    % User cancelled out
                otherwise
                    fullfilename=[pathname,filename];
                    fid=fopen(fullfilename,'w');
                    fprintf(fid,'File: %s\n',handles.meas.fullfile);
                    fprintf(fid,'Weighting: %s\n',handles.weighting);
                    fprintf(fid,'f (Hz) \tLpxx (dB)\n');
                    for I=1:length(handles.psd.f)
                        fprintf(fid,'%2.1f \t%2.1f\n',handles.psd.f(I),handles.psd.LPxx(I));
                    end
                    fclose(fid);
            end
            
        case 'lp'
            [filename,pathname]=uiputfile('vslmlp.txt','Save Plot File');
            switch filename
                case {0}
                    % User cancelled out
                otherwise
                    fullfilename=[pathname,filename];
                    fid = fopen(fullfilename, 'w');
                    fprintf(fid,'File: %s\n',handles.meas.fullfile);
                    fprintf(fid,'Weighting: %s\n',handles.weighting);
                    fprintf(fid,'Meter Speed: %s\n',handles.speed);
                    fprintf(fid,'t (s) \tLp (dB) \n');
                    for I=1:length(handles.lp.t)
                        fprintf(fid,'%3.2f \t%2.1f \n',handles.lp.t(I),handles.lp.Lpdt(I));
                    end
                    fclose(fid);
            end
    end
    
end

% --------------------------------------------------------------------
function xsave_settings_Callback(~, ~, ~)
% hObject    handle to xsave_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hgsave('vslmsave')


% --- Executes on key press with focus on Fast and none of its controls.
function Fast_KeyPressFcn(~, ~, ~)
% hObject    handle to Fast (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (impulse.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xbandmethod_Callback(~, ~, ~)
% hObject    handle to xbandmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xbandfft_Callback(hObject, ~, handles)
% hObject    handle to xbandfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.band.method='FFT';
set(handles.band.methodh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.band.methodh=hObject;  % update the selection handle
guidata(hObject,handles);

% --------------------------------------------------------------------
function xbandansi_Callback(hObject, ~, handles)
% hObject    handle to xbandansi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.band.method='ANSI';
set(handles.band.methodh,'Checked','off');  % turn off the check on the current selection
set(hObject,'checked','on');  % turn on the current checkmark
handles.band.methodh=hObject;  % update the selection handle
guidata(hObject,handles);


% --------------------------------------------------------------------
function xbandres_Callback(~, ~, ~)
% hObject    handle to xbandres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
