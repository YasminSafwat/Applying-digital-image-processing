function varargout = fig3two(varargin)
% FIG3TWO MATLAB code for fig3two.fig
%      FIG3TWO, by itself, creates a new FIG3TWO or raises the existing
%      singleton*.
%
%      H = FIG3TWO returns the handle to a new FIG3TWO or the handle to
%      the existing singleton*.
%
%      FIG3TWO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIG3TWO.M with the given input arguments.
%
%      FIG3TWO('Property','Value',...) creates a new FIG3TWO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fig3two_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fig3two_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fig3two

% Last Modified by GUIDE v2.5 19-May-2022 00:08:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fig3two_OpeningFcn, ...
                   'gui_OutputFcn',  @fig3two_OutputFcn, ...
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


% --- Executes just before fig3two is made visible.
function fig3two_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fig3two (see VARARGIN)
guidata(handles.axes3);
img=getappdata(0,'X');
axes(handles.axes3);
imshow(img,'Parent',handles.axes3);

% Choose default command line output for fig3two
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes fig3two wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fig3two_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Display_Histogram(handles)
img= getframe(handles.axes1); 
[X]=frame2im(img);
[rows,cols]=size(X);
counts1=zeros(1,256);
for i=1:rows
 for j=1:cols
    grayLevel=X(i,j);
    counts1(grayLevel+1)=counts1(grayLevel+1)+1;
 end
end
grayLevels1 = 0 : 255;
bar(handles.axes5,grayLevels1, counts1, 'BarWidth', 1, 'FaceColor', 'b');

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gray_image=getappdata(0,'gray_image');
axes(handles.axes1);
imshow(gray_image,'Parent',handles.axes1);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uiputfile({'*.jpg','JPG'},'Save Image As');
f=getframe(handles.axes1);
[x]=frame2im(f); 
imwrite(x,fullfile(path, file),'jpg');


% --- Executes on button press in upload_2_image.
function upload_2_image_Callback(hObject, eventdata, handles)
% hObject    handle to upload_2_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.*','')
originalImage=imread([pathname, filename]);
imshow(originalImage,'Parent',handles.axes4);






% --- Executes on button press in Subtraction.
function Subtraction_Callback(hObject, eventdata, handles)
% hObject    handle to Subtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
l=getframe(handles.axes3);
[img1]=frame2im(l);
d=getframe(handles.axes4);
[img2]=frame2im(d);
[m,n]=size(img1);
new=img1;
for i=1:m
    for j=1:n
 new(i,j)=img1(i,j)-img2(i,j);  
    end
end
imshow(new,'parent',handles.axes1);
Display_Histogram(handles);



% --- Executes on button press in Adittion.
function Adittion_Callback(hObject, eventdata, handles)
% hObject    handle to Adittion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
l=getframe(handles.axes3);
[img1]=frame2im(l);
d=getframe(handles.axes4);
[img2]=frame2im(d);
[m,n]=size(img1);
new=img1;
for i=1:m
    for j=1:n
 new(i,j)=img1(i,j)+img2(i,j);  
    end
end
imshow(new,'parent',handles.axes1);
Display_Histogram(handles);
