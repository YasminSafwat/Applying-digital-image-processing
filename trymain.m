function varargout = trymain(varargin)
% TRYMAIN MATLAB code for trymain.fig
%      TRYMAIN, by itself, creates a new TRYMAIN or raises the existing
%      singleton*.
%
%      H = TRYMAIN returns the handle to a new TRYMAIN or the handle to
%      the existing singleton*.
%
%      TRYMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRYMAIN.M with the given input arguments.
%
%      TRYMAIN('Property','Value',...) creates a new TRYMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trymain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trymain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trymain

% Last Modified by GUIDE v2.5 18-May-2022 23:48:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trymain_OpeningFcn, ...
                   'gui_OutputFcn',  @trymain_OutputFcn, ...
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

% --- Executes just before trymain is made visible.
function trymain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trymain (see VARARGIN)

% Choose default command line output for trymain
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trymain wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = trymain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function setGlobal(val)
global original;
original = val;
function r = getGlobal
global original
r = original;

function Display_Histogram(handles)
img= getframe(handles.axes3); 
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
bar(handles.axes4,grayLevels1, counts1, 'BarWidth', 1, 'FaceColor', 'b');

% --- Executes on button press in upload_image.
function upload_image_Callback(~, eventdata, handles)
% hObject    handle to upload_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.*','')
originalImage1=imread([pathname, filename]);
setGlobal(originalImage1)
gray_image1 =originalImage1;
imshow(gray_image1,'Parent',handles.axes1);
imshow(gray_image1,'Parent',handles.axes3);
setappdata(0,'gray_image1',gray_image1);
[rows,cols]=size(gray_image1);
counts1=zeros(1,256);
for i=1:rows
 for j=1:cols
    grayLevel1=gray_image1(i,j);
    counts1(grayLevel1+1)=counts1(grayLevel1+1)+1;
 end
end
grayLevels1 = 0 : 255;
bar(handles.axes5,grayLevels1, counts1, 'BarWidth', 1, 'FaceColor', 'b');
bar(handles.axes4,grayLevels1, counts1, 'BarWidth', 1, 'FaceColor', 'b');

% --- Executes on button press in contrast.
function contrast_Callback(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=getappdata(0,'gray_image1');
prompt = {'Enter Minimum Level:','Enter Maximum Level:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
r_min = min(X(:));
r_max = max(X(:)); 
s_min=str2double(answer(1));
s_max=str2double(answer(2));
y= (X-r_min).*((s_max - s_min)/(r_max-r_min))+s_min;
imshow(y,'parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in Brightness.
function Brightness_Callback(hObject, eventdata, handles)
% hObject    handle to Brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img);
prompt = {'Enter Number:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
L=str2double(answer(1));
AddImage = X+L;
imshow(AddImage,'Parent',handles.axes3);
Display_Histogram(handles);


% --- Executes on button press in Darkness.
function Darkness_Callback(hObject, eventdata, handles)
% hObject    handle to Darkness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3);
[X]=frame2im(img); 
prompt = {'Enter Number:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
L=str2double(answer(1));
SubtractImage = X-L;
imshow(SubtractImage,'Parent',handles.axes3);
Display_Histogram(handles);


% --- Executes on button press in Multiplication.
function Multiplication_Callback(hObject, eventdata, handles)
% hObject    handle to Multiplication (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img); 
prompt = {'Enter Number:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
L=str2double(answer(1));
scalingImage = X*L;
imshow(scalingImage,'Parent',handles.axes3);
Display_Histogram(handles);


% --- Executes on button press in division.
function division_Callback(hObject, eventdata, handles)
% hObject    handle to division (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img); 
prompt = {'Enter Number:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
L=str2double(answer(1));
scalingImage = X/L;
imshow(scalingImage,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in thresholding.
function thresholding_Callback(hObject, eventdata, handles)
% hObject    handle to thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[xx]=frame2im(img);
X=rgb2gray(xx);
s=X;
[m,n]=size(X);
k = 95;
for i=1:m
    for j=1:n
        if X(i,j)>= k
            s(i,j)=255;
        else 
            s(i,j)=0;  
   
        end
    end
end
imshow(s,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in gray_level_slicing1.
function gray_level_slicing1_Callback(hObject, eventdata, handles)
% hObject    handle to gray_level_slicing1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[img]=frame2im(img); 
X = rgb2gray(img);
newImage=X;
prompt = {'Enter Minimum Level:','Enter Maximum Level:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
pp=str2double(answer(1));
xx=str2double(answer(2));
[rows cols] = size(X);
for row_index=1:1:rows
    for col_index=1:1:cols
        if(X(row_index,col_index)>=pp && X(row_index,col_index)<=xx)
            newImage(row_index,col_index) = 255;
        else
             newImage(row_index,col_index) = 0;
        end
    end
end
imshow(newImage,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in gray_level_slicing2.
function gray_level_slicing2_Callback(hObject, eventdata, handles)
% hObject    handle to gray_level_slicing2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[img]=frame2im(img); 
X = rgb2gray(img);
newImage=X;
prompt = {'Enter Minimum Level:','Enter Maximum Level:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
pp=str2double(answer(1));
xx=str2double(answer(2));
[rows cols] = size(X);
for row_index=1:1:rows
    for col_index=1:1:cols
        if(X(row_index,col_index)>=pp && X(row_index,col_index)<=xx)
            newImage(row_index,col_index) = 255;
        else
             newImage(row_index,col_index) = X(row_index,col_index);
        end
    end
end
imshow(newImage,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in Bit_plane_slicing.
function Bit_plane_slicing_Callback(hObject, eventdata, handles)
% hObject    handle to Bit_plane_slicing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[XX]=frame2im(img);
X=rgb2gray(XX);
[rows ,cols] = size(X);
newImage = zeros(rows,cols,8);
prompt = {'Enter Number of Bits:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
L=str2double(answer(1));
figure;
for k=1:L
    for row_index=1:1:rows
        for col_index=1:1:cols
            newImage(row_index,col_index,k)=bitget(X(row_index,col_index),k);
        end
    end

    subplot(2, 5, k+1),
    imshow(newImage(:,:,k));
    title(['Image for bit number',num2str(k)]);

end


% --- Executes on button press in Log.
function Log_Callback(hObject, eventdata, handles)
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[x]=frame2im(img);
img=rgb2gray(x); 
im = im2double(img);
for i=1:size(im,1)
    for j=1:size(im,2)
        in = double(im(i,j)+1);
        L(i,j) = log2(in);
    end
end
logT = (L/(max(L(:))))*255;
imshow(logT,[],'Parent',handles.axes3);
Display_Histogram(handles);


% --- Executes on button press in inverse_Log.
function inverse_Log_Callback(hObject, eventdata, handles)
% hObject    handle to inverse_Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xx= getframe(handles.axes3); 
[x]=frame2im(xx);
img=rgb2gray(x); 
im = im2double(img);
for i=1:size(im,1)
    for j=1:size(im,2)
        in =im(i,j)+1;
        E(i,j) = exp(in);
    end
end
expT = (E/(max(E(:))))*255;
imshow(expT,[],'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in Identity.
function Identity_Callback(hObject, eventdata, handles)
% hObject    handle to Identity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img);
imshow(X,'Parent',handles.axes3);
Display_Histogram(handles);


% --- Executes on button press in Negative_image.
function Negative_image_Callback(hObject, eventdata, handles)
% hObject    handle to Negative_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[xx]=frame2im(img);
X=rgb2gray(xx);
L = 2 ^ 8;     
neg = (L - 1) - X; 
imshow(neg,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in Power.
function Power_Callback(hObject, eventdata, handles)
% hObject    handle to Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[xx]=frame2im(img);
X=rgb2gray(xx);
double_value = im2double(X);
prompt = {'Enter Power Value:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
c=1;
p=str2double(answer(1));
out1= c*(double_value.^p); 
axes(handles.axes3);
imshow((out1),[]);
Display_Histogram(handles);

% --- Executes on button press in Histogram_Equalization.

function Histogram_Equalization_Callback(hObject, eventdata, handles)
% hObject    handle to Histogram_Equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes1); 
GIm=frame2im(img);
numofpixels=size(GIm,1)*size(GIm,2);
HIm=uint8(zeros(size(GIm,1),size(GIm,2)));
freq=zeros(256,1);
probc=zeros(256,1);
output=zeros(256,1);
for i=1:size(GIm,1)
    for j=1:size(GIm,2)
        value=GIm(i,j);
        freq(value+1)=freq(value+1)+1;
    end
end
sum=0;
no_bins=255;
for i=1:size(freq)
   sum=sum+freq(i);
   probc(i)=sum/numofpixels;
   output(i)=round(probc(i)*no_bins);
end
for i=1:size(GIm,1)
    for j=1:size(GIm,2)
            HIm(i,j)=output(GIm(i,j)+1);
    end
end
imshow(HIm,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in second_derivative.
function second_derivative_Callback(hObject, eventdata, handles)
% hObject    handle to second_derivative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=getframe(handles.axes3);
[xx]=frame2im(f);
input_image=rgb2gray(xx);
A = input_image(:,:,1);
Original_image=A;
Filtered_Image=zeros(size(A));
F=[1 1 1;1 -8 1; 1 1 1];
A=double(A);
for k=1:size(A,1)-2
    for input_image=1:size(A,2)-2
        Filtered_Image(k,input_image)=sum(sum(F.*A(k:k+2,input_image:input_image+2)));

     end
end
Filtered_Image= uint8(Filtered_Image);
Deblurred_image=Original_image-2*Filtered_Image;
axes(handles.axes3);
imshow(Deblurred_image);


% --- Executes on button press in Sharping_first.
function Sharping_first_Callback(hObject, eventdata, handles)
% hObject    handle to Sharping_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img);
grayimage1=rgb2gray(X);
gray_image = double(grayimage1);
[rows,cols]=size(gray_image);
mask = [0,1,0;1,-4,1;0,1,0];
out = gray_image;
for i=2:rows-1
    for j=2:cols-1
     temp = mask.*gray_image(i-1:i+1,j-1:j+1);
     value = sum(temp(:));
     out(i, j)= value;
    end
end
out = uint8(out);
axes(handles.axes3);
imshow(out);
Display_Histogram(handles);



% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=getappdata(0,'gray_image1');
prompt = {'Enter Minimum Level:','Enter Maximum Level:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);

% --- Executes on button press in Median_Filter.
function Median_Filter_Callback(hObject, eventdata, handles)
% hObject    handle to Median_Filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img);
input_image=rgb2gray(X);
input=inputdlg({'Filter Size: '},'Median Filter',1,{'3'});
input_array=cell2mat(input);
order=str2double(input_array);
[m n]=size(input_image);
p=(order-1)/2;
q=(order-1)/2;
output_image=input_image;
for i=p+1:m-p
	for j=q+1:n-q
		mask_array=zeros(order);
		for k=-p: p
			for l=-q:q
				mask_array(p+1-k,q+1-l)=input_image(i-k,j-l);
			end
		end
		Array=mask_array( : );
		median_value=median(Array);
		output_image(i,j)=median_value;
	end
end
imshow(output_image,'Parent',handles.axes3);
Display_Histogram(handles);

function and_Masks_Callback(hObject, eventdata, handles)
% hObject    handle to and_Masks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x= getframe(handles.axes3); 
[img]=frame2im(x);
img=rgb2gray(img);
prompt = {'Enter Row start:','Enter Row end :','Enter column start:','Enter column end:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
row_s=str2double(answer(1));
row_e=str2double(answer(2));
col_s=str2double(answer(3));
col_e=str2double(answer(4));
mask = img;
[rows cols] = size(img);
mask(row_s:row_e,col_s:col_e,:)= 255;
mask(1:row_s,:,:)= 0;
mask(row_e:rows,:,:)= 0;
mask(:,col_e:cols,:)= 0;
mask(:,1:col_s,:)= 0;
new=img;
for row_index=1:1:rows
    for col_index=1:1:cols
      if mask(row_index,col_index)==255
         new(row_index,col_index)=img(row_index,col_index) ;
      else
          new(row_index,col_index)=mask(row_index,col_index)-img(row_index,col_index); 
      end
      
    end
end
imshow(new,'parent',handles.axes3);
Display_Histogram(handles);
figure
mask=imresize(mask,0.25);
imshow(mask);
title('AND mask');

function Or_Masks_Callback(hObject, eventdata, handles)
% hObject    handle to Or_Masks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x= getframe(handles.axes3); 
[img]=frame2im(x);
img=rgb2gray(img);
prompt = {'Enter Row start:','Enter Row end :','Enter column start:','Enter column end:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
row_s=str2double(answer(1));
row_e=str2double(answer(2));
col_s=str2double(answer(3));
col_e=str2double(answer(4));
mask = img;
[rows cols] = size(img);
mask(row_s:row_e,col_s:col_e,:)= 0;
mask(1:row_s,:,:)= 255;
mask(row_e:rows,:,:)= 255;
mask(:,col_e:cols,:)= 255;
mask(:,1:col_s,:)= 255;
new=img;
for row_index=1:1:rows
    for col_index=1:1:cols
 new(row_index,col_index)=img(row_index,col_index)+mask(row_index,col_index);  
    end
end
imshow(new,'parent',handles.axes3);
Display_Histogram(handles);
figure
mask=imresize(mask,0.25);
imshow(mask);
title('OR mask');

% --- Executes on button press in Subtraction.
function Subtraction_Callback(hObject, eventdata, handles)
% hObject    handle to Subtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=getGlobal;
guidata(handles.axes3,img);
setappdata(0,'X',img);
fig3two;

% --- Executes on button press in Addition.
function Addition_Callback(hObject, eventdata, handles)
% hObject    handle to Addition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = getGlobal ;
guidata(handles.axes3,img);
setappdata(0,'X',img);
fig3two;

% --- Executes on button press in Reset2.
function Reset2_Callback(hObject, eventdata, handles)
% hObject    handle to Reset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gray_image1=getappdata(0,'gray_image1');
imshow(gray_image1,'Parent',handles.axes3);
Display_Histogram(handles);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uiputfile({'*.jpg','JPG'},'Save Image As');
f=getframe(handles.axes3);
[x]=frame2im(f); 
imwrite(x,fullfile(path, file),'jpg');

% --- Executes on button press in Save_Output_Histogram.
function Save_Output_Histogram_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Output_Histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uiputfile({'*.jpg','JPG'},'Save Image As');
f=getframe(handles.axes4);
[x]=frame2im(f); 
imwrite(x,fullfile(path, file),'jpg');

% --- Executes on button press in Save_Input_Histogram.
function Save_Input_Histogram_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Input_Histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uiputfile({'*.jpg','JPG'},'Save Image As');
f=getframe(handles.axes5);
[x]=frame2im(f); 
imwrite(x,fullfile(path, file),'jpg');


function crop_Callback(hObject, eventdata, handles)
% hObject    handle to crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x= getframe(handles.axes3); 
[img]=frame2im(x);
prompt = {'Enter Row start:','Enter Row end :','Enter column start:','Enter column end:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
row_s=str2double(answer(1));
row_e=str2double(answer(2));
col_s=str2double(answer(3));
col_e=str2double(answer(4));
crop=img(row_s:row_e,col_s:col_e,:);
imshow(crop,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in avg.
function avg_Callback(hObject, eventdata, handles)
% hObject    handle to avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[xx]=frame2im(img);
an=rgb2gray(xx);
prompt = {'Enter Averaging Mask size:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
b=str2double(answer(1));
[m,n]=size(an);
z=ones(b);
[p,q]=size(z);
w=1:p;
x=round(median(w));
anz=zeros(m+2*(x-1),n+2*(x-1));
for i=x:(m+(x-1))
    for j=x:(n+(x-1))
        anz(i,j)=an(i-(x-1),j-(x-1));
    end
end
sum=0;
x=0;
y=0;
for i=1:m
    for j=1:n
        for k=1:p
            for l=1:q 
                sum= sum+anz(i+x,j+y)*z(k,l);
                y=y+1;
            end
            y=0;
            x=x+1;
        end
        x=0;
        ansi(i,j)=(1/(p*q))*(sum);
        sum=0;
    end
end
axes(handles.axes3);
imshow(uint8(ansi));
Display_Histogram(handles);

% --- Executes on button press in Routate.
function Routate_Callback(hObject, eventdata, handles)
% hObject    handle to Routate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=getframe(handles.axes1);
[x]=frame2im(f);
R=rot90(x(:,:,1),1);
imshow(R,'Parent',handles.axes1);
imshow(R,'Parent',handles.axes3);



% --- Executes on button press in Maxfilter.
function Maxfilter_Callback(hObject, eventdata, handles)
% hObject    handle to Maxfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img);
input_image=rgb2gray(X);
input=inputdlg({'Filter Size: '},'Maximum Filter',1,{'3'});
input_array=cell2mat(input);
order=str2double(input_array);
[m n]=size(input_image);
p=(order-1)/2;
q=(order-1)/2;
output_image=input_image;
for i=p+1:m-p
	for j=q+1:n-q
		mask_array=zeros(order);
		for k=-p: p
			for l=-q:q
				mask_array(p+1-k,q+1-l)=input_image(i-k,j-l);
			end
		end
		Array=mask_array( : );
		maximum_value=max(Array);
		output_image(i,j)=maximum_value;
	end
end
axes(handles.axes3);
imshow(output_image);
Display_Histogram(handles);

% --- Executes on button press in Min_filter.
function Min_filter_Callback(hObject, eventdata, handles)
% hObject    handle to Min_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[X]=frame2im(img);
input_image=rgb2gray(X);
input=inputdlg({'Filter Size: '},'Maximum Filter',1,{'3'});
input_array=cell2mat(input);
order=str2double(input_array);
[m n]=size(input_image);
p=(order-1)/2;
q=(order-1)/2;
output_image=input_image;
for i=p+1:m-p
	for j=q+1:n-q
		mask_array=zeros(order);
			for k=-p: p
				for l=-q:q
                    mask_array(p+1-k,q+1-l)=input_image(i-k,j-l);
			end
		end
		Array=mask_array( : );
		minimum_value=min(Array);
		output_image(i,j)=minimum_value;
	end
end
axes(handles.axes3);
imshow(output_image);
Display_Histogram(handles);

% --- Executes on button press in Reset_main.
function Reset_main_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gray_image1=getappdata(0,'gray_image1');
imshow(gray_image1,'Parent',handles.axes1);


% --- Executes on button press in RGB_To_Gray.
function RGB_To_Gray_Callback(hObject, eventdata, handles)
% hObject    handle to RGB_To_Gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[input_image]=frame2im(img);
gray_image1 = rgb2gray(input_image);
imshow(gray_image1,'Parent',handles.axes3);



% --- Executes on button press in low_pass_ideal.
function low_pass_ideal_Callback(hObject, eventdata, handles)
% hObject    handle to low_pass_ideal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[x]=frame2im(img);
input_image=rgb2gray(x);
[M, N] = size(input_image);
FT_img = fft2(double(input_image));
D0 = 15; 
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2+V.^2);
H = double(D <= D0);
G = H.*FT_img;
output_image = real(ifft2(double(G)));
axes(handles.axes3);
imshow(output_image, [ ]);
Display_Histogram(handles);


% --- Executes on button press in low_pass_BLPF.
function low_pass_BLPF_Callback(hObject, eventdata, handles)
% hObject    handle to low_pass_BLPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
img= getframe(handles.axes3); 
[x]=frame2im(img);
input_image=rgb2gray(x);
[M, N] = size(input_image);
FT_img = fft2(double(input_image));
D0 = 15; 
n=2*2;
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2+V.^2);
D = D./ D0;
H = 1./((1+D).^n);
G = H.*FT_img;
output_image = real(ifft2(double(G)));
axes(handles.axes3);
imshow(output_image, [ ]);
Display_Histogram(handles);


% --- Executes on button press in low_pass_GLPF.
function low_pass_GLPF_Callback(hObject, eventdata, handles)
% hObject    handle to low_pass_GLPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
img= getframe(handles.axes3); 
[x]=frame2im(img);
input_image=rgb2gray(x);
[M, N] = size(input_image);
FT_img = fft2(double(input_image));
D0 = 15; 
D0 = (D0^2)*2;
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2+V.^2);
D = -D.^2;
H = exp(D/D0);
G = H.*FT_img;
output_image = real(ifft2(double(G)));
axes(handles.axes3);
imshow(output_image, [ ]);
Display_Histogram(handles);


% --- Executes on button press in high_pass_ideal.
function high_pass_ideal_Callback(hObject, eventdata, handles)
% hObject    handle to high_pass_ideal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img= getframe(handles.axes3); 
[x]=frame2im(img);
input_image=rgb2gray(x);
[M, N] = size(input_image);
FT_img = fft2(double(input_image));
D0 = 10; 
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2+V.^2);
H = double(D > D0);
G = H.*FT_img;
output_image = real(ifft2(double(G)));
axes(handles.axes3);
 imshow(output_image, [ ]);
 Display_Histogram(handles);



% --- Executes on button press in high_pass_BHPF.
function high_pass_BHPF_Callback(hObject, eventdata, handles)
% hObject    handle to high_pass_BHPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% MATLAB Code | Butterworth High Pass Filter

img= getframe(handles.axes3); 
[x]=frame2im(img);
input_image=rgb2gray(x);
[M, N] = size(input_image);
FT_img = fft2(double(input_image));
n = 2; 
D0 = 10; 
u = 0:(M-1);
v = 0:(N-1);
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);
H = 1./(1 + (D0./D).^(2*n));
G = H.*FT_img;
output_image = real(ifft2(double(G)));
axes(handles.axes3);
 imshow(output_image, [ ]);
 Display_Histogram(handles);



% --- Executes on button press in high_pass_GHPF.
function high_pass_GHPF_Callback(hObject, eventdata, handles)
% hObject    handle to high_pass_GHPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img= getframe(handles.axes3); 
[x]=frame2im(img);
a=rgb2gray(x); 
[m n]=size(a);
f_transform=fft2(a);
f_shift=fftshift(f_transform);
p=m/2;
q=n/2;
d0=70;
for i=1:m
for j=1:n
distance=sqrt((i-p)^2+(j-q)^2);
low_filter(i,j)=1-exp(-(distance)^2/(2*(d0^2)));
end
end
filter_apply=f_shift.*low_filter;
image_orignal=ifftshift(filter_apply);
image_filter_apply=abs(ifft2(image_orignal));
axes(handles.axes3);
imshow(image_filter_apply,[])
Display_Histogram(handles);



% --- Executes on button press in subsampling.
function subsampling_Callback(hObject, eventdata, handles)
% hObject    handle to subsampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=getframe(handles.axes3);
[originalImage]=frame2im(f);
[rows cols matricesNo] = size(originalImage);
input=inputdlg({'SamplingFactor: '},'SamplingFactor',1,{'3'});
input_array=cell2mat(input);
SamplingFactor=str2double(input_array);
for metricesIndex=1:1:matricesNo
    resizedImage(:,:,metricesIndex) = subSampling(originalImage(:,:,metricesIndex),SamplingFactor);
end
axes(handles.axes3);
imshow(resizedImage);
imwrite(resizedImage,'resizedImage.png');
Display_Histogram(handles);



% --- Executes on button press in upSampling.
function upSampling_Callback(hObject, eventdata, handles)
% hObject    handle to upSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=getframe(handles.axes3);
[originalImage2]=frame2im(f);
[rows cols matricesNo] = size(originalImage2);
input=inputdlg({'SamplingFactor: '},'SamplingFactor',1,{'3'});
input_array=cell2mat(input);
SamplingFactor=str2double(input_array);
for metricesIndex=1:1:matricesNo
    resizedImage2(:,:,metricesIndex) = upSampling(originalImage2(:,:,metricesIndex),SamplingFactor);
end
figure;
imshow(resizedImage2);
Display_Histogram(handles);
imwrite(resizedImage2,'resizedImageSamplingUp.jpg');


% --- Executes on button press in Binary_Image.
function Binary_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Binary_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[X]=frame2im(image);
img=rgb2gray(X);
[r,c]=size(img);
for i=1:r;
    for j=1:c;
        f=img(i,j);
        if f<=127;
            img(i,j)=0;
        elseif f>127 && f<256;
            img(i,j)=255;
        end
    end
end
imshow(img,'Parent',handles.axes3);
Display_Histogram(handles);

% --- Executes on button press in Edge_Detection.
function Edge_Detection_Callback(hObject, eventdata, handles)
% hObject    handle to Edge_Detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=getframe(handles.axes1);
[originalImage]=frame2im(f);
gray_image1 = rgb2gray(originalImage);  
gray_image = double(gray_image1);
[rows,cols]=size(gray_image);
mask = [-1 0 1;-2 0 2;-1 0 1];
out = gray_image;
for i=2:rows-1
 for j=2:cols-1
     temp = mask.*gray_image(i-1:i+1,j-1:j+1);
     value = sum(temp(:));
     out(i, j)= value;
end
end
out = uint8(out);
axes(handles.axes3);
imshow(out);
Display_Histogram(handles);


% --- Executes on button press in gray_scale.
function gray_scale_Callback(hObject, eventdata, handles)
% hObject    handle to gray_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[X]=frame2im(image);
img=rgb2gray(X);
k = 8;
figure
while (k > 0)
 target_levels = 2^k;
 target_compr_factor = 256 / target_levels;
 reduced_image = uint8(floor(double(img)/256 * target_levels) * target_compr_factor);
 subplot(3, 3, k); 
 imshow(reduced_image, [0 255]);
 if (k==1)
      title('Black & White');
 else
      title(['Grey-level resolution 2^',num2str(k)]);
 end
 k = k - 1;
end


% --- Executes on button press in Midpoint.
function Midpoint_Callback(hObject, eventdata, handles)
% hObject    handle to Midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[xx]=frame2im(image);
input_image=rgb2gray(xx);
[m n]=size(input_image);
input=inputdlg({'Order of Filter: '},'Maximum Filter',1,{'3'});
input_array=cell2mat(input);
order=str2double(input_array);
p=(order-1)/2;
q=(order-1)/2;
output_image=input_image;
for i=p+1:m-p
	for j=q+1:n-q
		mask_array=zeros(order);
		for k=-p: p
			for l=-q: q
				mask_array(p+1-k,q+1-l)=input_image(i-k,j-l);
			end
		end
		Array=mask_array(: );
		middle_value=(min(Array)+max(Array))/2;
		output_image(i,j)=middle_value;
	end
end
axes(handles.axes3);
imshow(output_image);
Display_Histogram(handles);


% --- Executes on button press in geo_mean.
function geo_mean_Callback(hObject, eventdata, handles)
% hObject    handle to geo_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[xx]=frame2im(image);
input_image=rgb2gray(xx);
input=inputdlg({'Order of Filter: '},'Geometeric Mean Filter',1,{'1'});
input_array=cell2mat(input);
mask=floor(str2num(input_array(1))/2);
[m n]=size(input_image);
output_image=input_image;
input_image=double(input_image);
for i = 1:m
    for j = 1:n 
    array_number=0; 
    array_multiplied=1;
        for k1 = i-mask:i+mask
            for p1 = j-mask:j+mask
                if ((k1 > 0 && p1 > 0) && (k1 < m && p1 < n ))
                    array_number = array_number+1;
                    array_multiplied=array_multiplied*input_image(k1,p1);
                end
            end
        end
        output_image(i,j)=array_multiplied^(1/array_number);
    end 
end
axes(handles.axes3);
imshow(output_image);
Display_Histogram(handles);

% --- Executes on button press in contra_mean.
function contra_mean_Callback(hObject, eventdata, handles)
% hObject    handle to contra_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[xx]=frame2im(image);
input_image=rgb2gray(xx);
input=inputdlg({'Mask maskze of filter (1-3): ','R of filter: '},'Contraharmonic Filter',1,{'1','0'});
input_array_1={'0','0'};
input_array=cell2mat(input(1));
input_array_1=cell2mat(input(2));
mask=str2num(input_array(1));

Q=str2num(input_array_1(1));
if (numel(input_array_1)==2)
    R=str2num(input_array_1(2));
    Q=-R;
end
input_image=im2double(input_image);
% Order of the filter
[m n]=size(input_image);
for i = 1:m
    for j = 1:n
        denominator=0; numerator=0;
        for k1 = i-mask:i+mask
            for p1 = j-mask:j+mask
                if ((k1>0 && p1 >0) && (k1<m && p1<n))
                    denominator=denominator+(input_image(k1,p1)^Q); 
                    numerator=numerator+(input_image(k1,p1)^(Q+1));
                end
            end
        end
        output_image(i,j)=numerator/denominator;
    end
end
axes(handles.axes3);
imshow(output_image);
Display_Histogram(handles);

% --- Executes on button press in alpha_trim.
function alpha_trim_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[xx]=frame2im(image);
input_image=rgb2gray(xx);
input=inputdlg({'Order of Filter: ','Trim Size : '},'Alpha Trim Filter',1,{'3','2'});
input_array=cell2mat(input);
masksize=str2num(input_array(1));
D=str2num(input_array(2))*2;
[m n]=size(input_image);
output_image=input_image;
d = D/2;
p=(masksize-1)/2;
q=(masksize-1)/2;
for i=p+1:m-p
	for j=q+1:n-q
		mask_array=zeros(masksize);
        for k=-p: p
            for l=-q:q
                mask_array(p+1-k,q+1-l)=input_image(i-k,j-l);
            end
        end
        array=mask_array(:);
        array = sort(array);
        r = (size(array)-d);
        s1= array(d+1:r);
        output_image(i,j)= sum(s1)/(size(mask_array,1)*(size(mask_array,2))-D);
    end
end
axes(handles.axes3);
imshow(output_image);
Display_Histogram(handles);



% --- Executes on button press in Arth_mean.
function Arth_mean_Callback(hObject, eventdata, handles)
% hObject    handle to Arth_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[ima]=frame2im(image)
x2=rgb2gray(ima);
[m,n]=size(x2);
for i=1:m-3
    for j=1:n-3
        a=ima(i:i+3,j:j+3);
        v(i,j)=sum(sum(a));
    end
end
y=mat2gray(v);
axes(handles.axes3);
imshow(y);
Display_Histogram(handles);




% --- Executes on button press in Harmonic_mean.
function Harmonic_mean_Callback(hObject, eventdata, handles)
% hObject    handle to Harmonic_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[Im]=frame2im(image);
x=rgb2gray(Im);
f = im2double(x);
[m n]=size(f);
si=1;
for i = 1:m
    for j = 1:n
    con=0; s1=0;
        for k1 = i-si:i+si
            for p1 = j-si:j+si
                if ((k1 > 0 && p1 > 0) && (k1 < m && p1 < n))
                    con = con+1;
                    if f(k1,p1)==0
                        s1 = s1+0;
                    else
                         s1=s1+(1/f(k1,p1));
                    end
                end
            end
        end
            b1(i,j)=con/s1;
    end
end
axes(handles.axes3);
imshow(b1);
Display_Histogram(handles);


% --- Executes on button press in Bandreject_Butterworth.
function Bandreject_Butterworth_Callback(hObject, eventdata, handles)
% hObject    handle to Bandreject_Butterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[aa]=frame2im(image);
micro=rgb2gray(aa);
micro = double(micro);
[nx ny] = size(micro);
u = micro;
micro = uint8(u);
imwrite(micro,'grass5.jpg');
fftu = fft2(u,2*nx-1,2*ny-1);
fftu = fftshift(fftu);
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);
n = 4;
for i = 1:2*nx-1
for j =1:2*ny-1
dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
% Use Butterworth filter.
filter1(i,j)= 1/(1 + (dist/120)^(2*n));
filter2(i,j) = 1/(1 + (dist/30)^(2*n));
filter3(i,j)= 1.0 - filter2(i,j);
filter3(i,j) = filter1(i,j).*filter3(i,j);
 filter3(i,j) = 1-filter3(i,j);
end
end
fil_micro = fftu + filter3.*fftu;
fil_micro = ifftshift(fil_micro);
fil_micro = ifft2(fil_micro,2*nx-1,2*ny-1);
fil_micro = real(fil_micro(1:nx,1:ny));
% fil_micro = uint8(fil_micro);
axes(handles.axes3);
imshow(fil_micro,[])
Display_Histogram(handles);



% --- Executes on button press in PandPassButter.
function PandPassButter_Callback(hObject, eventdata, handles)
% hObject    handle to PandPassButter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[ima]=frame2im(image);
micro = rgb2gray(ima);
micro = double(micro);
[nx ny] = size(micro);
u = micro;
micro = uint8(u);
imwrite(micro,'grass5.jpg');
fftu = fft2(u,2*nx-1,2*ny-1);
fftu = fftshift(fftu);
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);
n = 4;
for i = 1:2*nx-1
for j =1:2*ny-1
dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
% Use Butterworth filter.
filter1(i,j)= 1/(1 + (dist/120)^(2*n));
filter2(i,j) = 1/(1 + (dist/30)^(2*n));
filter3(i,j)= 1.0 - filter2(i,j);
filter3(i,j) = filter1(i,j).*filter3(i,j);
end
end
fil_micro = fftu + filter3.*fftu;
fil_micro = ifftshift(fil_micro);
fil_micro = ifft2(fil_micro,2*nx-1,2*ny-1);
fil_micro = real(fil_micro(1:nx,1:ny));
% fil_micro = uint8(fil_micro);
axes(handles.axes3);
imshow(fil_micro,[])
Display_Histogram(handles);


% --- Executes on button press in BandPass_ideal.
function BandPass_ideal_Callback(hObject, eventdata, handles)
% hObject    handle to BandPass_ideal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[ima]=frame2im(image);
I = rgb2gray(ima);
ima = double(I);
imafft = fftshift(fft2(fftshift(ima)));
imafft2 = fft2(ima);
imafft3 = fftshift(imafft2);
s = size(ima);
ma=max(max((imafft))); 
maxr = 0.5*sqrt(s(1)^2+s(2)^2); 
cutoff1 = maxr*30;
cutoff2 = maxr*120;
c=1;
for i = 1 :s(1)
    for j = 1 : s(2)
        r = sqrt((i-1-s(1)/2)^2+(j-1-s(2)/2)^2);
        if( r < 30) z(i,j) = 0;
        else if( r > 120) z(i,j) = 0; 
            else z(i,j) =511;
            end
        end
    end
end
% Plots
imafft=imafft.*z/255;
ima_out = fftshift(ifft2(fftshift(imafft)));
ima_out =ima_out-ima; ima_out = 1-ima_out;
axes(handles.axes3);
imshow(ima_out,[]);
Display_Histogram(handles);

% --- Executes on button press in BandPass_gaussian.
function BandPass_gaussian_Callback(hObject, eventdata, handles)
% hObject    handle to BandPass_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image= getframe(handles.axes3); 
[ima]=frame2im(image);
ima = rgb2gray(ima);
filtered_image = gaussianbpf2(ima,30,120);
axes(handles.axes3);
imshow(filtered_image,[]);
Display_Histogram(handles);
