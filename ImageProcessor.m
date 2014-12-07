%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name		: Jayaweera W. J. A. I. U.
%Index No	: 100227D
%Project - CS4722 - Computer Vision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X - Original image
%Y - Copy of X on which operations are performed

function varargout = ImageProcessor(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ImageProcessor_OpeningFcn, ...
    'gui_OutputFcn',  @ImageProcessor_OutputFcn, ...
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
end

% --- Executes just before ImageProcessor is made visible.
function ImageProcessor_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ImageProcessor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImageProcessor wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%clear global variables
clear global X;
clear global Y;
clear global N;
clear global M;
clear global Bt_Level;
clear isSegmented;
clear isEdgeDetected;
clear isNRed;
end

% --- Outputs from this function are returned to the command line.
function varargout = ImageProcessor_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- OnClickListners for Buttons --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hObject    handle to bPushButton/<Tag of button> (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in bBrowse.
function bBrowse_Callback(hObject, eventdata, handles)

%define global variables
global X;
global Y;
global isSegmented;
global isEdgeDetected;
global isNRed;
global File_Name;
global Path_Name;
global Bt_Level;

[File_Name, Path_Name] = uigetfile('C:\Users\Public\Pictures\*.jpg;*.png;*.jpeg');
if File_Name == 0
    set(handles.tvPath,'String','Select an image using browse  button');
else
    set(handles.tvPath,'String',[Path_Name,File_Name]);
    X = imread([Path_Name,File_Name]);
    Y = X;
    Bt_Level = 50;
    axes(handles.ivOriginal);
    imshow(X);
    generateHistogram(X,handles.pOriginal);
    
   isSegmented = false;
   isEdgeDetected = false;
   isNRed = false;

    showImageCopy(Y,handles.ivCopy,handles.pCopy);

    elementHandles = [handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
    activateGUIElements(elementHandles);

    set(handles.sBrightness, 'Value', Bt_Level);
    
    cla(handles.ivOrgCompare);    
    cla(handles.ivCopyCompare);     
    set(handles.tvCompInv,'String','......'); 
    set(handles.tvOrgInv,'String','......');
    set(handles.tvComparePath,'String','Select image to compare with above image using a "Compare" button','HorizontalAlignment','Left');
    set(handles.tvResults,'String','......');
end
end


% --- Executes on button press in bReset.
function bReset_Callback(hObject, eventdata, handles)

%define global variables
global X;
global Y;
global isSegmented;
global isEdgeDetected;
global isNRed;
global Bt_Level;

if ~isempty(X)
    Y = X;
    Bt_Level = 50;
    
    isSegmented = false;
    isEdgeDetected = false;
    isNRed = false;

    showImageCopy(Y,handles.ivCopy,handles.pCopy);

    set(handles.sBrightness, 'Value', Bt_Level);
    
    cla(handles.ivOrgCompare);    
    cla(handles.ivCopyCompare);     
    set(handles.tvCompInv,'String','......'); 
    set(handles.tvOrgInv,'String','......');
    set(handles.tvComparePath,'String','Select image to compare with above image using a "Compare" button','HorizontalAlignment','Left');
    set(handles.tvResults,'String','......');
end
end

% --- Executes on button press in bResetBrightness.
function bResetBrightness_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Bt_Level;

if ~isempty(Y)
    Bt_Level = 50;

    showImageCopy(Y,handles.ivCopy,handles.pCopy);

    set(handles.sBrightness, 'Value', Bt_Level);
end
end

% --- Executes on button press in bSave.
function bSave_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global File_Name;
global Path_Name;

if ~isempty(Y)
    imwrite(set_Brightness(Y),[Path_Name,strcat('Modified_',File_Name)]);
end
end

% --- Executes on button press in bTranspose.
function bTranspose_Callback(hObject, eventdata, handles)

%define global variables
global Y;

if ~isempty(Y)
    [height, width, dim] = size(Y);
    Z = im2double(Y);
    A = zeros( width, height, dim);
    for  k=1:dim
        for i=1:height
            for j=1:width
                A(j,i,k) = Z(i,j,k);
            end
        end
    end

    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
    
    if(~strcmp(get(handles.tvOrgInv,'String'),'......'))
        C = invariantMoments(parseToGray(Y));
        set(handles.tvOrgInv,'String',C);
        
        if(~strcmp(get(handles.tvResults,'String'),'......'))
            set(handles.tvResults,'String','......');
        end
    end
end
end

% --- Executes on button press in bVFlip.
function bVFlip_Callback(hObject, eventdata, handles)

%define global variables
global Y;

if ~isempty(Y)
    %     disp('Older Image');
    %     disp(Y);
    [height, width, dim] = size(Y);
    Z = im2double(Y);
    A = zeros(height, width, dim);
    for  k=1:dim
        for i=1:height
            for j=1:width
                A(i,width+1-j,k) = Z(i,j,k);
            end
        end
    end

    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
    
    if(~strcmp(get(handles.tvOrgInv,'String'),'......'))
        C = invariantMoments(parseToGray(Y));
        set(handles.tvOrgInv,'String',C);
        
        if(~strcmp(get(handles.tvResults,'String'),'......'))
            set(handles.tvResults,'String','......');
        end
    end
end
end

% --- Executes on button press in bNegative.
function bNegative_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Bt_Level;

if ~isempty(Y)
    Z = set_Brightness(Y);
    [height, width, dim] = size(Z);
    A = zeros(height, width, dim);
    A = 255 - Z;

    Y = A;
    showImageCopyNoBtChange(Y,handles.ivCopy,handles.pCopy);

    % Remove brightness built into the image, then it is added through mask
    % with next operation
    amt = 50 - Bt_Level;
    Y = shiftBrightnessAmountToMask(Y, amt);
    
    if(~strcmp(get(handles.tvOrgInv,'String'),'......'))
        C = invariantMoments(parseToGray(Y));
        set(handles.tvOrgInv,'String',C);
        
        if(~strcmp(get(handles.tvResults,'String'),'......'))
            set(handles.tvResults,'String','......');
        end
    end
end
end

% --- Executes on button press in bContrastStretch.
function bContrastStretch_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Bt_Level;

set(handles.bContrastStretch,'Enable','off');

if ~isempty(Y)
    Z_temp = set_Brightness(Y);
      
    A = stretchContrast(Z_temp);
    
    Y = im2uint8(A);
    showImageCopyNoBtChange(Y,handles.ivCopy,handles.pCopy);

    % Remove brightness built into the image, then it is added through mask
    % with next operation
    amt = 50 - Bt_Level;
    Y = shiftBrightnessAmountToMask(Y, amt);
end

set(handles.bContrastStretch,'Enable','on');

end

% --- Executes on button press in bNormalize.
function bNormalize_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Bt_Level;

set(handles.bNormalize,'Enable','off');

if ~isempty(Y)
    Z = set_Brightness(Y);
    A = normalize(Z);
    
    Y = im2uint8(A);
    showImageCopyNoBtChange(Y,handles.ivCopy,handles.pCopy);

    % Remove brightness built into the image, then it is added through mask
    % with next operation
    amt = 50 - Bt_Level;
    Y = shiftBrightnessAmountToMask(Y, amt);
end

set(handles.bNormalize,'Enable','on');
end

% --- Executes on button press in bBSliceWBG.
function bBSliceWBG_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Y_temp;
global Bt_Level;

if ~isempty(Y)

    point = ginput(2);

    disp('Input Locations');
    disp(point);

    Y_temp = Y;
    Z = double(set_Brightness(Y));
    x1 = Z(round(point(1,2)),round(point(1,1)),1);
    x2 = Z(round(point(2,2)),round(point(2,1)),1);

    disp('Input Points');
    disp([x1, x2]);

    % mask images for colour intensity regions
    mask_1 = double(Z<=x1);
    mask_2 = double((Z>x1)&(Z<x2));
    mask_3 = double(Z>=x2);

    % contrast stretching in regions
    im1 = mask_1.*floor(Z);
    im2 = mask_2.*floor(255*Z);
    im3 = mask_3.*floor(Z);

    % concatance of output image
    A = cast(im1+im2+im3,'uint8');

    Y = im2uint8(A);
    showImageCopyNoBtChange(Y,handles.ivCopy,handles.pCopy);

    % Remove brightness built into the image, then it is added through mask
    % with next operation
    amt = 50 - Bt_Level;
    Y = shiftBrightnessAmountToMask(Y, amt);

    elementHandles = [handles.bSliceUndo];
    activateGUIElements(elementHandles);

end

end

% --- Executes on button press in bBSliceWOBG.
function bBSliceWOBG_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Y_temp;
global Bt_Level;

if ~isempty(Y)
    point = ginput(2);

    disp('Input Locations');
    disp(point);

    Y_temp = Y;
    Z = double(rgb2gray(set_Brightness(Y)));
    [height, width, dim] = size(Z);
    A = zeros(height, width, dim);

    x1 = Z(round(point(1,2)),round(point(1,1)),1);
    x2 = Z(round(point(2,2)),round(point(2,1)),1);

    disp('Input Points');
    disp([x1, x2]);

    mask = double((Z>x1)&(Z<x2));
    A = A + mask*255;

    Y = A;
    showImageCopyNoBtChange(Y,handles.ivCopy,handles.pCopy);

    % Remove brightness built into the image, then it is added through mask
    % with next operation
    amt = 50 - Bt_Level;
    Y = shiftBrightnessAmountToMask(Y, amt);

    elementHandles = [handles.bSliceUndo];
    activateGUIElements(elementHandles);
end

end

% --- Executes on button press in bSliceUndo.
function bSliceUndo_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Y_temp;

if ~isempty(Y_temp)
    Y = Y_temp;
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
end
set(handles.bSliceUndo,'Enable','off');
end

% --- Executes on button press in bNNSampling.
function bNNSampling_Callback(hObject, eventdata, handles)

%define global variables
global Y;

if ~isempty(Y)    
    sm_value = 2;   
    
    [height, width, dim] = size(Y);
    Z = im2double(Y);
    A = zeros(height/sm_value, width/sm_value, dim);
    
    for  k=1:dim
        for i=1:height/sm_value
            for j=1:width/sm_value
                A(i,j,k) = Z(i*sm_value,j*sm_value,k);
            end
        end
    end

    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
end
end

% --- Executes on button press in bInvMoments.
function bInvMoments_Callback(hObject, eventdata, handles)

%define global variables
global Y;

if ~isempty(Y)   
    
    G = parseToGray(Y);
        
    A = invariantMoments(G);
    set(handles.tvOrgInv,'String',A);
end  
end

% --- Executes on button press in bGausian.
function bGausian_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global isNRed;

if ~isempty(Y)     
    A = addGaussian(Y,1);
    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
    isNRed = true;
end

end

% --- Executes on button press in bSegmentationGT.
function bSegment_Callback(hObject, eventdata, handles)

%define global variables
global X;
global Y;
global isSegmented;
global isEdgeDetected;

if(isEdgeDetected || isSegmented)
    Y = X;
end

elementHandles = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bSliceUndo, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
set(elementHandles,'Enable','off');

if ~isempty(Y)    
    U = normalize(Y);    
    V = addGaussian(U,1); 
    W = normalize(V);
   
    [height, width] = size(W);
    A = zeros(height, width);
    
    opt = optThresholding(W);
    mask = double(W<opt);
    A = A + mask*255;

    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
    
    if(~strcmp(get(handles.tvOrgInv,'String'),'......'))
        C = invariantMoments(A);
        set(handles.tvOrgInv,'String',C);
        
        if(~strcmp(get(handles.tvResults,'String'),'......'))
            set(handles.tvResults,'String','......');
        end
    end

    isSegmented = true;
end

elementHandlesAct = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
activateGUIElements(elementHandlesAct);

end

% --- Executes on button press in bSegmentationLT.
function bSegmentationLT_Callback(hObject, eventdata, handles)

%define global variables
global X;
global Y;
global isSegmented;
global isEdgeDetected;

if(isEdgeDetected || isSegmented)
    Y = X;
end

elementHandles = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bSliceUndo, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
set(elementHandles,'Enable','off');

if ~isempty(Y)    
    U = normalize(Y);    
    V = addGaussian(U,1); 
    W = normalize(V);
    
    A = localThresholding(W);

    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
    
    if(~strcmp(get(handles.tvOrgInv,'String'),'......'))
        C = invariantMoments(A);
        set(handles.tvOrgInv,'String',C);
        
        if(~strcmp(get(handles.tvResults,'String'),'......'))
            set(handles.tvResults,'String','......');
        end
    end

    isSegmented = true;
end

elementHandlesAct = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
activateGUIElements(elementHandlesAct);

end

% --- Executes on button press in bDetectEdges.
function bCanny_Callback(hObject, eventdata, handles)

%define global variables
global X;
global Y;
global isEdgeDetected;
global isSegmented;

if(isEdgeDetected || isSegmented)
    Y = X;
end

elementHandles = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bSliceUndo, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
set(elementHandles,'Enable','off');

if ~isempty(Y)    
    W = parseToGray(Y);
    
    A = cannyEdge(W);
    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
        
    if(~strcmp(get(handles.tvOrgInv,'String'),'......'))
        C = invariantMoments(A);
        set(handles.tvOrgInv,'String',C);
        
        if(~strcmp(get(handles.tvResults,'String'),'......'))
            set(handles.tvResults,'String','......');
        end
    end
    
    isEdgeDetected = true;
end
 
elementHandlesAct = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
activateGUIElements(elementHandlesAct);

end


% --- Executes on button press in bCompareGT.
function bCompare_Callback(hObject, eventdata, handles)

elementHandles = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bSliceUndo, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
set(elementHandles,'Enable','off');

%define global variables
global X;
global Y;
global N;
global M;
global File_Name_Compare;
global Path_Name_Compare;
global isSegmented;

Y = X;
C1 = zeros(1, 7);
C2 = zeros(1, 7);
Results = zeros(1, 7);
resultsCount = 0;

[File_Name_Compare, Path_Name_Compare] = uigetfile('C:\Users\Public\Pictures\Sample Pictures\*.jpg');
if File_Name_Compare == 0
    set(handles.tvComparePath,'String','Select image to compare with above image using a "Compare" button','HorizontalAlignment','Left');
else
    set(handles.tvComparePath,'String',[Path_Name_Compare,File_Name_Compare],'HorizontalAlignment','Left');
    N = imread([Path_Name_Compare,File_Name_Compare]);
    M = N;
    axes(handles.ivOrgCompare);
    imshow(N);
end


h = waitbar(0, 'Image comparison is in progress...');
set(h, 'WindowStyle','modal', 'CloseRequestFcn','');

if ~isempty(M)
    P = normalize(M);
    Q = addGaussian(P,1);
    A = normalize(Q);
    
    [height, width] = size(A);
    B = zeros(height, width);    
        
    opt = optThresholding(A);
    mask = double(A<opt);
    B = B + mask*255;
    M = im2uint8(B);  
    axes(handles.ivCopyCompare);
    imshow(M);    
    C1 = invariantMoments(B); 
    set(handles.tvCompInv,'String',C1);
end

if ~isempty(Y)
    V = normalize(Y);
    U = addGaussian(V, 1);
    A = normalize(U);
    
    [height, width] = size(A);
    B = zeros(height, width);
    
    opt = optThresholding(A);
    mask = double(A<opt);
    B = B + mask*255;
    Y = im2uint8(B);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);        
    C2 = invariantMoments(B);
    set(handles.tvOrgInv,'String',C2);
    
    isSegmented = true;
end

if (~isempty(Y) && ~isempty(M))
    for  i=1:7
        result1 = abs((C1(1,i)-C2(1,i))/C1(1,i))*100;
        result2 = abs((C1(1,i)-C2(1,i))/C2(1,i))*100;
        if(result1<10 || result2<10)
            Results(1,i)=1;            
            resultsCount = resultsCount + 1;
        else
            Results(1,i)=0;
        end
        waitbar(i/7,h);
    end
end

if(resultsCount>0)
    if(resultsCount>6)          
        set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are similar.']);
    else        
        set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are somewhat similar.']);
    end
else        
    set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are not similar.']);
end

delete(h);

elementHandlesAct = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
activateGUIElements(elementHandlesAct);
end

% --- Executes on button press in bCompareLT.
function bCompareLT_Callback(hObject, eventdata, handles)

elementHandles = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bSliceUndo, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
set(elementHandles,'Enable','off');

%define global variables
global X;
global Y;
global N;
global M;
global File_Name_Compare;
global Path_Name_Compare;
global isSegmented;

Y = X;
C1 = zeros(1, 7);
C2 = zeros(1, 7);
Results = zeros(1, 7);
resultsCount = 0;

[File_Name_Compare, Path_Name_Compare] = uigetfile('C:\Users\Public\Pictures\Sample Pictures\*.jpg');
if File_Name_Compare == 0
    set(handles.tvComparePath,'String','Select image to compare with above image using a "Compare" button','HorizontalAlignment','Left');
else
    set(handles.tvComparePath,'String',[Path_Name_Compare,File_Name_Compare],'HorizontalAlignment','Left');
    N = imread([Path_Name_Compare,File_Name_Compare]);
    M = N;
    axes(handles.ivOrgCompare);
    imshow(N);
end


h = waitbar(0, 'Image comparison is in progress...');
set(h, 'WindowStyle','modal', 'CloseRequestFcn','');

if ~isempty(M)
    P = normalize(M);
    Q = addGaussian(P,1);
    A = normalize(Q); 
        
    B = localThresholding(A);
    
    M = im2uint8(B);  
    axes(handles.ivCopyCompare);
    imshow(M);     
    C1 = invariantMoments(B);
    set(handles.tvCompInv,'String',C1);
end

if ~isempty(Y)
    V = normalize(Y);
    U = addGaussian(V, 1);
    A = normalize(U);
    
    B = localThresholding(A);
    
    Y = im2uint8(B); 
    showImageCopy(Y,handles.ivCopy,handles.pCopy);       
    C2 = invariantMoments(B);
    set(handles.tvOrgInv,'String',C2);
    
    isSegmented = true;
end

if (~isempty(Y) && ~isempty(M))
    for  i=1:7
        result1 = abs((C1(1,i)-C2(1,i))/C1(1,i))*100;
        result2 = abs((C1(1,i)-C2(1,i))/C2(1,i))*100;
        if(result1<10 || result2<10)
            Results(1,i)=1;            
            resultsCount = resultsCount + 1;
        else
            Results(1,i)=0;
        end
        waitbar(i/7,h);
    end
end

if(resultsCount>0)
    if(resultsCount>6)          
        set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are similar.']);
    else        
        set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are somewhat similar.']);
    end
else        
    set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are not similar.']);
end

delete(h);

elementHandlesAct = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
activateGUIElements(elementHandlesAct);
end


% --- Executes on button press in bCompareEdges.
function bCompareEdge_Callback(hObject, eventdata, handles)

elementHandles = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bSliceUndo, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
set(elementHandles,'Enable','off');

%define global variables
global X;
global Y;
global N;
global M;
global File_Name_Compare;
global Path_Name_Compare;
global isEdgeDetected;

Y = X;
C1 = zeros(1, 7);
C2 = zeros(1, 7);
Results = zeros(1, 7);
resultsCount = 0;

[File_Name_Compare, Path_Name_Compare] = uigetfile('C:\Users\Public\Pictures\Sample Pictures\*.jpg');
if File_Name_Compare == 0
    set(handles.tvComparePath,'String','Select image to compare with above image using a "Compare" button','HorizontalAlignment','Left');
else
    set(handles.tvComparePath,'String',[Path_Name_Compare,File_Name_Compare],'HorizontalAlignment','Left');
    N = imread([Path_Name_Compare,File_Name_Compare]);
    M = N;
    axes(handles.ivOrgCompare);
    imshow(N);
end

if ~isempty(M)    
    Q = parseToGray(M);
    
    A = cannyEdge(Q);
    
    M = im2uint8(A);
    axes(handles.ivCopyCompare);
    imshow(M);   
    C1 = invariantMoments(A);
    set(handles.tvCompInv,'String',C1);
end

if ~isempty(Y)    
    W = parseToGray(Y);
    
    A = cannyEdge(W);
    
    Y = im2uint8(A);
    showImageCopy(Y,handles.ivCopy,handles.pCopy); 
    
    C2 = invariantMoments(A);
    set(handles.tvOrgInv,'String',C2);
    isEdgeDetected = true;
end

h = waitbar(0, 'Image comparison is in progress...');
set(h, 'WindowStyle','modal', 'CloseRequestFcn','');

if (~isempty(Y) && ~isempty(M))
    for  i=1:7
        result1 = abs((C1(1,i)-C2(1,i))/C1(1,i))*100;
        result2 = abs((C1(1,i)-C2(1,i))/C2(1,i))*100;
        if(result1<10 || result2<10)
            Results(1,i)=1;            
            resultsCount = resultsCount + 1;
        else
            Results(1,i)=0;
        end
        waitbar(i/7,h);
    end
end

if(resultsCount>0)
    if(resultsCount>6)          
        set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are similar.']);
    else        
        set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are somewhat similar.']);
    end
else        
    set(handles.tvResults,'String',[num2str(resultsCount) ' out of 7 invariant moments are matching. Images are not similar.']);
end

delete(h);

elementHandlesAct = [handles.bBrowse, handles.bReset, handles.bResetBrightness, handles.bSave, handles.bTranspose, handles.bVFlip, handles.bNegative, handles.bContrastStretch, handles.bNormalize, handles.bBSliceWBG, handles.bBSliceWOBG, handles.bNNSampling, handles.sBrightness, handles.bInvMoments, handles.bGausian, handles.bSegment, handles.bSegmentationLT, handles.bCanny, handles.bCompare, handles.bCompareLT, handles.bCompareEdge];
activateGUIElements(elementHandlesAct);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- OnSliderMovementListners for Sliders --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hObject    handle to sBrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider

% --- Executes on slider movement.
function sBrightness_Callback(hObject, eventdata, handles)

%define global variables
global Y;
global Bt_Level;

if ~isempty(Y)
    Bt_Level = get(handles.sBrightness, 'Value');
    disp('Slider');
    disp(Bt_Level);
    showImageCopy(Y,handles.ivCopy,handles.pCopy);
else
    Bt_Level = 50;
    set(handles.sBrightness, 'Value', Bt_Level);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- OnCreate Function for Sliders --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hObject    handle to sBrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function sBrightness_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- User Defined Custom Functions --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Generate grayscale histograms for given image.
function h = generateHistogram(img, dest)

h = zeros(1,256);
if ~isempty(img)
    m = parseToGray(img);    
    [height, width] = size(m);
    n = im2uint8(m);
    for i=1:height
        for j=1:width
            h(1,n(i,j)+1)= h(1,n(i,j)+1)+1;
        end
    end
else
end
axes(dest);
axis([0 256 0 max(h)]);
area(h);
end

% --- Normalizing functionality
function normedImg = normalize(img)
    mx = max(img(:));
    mn = min(img(:));

    normedImg = round((img-mn)*(255/(mx-mn)));
end

% --- Contrast stretching functionality
function stretchedImg = stretchContrast(img)
    V = double(img(:));
    Z = double(img);
    x1 = prctile(V,5);
    x2 = prctile(V,95);

    [y1 y2 x1 x2] = getStretchedValues(x1, x2);

    a1 = y1/x1;
    a2 = (y2-y1)/(x2-x1);
    a3 = (255-y2)/(255-x2);

    % mask images for colour intensity regions
    mask_1 = double(Z<=x1);
    mask_2 = double((Z>x1)&(Z<x2));
    mask_3 = double(Z>=x2);

    % contrast stretching in regions
    im1 = mask_1.*floor(a1*Z);
    im2 = mask_2.*floor(y1 + (a2*(Z-x1)));
    im3 = mask_3.*floor(y2 + (a3*(Z-x2)));

    % concatance of output image
    stretchedImg = cast(im1+im2+im3,'uint8');
end

% --- Generate stretched values for given values.
function [y1, y2, x1, x2] = getStretchedValues(x1, x2)

if x1 == 0
    x1 = 1;
    y1 = x1;
else
    y1 = x1-floor((x1-0)*0.9);
end

if x2 == 255
    x2 = 254;
    y2 = x2;
else
    y2 = x2+floor((255-x2)*0.9);
end
end

% --- Set brightness of the image
function btm = set_Brightness(m)

%define global variables
global Bt_Level;

if ~isempty(m)
    l = m;

    if Bt_Level == 50 || isempty(Bt_Level)
        btm = l;
    else if Bt_Level < 50
            bt = round((50 - Bt_Level)*255/100);
            btm = l - bt;
        else
            bt = round((Bt_Level - 50)*255/100);
            btm = l + bt;
        end
    end
end
end

% --- Remove brightness from image from which brightness merged into
% --- If current brightness > 50, amt < 0; otherwise amt > 0
function btr_temp = shiftBrightnessAmountToMask(m, amt)

if ~isempty(m)
    btr_temp = m + amt;
end
end

% --- Show modified image
function showImageCopy(m,dest,plot)
Z = set_Brightness(m);
axes(dest);
imshow(Z);
generateHistogram(Z,plot);
end

% --- Show modified image without adding brightness mask
function showImageCopyNoBtChange(m,dest,plot)
Z = m;
axes(dest);
imshow(Z);
generateHistogram(Z,plot);
end

% --- Activate all the GUI elements
function activateGUIElements(eH)

elementCount = size(eH);

for i = 1:elementCount(2)
    set(eH(i),'Enable','on');
end
end

% --- Convert a RGB image to grayscale
function grayImg = parseToGray(img)
[height, width, dim] = size(img);
img = im2double(img);
if(dim == 3)
    grayImg = 0.299*img(:,:,1) + 0.587*img(:,:,2) + 0.114*img(:,:,3);
else
    grayImg = img;
end
end

% Calculating 7 invariant moments of the image, img (ref.[5])
function im = invariantMoments(img)

A = im2bw(img);

[H, W] = size(A);
[x, y] = meshgrid(1:W, 1:H);
m = zeros(1,10);
e = zeros(1,7);
im = zeros(1,7);

x = x(:);
y = y(:);
A = A(:);

% Calculating image moments
m(1,1) = sum(A);              %m00
if (m(1,1) == 0)
   m(1,1) = eps;
end 
m(1,2)  = sum(x .* A);        %m10
m(1,3)  = sum(y .* A);        %m01
m(1,4)  = sum(x .* y .* A);   %m11
m(1,5)  = sum(x.^2 .* A);     %m20
m(1,6)  = sum(y.^2 .* A);     %m02
m(1,7)  = sum(x .* y.^2 .* A);%m12
m(1,8)  = sum(x.^2 .* y .* A);%m21
m(1,9)  = sum(x.^3 .* A);     %m30
m(1,10) = sum(y.^3 .* A);     %m03

x_bar = m(1,2) / m(1,1);
y_bar = m(1,3) / m(1,1);

% Calculating normalized central moments
e(1,1) = (m(1,4) - y_bar*m(1,2)) / m(1,1)^2;   %eta11
e(1,2) = (m(1,5) - x_bar*m(1,2)) / m(1,1)^2;   %eta20
e(1,3) = (m(1,6) - y_bar*m(1,3)) / m(1,1)^2;   %eta02
e(1,4) = (m(1,9) - 3 * x_bar * m(1,5) + 2 * x_bar^2 * m(1,2)) / m(1,1)^2.5;    %eta30
e(1,5) = (m(1,10) - 3 * y_bar * m(1,6) + 2 * y_bar^2 * m(1,3)) / m(1,1)^2.5;   %eta03
e(1,6) = (m(1,8) - 2 * x_bar * m(1,4) - y_bar * m(1,5) + 2 * x_bar^2 * m(1,3)) / m(1,1)^2.5;   %eta21
e(1,7) = (m(1,7) - 2 * y_bar * m(1,4) - x_bar * m(1,6) + 2 * y_bar^2 * m(1,2)) / m(1,1)^2.5;   %eta12

% Calculating invariant moments
im(1,1) = e(1,2) + e(1,3);                                  %phi1
im(1,2) = (e(1,2) - e(1,3))^2 + 4*e(1,1)^2;                 %phi2
im(1,3) = (e(1,4) - 3*e(1,7))^2 + (3*e(1,6) - e(1,5))^2;    %phi3
im(1,4) = (e(1,4) + e(1,7))^2 + (e(1,6) + e(1,5))^2;        %phi4
im(1,5) = (e(1,4) - 3*e(1,7)) * (e(1,4) + e(1,7)) *((e(1,4) + e(1,7))^2 - 3*(e(1,6) + e(1,5))^2 ) + (3*e(1,6) - e(1,5)) * (e(1,6) + e(1,5)) * ( 3*(e(1,4) + e(1,7))^2 - (e(1,6) + e(1,5))^2 );  %phi5
im(1,6) = (e(1,2) - e(1,3)) * ( (e(1,4) + e(1,7))^2 - (e(1,6) + e(1,5))^2 ) + 4 * e(1,1) * (e(1,4) + e(1,7)) * (e(1,6) + e(1,5));                                                               %phi6
im(1,7) = (3*e(1,6) - e(1,5)) * (e(1,4) + e(1,7)) * ((e(1,4) + e(1,7))^2 - 3*(e(1,6) + e(1,5))^2 ) + (3*e(1,7) - e(1,4)) * (e(1,6) + e(1,5)) * ( 3*(e(1,4) + e(1,7))^2 - (e(1,6) + e(1,5))^2 ); %phi7
end

% Reduce noise by convolving image, img, with the Gaussian mask
function nr = addGaussian(img, type)
grayImg = parseToGray(img);
if(type == 1)
    G = double([1 2 3 2 1; 2 7 11 7 2; 3 11 17 11 3; 2 7 11 7 2; 1 2 3 2 1]*(1/121));
else
    if(type == 2)
        G = double([2 4 5 4 2; 4 9 12 9 4; 5 12 15 12 5; 4 9 12 9 4; 2 4 5 4 2]*(1/115));
    else
        G = double([1 4 7 4 1; 4 16 26 16 4; 7 26 41 26 7; 4 16 26 16 4; 1 4 7 4 1]*(1/273));
    end
end
nr = conv2(grayImg,G,'same');
end

% Find the optimal threshold in order to do segmentation of the image, img (ref.[1]) 
function optimalThreshold = optThresholding(img)

wb = waitbar(0, 'Segmentation is in progress...');
set(wb, 'WindowStyle','modal', 'CloseRequestFcn','');

% Inorder to store the criterion function values
J = Inf * ones(255, 1);

% Computing the relative histogram
h = double(histc(img(:), 0:255)) / size(img(:),1);

% Consider each intensity value as a potential threshold
for T = 1:255

    % Split the hostogram at the threshold T
    h1 = h(1:T);
    h2 = h((T+1):end);

    % Calculate portions of pixels in the two classes induced by T
    P1 = sum(h1);
    P2 = sum(h2);

    % Continue when both portions are non-empty
    % In order to avoid the error of division by zero
    if ((P1 > 0) && (P2 > 0))

        % Calculate the mean of each of the two sections induced by 
        % threshold T
        m1 = sum(h1 .* (1:T)') / P1;
        m2 = sum(h2 .* (1:(256-T))') / P2;
        
        % Calculate the standard deviation of each of the two sections
        % induced by threshold T
        sd1 = sqrt(sum(h1 .* ((transpose(1:T) - m1) .^2) ) / P1);
        sd2 = sqrt(sum(h2 .* ((transpose(1:(256-T)) - m2) .^2) ) / P2);

        % T becomes a potential threshold only if both standard deviation
        % are greater than zero
        if (sd1 > 0) && (sd2 > 0)

            % Criterion function calculation according to ref.[3]
            J(T) = 1 + 2 * (P1 * log(sd1) + P2 * log(sd2)) - 2 * (P1 * log(P1) + P2 * log(P2));

        end
    end
    waitbar(T/255, wb);
end

% The threshold which induced the minimum value for the criterion function
% becomes the optimal threshold
[minCount,optimalThreshold] = min(J);
optimalThreshold = optimalThreshold - 0.5;

delete(wb);

end

% Segment the image using the local thresholding
function segmentedMat = localThresholding(img)

    wb = waitbar(0, 'Segmentation is in progress...');
    set(wb, 'WindowStyle','modal', 'CloseRequestFcn','');

    [h, w] = size(img);
    segmentedMat = zeros(h, w);

    % Calculate average intensity of each pixel's 250x250 neighborhood 
    avgMat=imfilter(img,fspecial('average',250),'replicate');

    % Use average values caculated as thresholds to decide whether the 
    % pixel is a background or a foreground pixel 
    for i=1:h
        for j=1:w
            if(img(i,j)<=(avgMat(i,j)-0.058))
                segmentedMat(i,j)=1;
            else
                segmentedMat(i,j)=0; 
            end
        end        
        waitbar(i/h, wb);
    end

    delete(wb);
end

% Detect and mark edges in the image using Canny edge detection algorithm
% (ref.[4])
function edgeMap = cannyEdge(img)

    wb = waitbar(0, 'Edge detection is in progress...');
    set(wb, 'WindowStyle','modal', 'CloseRequestFcn','');
        
	% Convolve with a Gaussian mask to filter out noises
    img = addGaussian(img,3);
    
    % Sobel masks to calculate gradients in X and Y directions 
    SobelX = [-1,0,+1;-2,0,+2;-1,0,+1];
    SobelY = [-1,-2,-1;0,0,0;+1,+2,+1];
    
    % Convolving the image with the Sobel masks
    Gx = conv2(img,SobelX,'same');
    Gy = conv2(img,SobelY,'same'); 
    
	% Obtaining gradient matrix by calculating the magnitudes of gradient
	% vectors resulted from gradients in X and Y directions
	Gm = sqrt(Gx.^2+Gy.^2);

	% Defining threshold for non-maximum suppression
	Gm_max = max(max(Gm));
	Gm_min = min(min(Gm));
	T_nonMax = 0.07*(Gm_max-Gm_min) + Gm_min;
    T_Max = 0.1*(Gm_max-Gm_min) + Gm_min;
    
    % Performing non-maximum suppression 
	Gm_nonMax = max(Gm,T_nonMax.*ones(size(Gm)));

	% Interpolating to find the pixels where the magnitudes of gradient are 
    % locally maximum
	[h,w] = size(Gm_nonMax);
    edgeMap = zeros(h,w);
    for i = 2:h-1,
        for j = 2:w-1,
            if (Gm_nonMax(i,j) > T_nonMax)
                X=[-1,0,+1;-1,0,+1;-1,0,+1];
                Y=[-1,-1,-1;0,0,0;+1,+1,+1];
                T=[Gm_nonMax(i-1,j-1),Gm_nonMax(i-1,j),Gm_nonMax(i-1,j+1);
                   Gm_nonMax(i,j-1),Gm_nonMax(i,j),Gm_nonMax(i,j+1);
                   Gm_nonMax(i+1,j-1),Gm_nonMax(i+1,j),Gm_nonMax(i+1,j+1)];
                XI=[Gx(i,j)/Gm(i,j), -Gx(i,j)/Gm(i,j)];
                YI=[Gy(i,j)/Gm(i,j), -Gy(i,j)/Gm(i,j)];
                TI=interp2(X,Y,T,XI,YI);
                if (Gm_nonMax(i,j) >= TI(1) && Gm_nonMax(i,j) >= TI(2))
                    if(Gm_nonMax(i,j) > T_Max)                        
                        edgeMap(i,j) = 255;
                    else                        
                        edgeMap(i,j)=T_nonMax;
                    end
                end
            end
        end
        waitbar(i/h, wb);
    end 
    
    % Performing hysteresis
    looping = true;
    while(looping)
        looping = false;
        for i = 2:h-1,
            for j = 2:w-1,
                if (edgeMap(i,j) < T_Max && edgeMap(i,j) > 0)
                    if (edgeMap(i-1,j-1) > T_Max || edgeMap(i-1,j) > T_Max||edgeMap(i-1,j+1) > T_Max|| edgeMap(i,j-1) > T_Max||edgeMap(i,j) > T_Max||edgeMap(i,j+1) > T_Max|| edgeMap(i+1,j-1) > T_Max||edgeMap(i+1,j) > T_Max|| edgeMap(i+1,j+1) > T_Max)                    
                        edgeMap(i,j) = 255;
                        looping = true;
                    end
                end
            end
        end
    end
    
    for i = 2:h-1,
        for j = 2:w-1,
            if (edgeMap(i,j) < T_Max && edgeMap(i,j) > 0)
                edgeMap(i,j) = 0;
            end
        end
    end
    
    delete(wb);
end