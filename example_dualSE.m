% example_dualSE.m
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de
clear;clc;clear global;
global FOV;
global Mat;
global Kmax2Pi;     %Skope of 3D-K space
global TDImage;     %Output 3D image or 1D signal intensity
TDImage=true;
global plot3proj;
plot3proj=false;

%==========================================================================
%% sequence preparation
%==========================================================================
BasicPulse=4;
gamma = 42.5756; %MHz/T
FOV = 192; %mm
Mat = 64;
Kmax2Pi= ceil(Mat/FOV*2*10)/10; %default Kmax2Pi=3;
Position = [];%phantomP;
RF_BandWidth = 4000; %Hz
Gs = RF_BandWidth/gamma/FOV;%mT/m
Gs_Duration = 4; %ms
KGs = Gs*Gs_Duration*gamma;
TE=76;
TE1hf = 2.11; %ms
TE2hf = TE/2; %ms
TADCcenter=TE2hf-TE1hf;
TR=250;
KSpoil1 = ceil(1*Mat/FOV*1e3);
KSpoil2 = ceil(1*Mat/FOV*1e3);
KSpoil3 = 25*Mat/FOV*1e3;
Gx = [-KSpoil1/TE1hf;  -(KSpoil1+KSpoil2)/TE2hf;   -KSpoil2/TADCcenter; -KSpoil3/(TR-TE)]./gamma;
Gy = [-KSpoil2/TE1hf;  -(KSpoil2+KSpoil1)/TE2hf;   -KSpoil1/TADCcenter; -KSpoil3/(TR-TE)]./gamma;
Gz = [KSpoil1/TE1hf;    (KSpoil1+KSpoil2)/TE2hf;    KSpoil2/TADCcenter; -KSpoil3/(TR-TE)]./gamma;
dG = 40/FOV/gamma.*[1,2,3]';
Gx = Gx+dG(1);  Gy = Gy+dG(2);  Gz = Gz+dG(3);

angle = [33;150;150;0];
axes = [90;0;180;0];
time = [TE1hf;TE2hf;TADCcenter;TR-TE];
echo = [TE1hf;TE2hf;TADCcenter;TR-TE];
%==========================================================================
%% simulation preparation
%==========================================================================
NPulse=BasicPulse*4;
sequence.npulse=NPulse;
sequence.time=col(repmat(time,[1,NPulse/BasicPulse]));% ms
sequence.echo=col(repmat(echo,[1,NPulse/BasicPulse]));% ms
sequence.gradx=col(repmat(Gx,[1,NPulse/BasicPulse]));
sequence.grady=col(repmat(Gy,[1,NPulse/BasicPulse]));
sequence.gradz=col(repmat(Gz,[1,NPulse/BasicPulse]));
sequence.velocityx=zeros(1,NPulse);     % transport velocity in mm/s
sequence.velocityy=zeros(1,NPulse);     % transport velocity in mm/s
sequence.velocityz=zeros(1,NPulse);     % transport velocity in mm/s
sequence.angle=col(repmat(angle,[1,NPulse/BasicPulse]));
sequence.axes=col(repmat(axes,[1,NPulse/BasicPulse]));  % rotation axis.
sequence.diffTensor=[0,0,0,0,0,0,0,0,0];
sequence.T1=1000;                       % ms, 0 equivalent to infinity
sequence.T2=100;                         % ms, 0 equivalent to infinity
sequence.Omega=0;                       % Hz, frequency offset
sequence.ktolerance=1.e-8;              % computational parameter:
sequence.ktolerance_phys=1.e-6;         % 3d k-vectors gridding size 1e-3 = 1 rad/mm
sequence.nOutput=1;
matlabseq = 'matlabseq_1.txt';
matlab2c(sequence, matlabseq);
if ismac
    string = ['./3D_SR-PG_ios.out ',matlabseq];
elseif isunix
    string = ['./3D_SR-PG.out ',matlabseq];
elseif ispc
   string = ['3D_SR-PG_tensor.exe ',matlabseq];
else
    disp('Platform not supported')
end

%==========================================================================
%% simulation execution
%==========================================================================
system(string);
logData=KmapSignal2(['output',int2str(sequence.nOutput),'.txt']);%faster, but need preallocate memory(!).
logData.sequence = sequence;
logData=IMG_DFT(logData,Position);
%==========================================================================
%% OUTPUT
%==========================================================================
%% prep 3D kspace graph
% post-gradding for 3d kmap (due to memory limitation)
Post_grid=[1,1,1]; %post-grid for dKx,dKy,dKz
dK = Post_grid/FOV*1e3;%1/m
[Kgrid3D, Bgrid3D]=prep_3dkmap(logData,dK);

%% plot 3D kspace graph && projected image
outecho=3:BasicPulse:NPulse;
for i=1:length(outecho)
    figure; set(gcf,'outerposition',[0,0,1000,500]);
    subplot(1,2,1);
    plot_3dkmap(squeeze((Kgrid3D(:,:,:,outecho(i)))),Bgrid3D(:,:,:,outecho(i)))
    subplot(1,2,2);
    imagesc(single(abs(logData.sig3d(:,:,outecho(i)))).*imresize(double(80*dicomread('Brain.dcm')),Mat/length(double(dicomread('Brain.dcm')))));
    colormap(gca,'gray')
    suptitle(['TR ',int2str(i)])
end

return
