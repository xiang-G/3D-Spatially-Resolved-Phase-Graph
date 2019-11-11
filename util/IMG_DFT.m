%% DFT
function slogData=IMG_DFT(varargin)
slogData = varargin{1};
global FOV;
global Mat;
global TDImage;
res = FOV/Mat;
%% simulate signal on which position
if nargin <2 || isempty(varargin{2})
        FOVx=-FOV/2+res:res:FOV/2;
        FOVy=-FOV/2+res:res:FOV/2;
        FOVz=-FOV/2+res:res:FOV/2;
else
    FOVp= varargin{2};
    if size(FOVp,1)==1
        FOVp=FOVp';
        FOVx=FOVp(1,:);FOVy=FOVp(2,:);FOVz=FOVp(3,:);
    else
        FOVx=FOVp(1,:);FOVy=FOVp(2,:);FOVz=FOVp(3,:);
    end
end
N=length(FOVx);
NE=length(slogData.NKperEcho);
%% signal indensity (x dim without deltaG, y dim with deltaG)
if ~TDImage
    disp(['K_filter = ',num2str(2*pi/(res*1e3)/2)]);
    pos = zeros(3,N);
    posx=pos;posx(1,:)=FOVx;   
    sigx = zeros(N,NE-1); 
    for i=1:NE-1
        [sigx(:,i)]=DFT(slogData,posx*1e3,res*1e3,slogData.NKperEcho(i:i+1)); %signal without motion
    end
    %write back
    slogData.sigxy = (sigx);
    slogData.sig3d = abs(sigx);
    slogData.pha3d = angle(sigx);
else
    %% 2D Image (x-y)
    pos = zeros(3,N*N);
    posxy=pos;
    posxy(1,:)=repmat(FOVx,[1,N]);
    posxy(2,:)=col(repmat(FOVy',[1,N])');
    posxy(3,:)=0;%ceil(FOV/2);
    disp(['K_filter = ',num2str(2*pi/(res*1e3)/2)]); % don't foget to divid by 2
    sigxy = zeros(N*N,NE-1);
    for i=1:NE-1
        [sigxy(:,i)]=DFT(slogData,posxy*1e3,res*1e3,slogData.NKperEcho(i:i+1)); %signal with motion
    end
    %% simulated signal
    ImgPerExy = abs(reshape(sigxy,[N,N,NE-1]));%
    PhasePerExy = angle(reshape(sigxy,[N,N,NE-1]));%
    
    %% 2D Image (y-z)
    pos = zeros(3,N*N);
    posxy=pos;
    posxy(1,:)=0;%ceil(FOV/2);
    posxy(2,:)=repmat(FOVy,[1,N]);
    posxy(3,:)=col(repmat(FOVz',[1,N])');
    disp(['K_filter = ',num2str(2*pi/(res*1e3)/2)]);
    sigxy = zeros(N*N,NE-1);
    for i=1:NE-1
        [sigxy(:,i)]=DFT(slogData,posxy*1e3,res*1e3,slogData.NKperEcho(i:i+1)); %signal with motion
    end
    %% simulated signal
    ImgPerEyz = abs(reshape(sigxy,[N,N,NE-1]));%
    PhasePerEyz = angle(reshape(sigxy,[N,N,NE-1]));%
    
    %% 2D Image (x-z)
    pos = zeros(3,N*N);
    posxy=pos;
    posxy(1,:)=repmat(FOVx,[1,N]);
    posxy(2,:)=0;%ceil(FOV/2);
    posxy(3,:)=col(repmat(FOVz',[1,N])');
    disp(['K_filter = ',num2str(2*pi/(res*1e3)/2)]);
    sigxy = zeros(N*N,NE-1);
    for i=1:NE-1
        [sigxy(:,i)]=DFT(slogData,posxy*1e3,res*1e3,slogData.NKperEcho(i:i+1)); %signal with motion
    end
    %% simulated signal
    ImgPerExz = abs(reshape(sigxy,[N,N,NE-1]));%
    PhasePerExz = angle(reshape(sigxy,[N,N,NE-1]));%
    
    %%
    slogData.sig3d = cat(3,ImgPerExy,ImgPerEyz,ImgPerExz);
    slogData.pha3d = cat(3,PhasePerExy,PhasePerEyz,PhasePerExz);
end
end

function [s]=DFT(slogData,pos,resseq,NKlength) %for KmapSignal
if slogData.KmapSignal == 1
    kvec = slogData.kvector(1+NKlength(1):NKlength(2)-1,:);
    sk = slogData.signal(1+NKlength(1):NKlength(2)-1,1)+1i*slogData.signal(1+NKlength(1):NKlength(2)-1,2);
elseif slogData.KmapSignal == 2
    kvec = slogData.kvector(1+NKlength(1):NKlength(2),:);
    sk = slogData.signal(1+NKlength(1):NKlength(2),1)+1i*slogData.signal(1+NKlength(1):NKlength(2),2);
end
kvecF = kvec((abs(kvec(:,1))<=2*pi/resseq/2)&(abs(kvec(:,2))<=2*pi/resseq/2)&(abs(kvec(:,3))<=2*pi/resseq/2),:);
skF = sk((abs(kvec(:,1))<=2*pi/resseq/2)&(abs(kvec(:,2))<=2*pi/resseq/2)&(abs(kvec(:,3))<=2*pi/resseq/2),:);

s = skF'*exp(-1i*kvecF*pos);
end
