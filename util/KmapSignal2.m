function logData = KmapSignal2( dataFile )
fid = fopen(dataFile,'r');
bufferSize = 1e20;
% eol = newline;
% tic;  A = textscan(fid,'%f %f %f %f %f %f');  toc;

% dataBatch = fread(fid,bufferSize,'uint8=>char')';
% dataIncrement = fread(fid,1,'uint8=>char');
% while ~isempty(dataIncrement) && (dataIncrement(end) ~= eol) && ~feof(fid)
%     dataIncrement(end+1) = fread(fid,1,'uint8=>char');  %This can be slightly optimized
% end
% data = [dataBatch dataIncrement];
% 
% while ~isempty(data)
%     A = sscanf(data,'%f %f %f %f %f %f');
%     scannedData = reshape(sscanf(data,'%f %f %f %f %f %f'),6,[])';
%     
%     dataBatch = fread(fid,bufferSize,'uint8=>char')';
%     dataIncrement = fread(fid,1,'uint8=>char');
%     while ~isempty(dataIncrement) && (dataIncrement(end) ~= eol) && ~feof(fid)
%         dataIncrement(end+1) = fread(fid,1,'uint8=>char');%This can be slightly optimized
%     end
%     data = [dataBatch dataIncrement];
% end

data = fread(fid,bufferSize,'uint8=>char')';
scannedData = textscan(data,'%f32 %f32 %f32 %f32 %f32 %f32'); %single instead of double
fclose(fid);

%% output nspins per echo
NKperEcho=[];
fid2 = fopen([dataFile(1:end-4),'_Nspins.txt'],'r');
eof               = false;
while ~eof
    line = fgetl(fid2);
    if line == -1
        % End of file.
        disp('Reached end of file');
        eof = true;
    elseif isempty(line)
        nEcho = nEcho+1;
        disp(['Echo = ',int2str(nEcho)]);
    else
        nSpins=textscan(line,'%d');
        NKperEcho = [NKperEcho;nSpins{1}];
    end
end

logData.kvector=[scannedData{1} scannedData{2} scannedData{3}];
logData.signal=[scannedData{4} scannedData{5} scannedData{6}];
logData.NKperEcho=NKperEcho;
logData.KmapSignal=2; %use this KmapSignal2 format
end