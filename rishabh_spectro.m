function STFTMdb = rishabh_spectro(y,Fs,WinD,shiftD,nfft)
%nfft = number of nfft points
%filename is the audio file input
%WinD is the window length in seconds
%shiftD is the hop size in seconds, overlap = WinD - shiftD
%[y,Fs] = audioread(filename);
WinL = floor(WinD*Fs); % in samples
shiftL = floor(shiftD*Fs); % hop size in samples
nFr  = round(length(y)/shiftL); %no., of frames
win = hamming(WinL);
STFT = [];
%STFTM = zeros((nfft/2)+1,nFr);
for c = 1:nFr - round(WinL/shiftL)
    FB = (c-1)*shiftL+1; % beginning of the frame in samples
    FE = FB + WinL -1;   % ending of the frame in samples
    wseg = y(FB:FE).*win; % selcting the appropriate window with FB and FE as window boundaries
    STFT(:,c) = fft(wseg,nfft); % taking nfft points FFT within wseg
end
STFTM = abs(STFT(1:((nfft/2)+1),:)); % as FFT is symmetrical ,considering first n points
STFTMdb = 20*log10(STFTM); % conversion to dB scale
faxis = (0:nfft/2)*Fs/nfft; 
naxis = (0:size(STFTM,2)-1)*shiftD; % in seconds
STFTMdbmax = max(STFTMdb(:));
dbdown = 60; %deciding the range of the plot
caxisl = [STFTMdbmax-dbdown STFTMdbmax];% limiting the range of STFT values
imagesc(naxis,faxis,STFTMdb,caxisl);
axis xy;
%ylabel('Frequency');
%xlabel('Duration of signal in Seconds');
colorbar;
end
