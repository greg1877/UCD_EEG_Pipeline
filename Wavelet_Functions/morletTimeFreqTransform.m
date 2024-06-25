function time_freq_output = morletTimeFreqTransform(data, min_freq, max_freq, num_frex,...
    fwhp_dt, sampleRate)

% Wavelet Parameters
dT = 1/sampleRate;
channelCount = size(data,1);

freqX = linspace(min_freq,max_freq,num_frex);
num_cycles = sqrt(2/log(2))*pi*fwhp_dt.*freqX;

% Get the Wavelet Data
for i=1:channelCount
    data4wavelet = data(i,:);
    if i==1
    [~,complexOutput,wavelet,waveletFFT,waveletTime] = ...
        wavelet_time_freq(freqX,data4wavelet,num_cycles,sampleRate);
    else 
        [~,complexOutput,~,~,~] = ...
        wavelet_time_freq(freqX,data4wavelet,num_cycles,sampleRate);
    end
%     filterDataCmplx(i,:,:) = complexOutput;
    filterData(i,:,:) = 2.*real(complexOutput); % The two is for correction.  See wavelet_filter_study
    filterDataPower(i,:,:) = abs(complexOutput).^2;
    dataPhase = angle(complexOutput);
    L = size(dataPhase,1);
    for j=1:L
        dataHold = dataPhase(j,:);
        dataHold(find(dataHold<0))=dataHold(find(dataHold<0))+2*pi;
        dataHold=unwrap(dataHold);
        filterDataPhase(i,j,:)=dataHold;
        clear dataHold;
    end
%     angleDataHold = squeeze(filterDataPhase(i,:,:));
%     omega = diff(angleDataHold,1,2)./dT; 
%     filterDataOmega(i,:,:) = omega;
    clear data4wavelet complexOutput omega angleDataHold dataPhase;
end

% output data structure
time_freq_output.real = filterData;
time_freq_output.power = filterDataPower;
time_freq_output.phase = filterDataPhase;
% time_freq_output.omega = filterDataOmega;
% time_freq_output.cmplx = filterDataCmplx;
time_freq_output.freq = freqX;
time_freq_output.wavelet = wavelet;
% time_freq_output.waveletFFT = waveletFFT;
time_freq_output.waveletTime = waveletTime;
time_freq_output.wavelet_hz = linspace( 0, sampleRate/2, floor(length(waveletTime)/2)+1 );
end

