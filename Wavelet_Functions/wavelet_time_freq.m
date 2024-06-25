function [tf,complexOutput,wavelet,waveletFFT,waveletTime] = wavelet_time_freq(freqX,X,num_cycles,sampleRate)

% other wavelet parameters
% freqX = linspace(min_freq,max_freq,num_frex);
waveletTime = -5:(1/sampleRate):5;
half_wave = (length(waveletTime)-1)/2;

% FFT parameters
nKern = length(waveletTime);
nData = length(X);
nConv = nKern+nData-1;

dataX = fft(X,nConv);
%% loop over cycles

for fi=1:length(freqX)
    % create wavelet and get its FFT
    s = num_cycles(fi)/(2*pi*freqX(fi));
    
    cmw  = exp(2*1i*pi*freqX(fi).*waveletTime) .* exp(-waveletTime.^2./(2*s^2));
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % run convolution, trim edges, and reshape to 2D (time X trials)
    as1 = ifft(cmwX.*dataX,nConv);
%     as = 2.*as1(half_wave+1:end-half_wave);
    as = as1(half_wave+1:end-half_wave); %DO NOT include the 2 in this function.
    %         as = as(1:length(X));
    as = reshape(as,nData,1);
    
    % put power data into big matrix
    tf(fi,:) = mean(abs(as).^2,2);
    complexOutput(fi,:) = as;
    waveletFFT(fi,:) = cmwX;
    wavelet(fi,:) = cmw;
    clear as as1 s cmw cmwX;
end


end