function Power = calculateWaveletPower(cmplxWaveletData)
%CALCULATEWAVELETPOWER
channels = size(cmplxWaveletData,1);
freqCount = size(cmplxWaveletData,2);

for ch=1:channels
    for fi=1:freqCount
        % put power data into big matrix
        as = squeeze(cmplxWaveletData(ch,fi,:));
        Power(ch,fi,:) = mean(abs(as).^2,2);
        clear as;
    end
end



