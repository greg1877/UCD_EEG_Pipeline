clearvars; close all;

%% NOTES:
% Author: Greg Bales
% email: glbales@ucdavis.edu
%
% 
%% START EEG DATA PIPELINE
workingDirectory = pwd;

%
%%%% Choose from the following data files %%%%
% 2024-03-04-AFP28_sample_EEG_LSL_Stream.mat
% 2024-03-19-AFP31_sample_EEG_LSL_Stream.mat
% 2024-03-20-AFP30_sample_EEG_LSL_Stream.mat
% 2024-04-03-AFP32_sample_EEG_LSL_Stream.mat
lsl_sampleFileName = '2024-04-03-AFP32_sample_EEG_LSL_Stream.mat';

load([workingDirectory,'\Sample_Data\',lsl_sampleFileName]);

%% %%%%%%% CREATE MONTAGE FROM EASYLOC AND GTEC FILES %%%%%%%%
% Here we collect the channel names, numbers and locations.

addpath(genpath([workingDirectory,'\CSDtoolbox']));
addpath([workingDirectory,'\Wavelet_Functions']);
addpath([workingDirectory,'\Support_Files']);

clear eeg timeHold easyLocData eegFileName Mon;

load CU_montage.mat;
eeg.montage_info = Mon;

% get the eeg data from the EasyLoc file.
eegFileName = 'CORRECT_EASYLOC_FILE_061824.txt';
easyLocData = readtable(eegFileName);

allchannels = easyLocData.Var1;

% The exclusion rule: All channels that begin with "Ch" should be removed
nullChannels=[];
for channelCount=1:length(easyLocData.Var1)
    clear startHold endHold;
    [startHold,endHold]=regexp(easyLocData.Var4{channelCount},'Ch');
    if startHold==1 & endHold==2
        nullChannels = [nullChannels,channelCount];
    end
end
activeChannels = setdiff(allchannels,nullChannels);

% Get the Channel Names and Numbers from the EasyLoc File.
channelNumbers = easyLocData.Var1(activeChannels);
channelNames = easyLocData.Var4(activeChannels); 

% Create the montage data from the g_tec montage file.
for i=1:length(channelNames)
    for j=1:length(eeg.montage_info.electrodename)
        if strcmpi(channelNames{i},eeg.montage_info.electrodename{j})
            eeg.channel_labels{i} = eeg.montage_info.electrodename{j};
            eeg.electrodeX(i) = eeg.montage_info.xposition(j);
            eeg.electrodeY(i) = eeg.montage_info.yposition(j);
            eeg.electrodeZ(i) = eeg.montage_info.zposition(j);
            eeg.electrodeX2d(i) = eeg.montage_info.xposition_2d(j);
            eeg.electrodeY2d(i) = eeg.montage_info.yposition_2d(j);
        end
    end
end

% HEADMAP SHOWING CHANNEL POSITIONS
headPlot(eeg.electrodeX2d,eeg.electrodeY2d,eeg.electrodeZ+1,eeg.channel_labels);

%% %%%%%%%%%%%%%% GET RAW EEG FROM STREAM DATA %%%%%%%%%%%%%%%

eeg.time = EEG.time_stamps;
eeg.data.raw = EEG.time_series(activeChannels,:)';
eeg.data.raw = cast(eeg.data.raw,'double'); %<- Data must be cast to double in order to use filtfilt


%% %%%%%%%%%%%%% ZERO PHASE FILTER THE DATA %%%%%%%%%%%%%%%%%%

bandpassFlag = true;    % If not bandpass, then highpass filter
butterWorthFlag = true;  % If not butterworth, then iir filter
filterState = 2*bandpassFlag + butterWorthFlag; %<-convert to decimal value

% You can select any frequencies you wish. I have not determined which are
% best at this point.
effSampleRate = round(1./mean(diff(eeg.time)));
c_freq_low = .5; 
c_freq_high = 55; 
pb_ripple = .1;
butter_f_order = 4;
iir_f_order = 8;
filterParams = [effSampleRate,c_freq_low,c_freq_high,pb_ripple,butter_f_order,iir_f_order];

% filter the raw voltage data
eeg.data.flt = filter_eeg_data(eeg.data.raw,filterParams,filterState);

%% %%%%%%%%%%%%% check the filtered output %%%%%%%%%%%%%%%%%%%
compareResample = true;
if compareResample
    clear clrMap pUnfilt fUnfilt pFilt fFilt;

    [pUnfilt,fUnfilt] = pwelch(eeg.data.raw,[],[],[],effSampleRate);
    [pFilt,fFilt] = pwelch(eeg.data.flt,[],[],[],effSampleRate);

    clrMap = colormap(jet(20));
    figure; tiledlayout('flow');
    nexttile; hold on; title('Raw EEG','Interpreter','latex');
    plot(eeg.time-eeg.time(1), eeg.data.raw,'LineWidth',2);
    xlabel('sec','Interpreter','latex'); ylabel('$\mu V$','Interpreter','latex');
    set(gca,'Box','on','TickLabelInterpreter','latex','ColorOrder',clrMap);
    l=legend(channelNames); l.Interpreter='latex'; l.FontSize=14;

    nexttile; hold on; title('EEG Detail','Interpreter','latex');
    plot(EEG.time_stamps - EEG.time_stamps(1), eeg.data.flt);
    xlabel('sec','Interpreter','latex');
    ylabel('$\mu V$','Interpreter','latex');
    set(gca,'Box','on','TickLabelInterpreter','latex','ColorOrder',clrMap);
    l=legend(channelNames); l.Interpreter='latex'; l.FontSize=14;

    nexttile; hold on; title('Spectrum of Raw EEG','Interpreter','latex');
    plot(fUnfilt,pUnfilt);
    set(gca,'XScale','log','YScale','log');
    xlabel('freq (Hz)','Interpreter','latex');
    ylabel('$\mu V^2/Hz$','Interpreter','latex');
    set(gca,'Box','on','TickLabelInterpreter','latex','ColorOrder',clrMap);
    xline(c_freq_low,'k','LineWidth',2);
    xline(60,'r','LineWidth',2);
    if bandpassFlag
        xline(c_freq_high,'k','LineWidth',2);
    end
    l=legend(channelNames); l.Interpreter='latex'; l.FontSize=14;

    nexttile; hold on; title('Spectrum of EEG After Filtering','Interpreter','latex');
    plot(fFilt,pFilt);
    set(gca,'XScale','log','YScale','log');
    xlabel('freq (Hz)','Interpreter','latex');
    ylabel('$\mu V^2/Hz$','Interpreter','latex');
    set(gca,'Box','on','TickLabelInterpreter','latex','ColorOrder',clrMap);
    xline(c_freq_low,'k','LineWidth',2);
    xline(60,'r','LineWidth',2);
    if bandpassFlag
        xline(c_freq_high,'k','LineWidth',2);
    end
    l=legend(channelNames); l.Interpreter='latex'; l.FontSize=14;
end

%% %%%%%%%%%%%%%%%%% SURFACE LAPLACIAN %%%%%%%%%%%%%%%%%%%%%%%
%%%%% CSD Spherical Spline %%%%%
[theta,phi,R] = cart2sph(eeg.electrodeX,eeg.electrodeY,eeg.electrodeZ);
theta = rad2deg(theta); phi=rad2deg(phi);  %<--correction for use of the CSD toolbox;
lab = eeg.channel_labels;
xy = [eeg.electrodeX(:),eeg.electrodeY(:)];
MA.phi = phi'; MA.theta = theta'; MA.xy = xy./2+.5; MA.lab = lab';
[G,H,D,pos_B] = GetGH(MA);

% save the parameters in the eeg structure.
lambda = 1.0e-5; eeg.CSD_info.lambda = lambda;
eeg.CSD_info.G = G; 
eeg.CSD_info.H = H;

%%%%% Take the surface lap of the raw EEG signal %%%%%
[eeg.data.fltLap,Y_Voltage_reconstruct] = CSD(eeg.data.flt', G, H, lambda);
eeg.data.fltLap = eeg.data.fltLap'; 
Y_Voltage_reconstruct = Y_Voltage_reconstruct';

%% %%%%%%%%%%%%% check the laplacian reconstruction %%%%%%%%%% 
lapDataCheck = true;
if lapDataCheck
    comparison_channel = 1;
    figure; hold on; title('Check Reconstruction of EEG Voltage From Calculated Current Source Density','Interpreter','latex');
    plot(eeg.time,eeg.data.flt(:,comparison_channel),'r');
    plot(eeg.time,Y_Voltage_reconstruct(:,comparison_channel),'b.');
    xlabel('sec','Interpreter','latex'); ylabel('$\mu V$','Interpreter','latex');
    set(gca,'Box','on','TickLabelInterpreter','latex');
    l=legend('Data','Reconstruction'); l.Interpreter='latex'; l.FontSize=14;
end

%% %%%%%%%%%%%%%%%%% WAVELET TRANSFORM %%%%%%%%%%%%%%%%%%%%%%%

clear time_freq_output wvltTransform;

% I would keep the FWHP dt the same.  However, you are free to select the
% frequency bandwidths and number of analysis frequencies as you see fit.
% These are fill in values for now.
fwhp_dt = .2;  % <-Full width half power time in seconds
min_freq = 4;   % <-Minimum frequency in Hz
max_freq = 40;  % <-Maximum frequency in Hz
num_frex = 10;  % <-Number of equally spaced frequecies in the transform
time_freq_output = morletTimeFreqTransform(eeg.data.fltLap', min_freq, max_freq, num_frex,...
    fwhp_dt, effSampleRate);

wvltTransform.info.fwhp_dt = fwhp_dt;
wvltTransform.info.min_freq = min_freq;
wvltTransform.info.max_freq = max_freq;
wvltTransform.info.num_frex = num_frex;

wvltTransform.data = time_freq_output;

%% %%%%%%%%%%%%% check the results of the wavelet transform %%
waveletDataCheck = true;
if waveletDataCheck

    % Select channels and frequecies to explore
    comparison_channel = [1 3 10]; %<- Enter any number of channels you wish to plot
    comparison_freq = [1 2 3 5 8 9 10]; %<- Enter any number of frequencies you wish to plot

    for channelCount = 1:length(comparison_channel)
        legendText{channelCount} = eeg.channel_labels{channelCount};
    end
    
    % Plot the wavelet filtered timeseries
    figure; tiledlayout('flow');
    for freqCount = 1:length(comparison_freq)
        nexttile; hold on; title(['Wavelet Filtered Data. Center Freq = ',...
            num2str(wvltTransform.data.freq(comparison_freq(freqCount))),' Hz'],'Interpreter','latex');
        plot(eeg.time,squeeze(wvltTransform.data.real(comparison_channel,comparison_freq(freqCount),:)));
        xlabel('sec','Interpreter','latex'); ylabel('$\mu V$','Interpreter','latex');
        set(gca,'Box','on','TickLabelInterpreter','latex');
        l=legend(legendText); l.Interpreter='latex'; l.FontSize=14;
    end

    % Plot the wavelet filtered power for all channels at selected frequencies
    figure; tiledlayout('flow');
    for freqCount = 1:length(comparison_freq)
        nexttile; hold on; title(['Wavelet Filtered Power. Center Freq = ',...
            num2str(wvltTransform.data.freq(comparison_freq(freqCount))),' Hz'],'Interpreter','latex');
        imagesc('XData',eeg.time,'YData',1:1:length(eeg.electrodeY),'CData',...
            log10(squeeze(wvltTransform.data.power(:,comparison_freq(freqCount),:))),'Interpolation','bilinear');
        xlabel('sec','Interpreter','latex'); ylabel('Channel','Interpreter','latex');
        set(gca,'Box','on','TickLabelInterpreter','latex',...
            'YLim',[0 length(eeg.electrodeY)],'XLim',[eeg.time(1) eeg.time(end)],...
            'YTick',1:1:length(eeg.electrodeY),'YTickLabel',eeg.channel_labels);
    end

    % Plot the wavelet filtered power for all frequencies at selected channels
    figure; tiledlayout('flow');
    for channelCount = 1:length(comparison_channel)
        nexttile; hold on; title(['Wavelet Filtered Power. Channel ',eeg.channel_labels(comparison_channel(channelCount))],'Interpreter','latex');
        imagesc('XData',eeg.time,'YData',wvltTransform.data.freq,'CData',...
            log10(squeeze(wvltTransform.data.power(comparison_channel(channelCount),:,:))),'Interpolation','bilinear');
        xlabel('sec','Interpreter','latex'); ylabel('Frequency (Hz)','Interpreter','latex');
        set(gca,'Box','on','TickLabelInterpreter','latex',...
            'YLim',[wvltTransform.data.freq(1) wvltTransform.data.freq(end)],'XLim',[eeg.time(1) eeg.time(end)],...
            'YTick',wvltTransform.data.freq,'YTickLabel',wvltTransform.data.freq);
    end

end 


%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputFigHandle = headPlot(channelX,channelY,channelZ,channelLabels)

[imA,mapA]=imread('brain_regions_colored_eeg_light.png');
imA = flipud(imA); mapA = flipud(mapA);

%%%%%%%%%%%% FROM HEADPLOT 
headrad = .5;
% Plot ears and nose
base  = headrad-.0046;
basex = 0.18*headrad;                   % nose width
tip   = 1.15*headrad;
tiphw = .04*headrad;                    % nose tip half width
tipr  = .01*headrad;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % headrad = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = 2.1;
circ = 0:.01:2*pi+.1;
rx = 1.05.*sin(circ);
ry = 1.05.*cos(circ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputFigHandle = figure; 
set(gcf,'Color','w');
tiledlayout('flow','Padding','compact','TileSpacing','none');

%%%%%%% This is a headplot

nexttile; hold on; %title(titleText, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
image(imA,'XData',[-1.05 1.05],'YData',[-1.05 1.05]);
plot3(channelX,channelY,channelZ,'o','MarkerSize',6,'MarkerFaceColor','k');
axis('off'); 
axis([-1.3 1.3 -1.3 1.3]);
set(gca,'Box','on');  axis square;
plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,2*ones(size([basex;tiphw;0;-tiphw;-basex])),'Color','k','LineWidth',2,'hittest','off');                 % plot nose
plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color','k','LineWidth',2,'hittest','off')    % plot left ear
plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color','k','LineWidth',2,'hittest','off')   % plot right ear
plot(rx,ry,'Color','k','LineWidth',2);
for channelCount=1:length(channelLabels)
    text(channelX(channelCount),channelY(channelCount),channelLabels(channelCount),...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontSize',10,'Interpreter','latex');
end


end % END FUNCTION headplot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputFilteredData = filter_eeg_data(rawData,filterParams,filterState)

% filterParams = [effSampleRate,c_freq_low,c_freq_high,pb_ripple,butter_f_order,iir_f_order];
effSampleRate = filterParams(1);
c_freq_low = filterParams(2); 
c_freq_high = filterParams(3); 
pb_ripple = filterParams(4);
butter_f_order = filterParams(5);
iir_f_order = filterParams(6);

switch filterState
    case 0 % highpass, iir filter
        flt = designfilt('highpassiir','FilterOrder',iir_f_order, ...
             'PassbandFrequency',c_freq,'PassbandRipple',pb_ripple, ...
             'SampleRate',effSampleRate);
        outputFilteredData= filtfilt(flt,rawData);
    case 1 % highpass, butterworth filter
        [bP,aP] = butter(butter_f_order, c_freq_low/(effSampleRate/2),'high');
        outputFilteredData = filtfilt(bP,aP,rawData);
    case 2 % bandpass, iir filter
        flt = designfilt('bandpassiir','FilterOrder',iir_f_order, ...
            'HalfPowerFrequency1',c_freq_low,'HalfPowerFrequency2',c_freq_high, ...
            'SampleRate',effSampleRate);
        outputFilteredData = filtfilt(flt,rawData);
    case 3 % bandpass, butterworth filter
        [bP,aP] = butter(butter_f_order, [c_freq_low c_freq_high]/(effSampleRate/2),'bandpass');
        outputFilteredData = filtfilt(bP,aP,rawData);
end
end % END FUNCTION filter_eeg_data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputSampleRateData = getSampleRateData(timeVector)

timeDiff = diff(timeVector); sampleRate = 1./timeDiff;
sampleRate_mean = mean(sampleRate); sampleRate_mode = mode(sampleRate);
sampleRate_var = var(sampleRate);
sampleRate_shiftPercent = 100*(sampleRate_mean-sampleRate_mode)/sampleRate_mode;

outputSampleRateData = [sampleRate_mean, sampleRate_mode, sampleRate_shiftPercent, sampleRate_var]';

end % END FUNCTION getSampleRateData



