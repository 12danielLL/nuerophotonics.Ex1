


function [dHbR, dHbO, fig, SNR] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
    %% CalcNIRS - calculate and plot HbR and HbO
    % Input:
    % dataFile - .mat file with intensity data.
    % SDS - Source-Detector Separation distance in cm (default = 3)
    % tissueType - one of the rows in DPFperTissueFile (for example 'adult_forearm' \ 'baby_head' \ 'adult_head' \ 'adult_leg' )
    % plotChannelIdx - vector with numbers in the range of [1-20] indicating channels to plot. If empty - none is plotted. (default = [])
    % extinctionCoefficientsFile - .csv file with the following columns : wavelength, Water, HbO2, HHb, FatSoybean
    % DPFperTissueFile - .txt file with two columns: Tissue and DPF (Tissue is tissue type, corresponding with tissueType input variable)
    % relDPFfile - relative DPF according to wavelength
    %
    % Output :
    % dHbR - HbR concentration change for all channels (nx20) where n is time vector length
    % dHbO - HbO concentration change for all channels (nx20) where n is time vector length
    % fig - handle to figure. Empty if plotChannelIdx == [].
    % SNR - Signal-to-Noise Ratio

    %default values for optional files
    if ~exist('extinctionCoefficientsFile', 'var') || isempty(extinctionCoefficientsFile)
        extinctionCoefficientsFile = '.\ExtinctionCoefficientsData.csv';
    end
    if ~exist('DPFperTissueFile', 'var') || isempty(DPFperTissueFile)
        DPFperTissueFile = '.\DPFperTissue.txt';
    end
    if ~exist('relDPFfile', 'var') || isempty(relDPFfile)
        relDPFfile = '.\RelativeDPFCoefficients.csv';
    end

    %input validation
    if ~ischar(dataFile) || ~exist(dataFile, 'file')
        error('dataFile must be a valid filename. File does not exist: %s', dataFile);
    end
    if ~exist('SDS', 'var') || isempty(SDS)
        SDS = 3;
    elseif ~isnumeric(SDS) || SDS <= 0
        error('SDS must be a positive number');
    end
    if ~ischar(tissueType)
        error('tissueType must be a string');
    end
    if ~exist('plotChannelIdx', 'var') || isempty(plotChannelIdx)
        plotChannelIdx = [];
    elseif ~isvector(plotChannelIdx) || any(plotChannelIdx < 1) || any(plotChannelIdx > 20)
        error('plotChannelIdx must be a vector with values in the range [1-20]');
    end
    if ~ischar(extinctionCoefficientsFile) || ~exist(extinctionCoefficientsFile, 'file')
        error('extinctionCoefficientsFile must be a valid filename. File does not exist: %s', extinctionCoefficientsFile);
    end
    if ~ischar(DPFperTissueFile) || ~exist(DPFperTissueFile, 'file')
        error('DPFperTissueFile must be a valid filename. File does not exist: %s', DPFperTissueFile);
    end
    if ~ischar(relDPFfile) || ~exist(relDPFfile, 'file')
        error('relDPFfile must be a valid filename. File does not exist: %s', relDPFfile);
    end

    %load data
    data = load(dataFile);
    if ~isfield(data, 'SD') || ~isfield(data.SD, 'Lambda') || ~isfield(data, 't') || ~isfield(data, 'd')
        error('dataFile must contain DS.Lambda, t, and d fields');
    end

    wavelengths = data.SD.Lambda;
    time = data.t;
    intensity = data.d;
    if size(intensity, 2) ~= 40
        error('Intensity data must have 40 columns');
    end

   
    intensity1 = intensity(:, 1:20); %columns 1-20 
    intensity2 = intensity(:, 21:40); %columns 21-40 

    extCoeff = readtable(extinctionCoefficientsFile);
    if ~any(extCoeff.wavelength == wavelengths(1)) || ~any(extCoeff.wavelength == wavelengths(2))
        error('Extinction coefficients file does not contain the provided wavelengths');
    end
    eHbO1 = extCoeff{extCoeff.wavelength == wavelengths(1), 'HbO2'} / log(10);
    eHbO2 = extCoeff{extCoeff.wavelength == wavelengths(2), 'HbO2'} / log(10);
    eHbR1 = extCoeff{extCoeff.wavelength == wavelengths(1), 'HHb'} / log(10);
    eHbR2 = extCoeff{extCoeff.wavelength == wavelengths(2), 'HHb'} / log(10);

    dpfData = readtable(DPFperTissueFile, 'ReadVariableNames', false);
    dpfRow = dpfData(strcmp(dpfData.Var1, tissueType), :);
    if isempty(dpfRow)
        error('tissueType not found in DPFperTissueFile');
    end
    dpf807 = dpfRow.Var2;

    relDPF = readtable(relDPFfile);
    if ~any(relDPF.wavelength == wavelengths(1)) || ~any(relDPF.wavelength == wavelengths(2))
        error('Relative DPF file does not contain the provided wavelengths');
    end
    relDPF1 = relDPF{relDPF.wavelength == wavelengths(1), 'relDPFcoeff'};
    relDPF2 = relDPF{relDPF.wavelength == wavelengths(2), 'relDPFcoeff'};

    %calculate DPF 
    dpf1 = dpf807 * relDPF1;
    dpf2 = dpf807 * relDPF2;

    %calculate HbO and HbR changes 
    dHbO = zeros(length(time), 20);
    dHbR = zeros(length(time), 20);
    for ch = 1:20
        A1 = -log(intensity1(:, ch) / mean(intensity1(:, ch)));
        A2 = -log(intensity2(:, ch) / mean(intensity2(:, ch)));
        A = [A1 / (SDS * dpf1), A2 / (SDS * dpf2)];
        E = [eHbO1, eHbO2; eHbR1, eHbR2];
        Hb = E \ A';
        dHbO(:, ch) = Hb(1, :)';
        dHbR(:, ch) = Hb(2, :)';
    end

    %plot
    fig = [];
    if ~isempty(plotChannelIdx)
        fig = figure;
        for idx = 1:length(plotChannelIdx)
            ch = plotChannelIdx(idx);
            subplot(length(plotChannelIdx), 1, idx);
            plot(time, dHbR(:, ch) * 1e6, 'b', time, dHbO(:, ch) * 1e6, 'r');
            title(sprintf('Channel %d', ch));
            xlabel('Time (s)');
            ylabel('Concentration Change (\muM)');
            legend({'HbR', 'HbO'});
        end
    end

    %fourier transform 
    signal = dHbR(:, 1);
    L = min(length(time), 1000); 
    Fs = 1 / (time(2) - time(1));
    f = Fs * (0:(L/2)) / L;  
    Y = fft(signal, L);
    P2 = abs(Y / L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    %plot the FFT
    if ~isempty(plotChannelIdx)
        figure;
        plot(f(2:end), P1(2:end)); 
        title('Single-Sided Amplitude Spectrum of Channel 1');
        xlabel('Frequency (Hz)');
        ylabel('|P1(f)|');
        xlim([0.65, 2]);
    end

   
    %calculate SNR for Channel 1 
    [~, idx1_5Hz] = min(abs(f - 1.5));
    signalPower1_5Hz = P1(idx1_5Hz); 
    noisePowerAbove2_5Hz = mean(P1(f > 2.5));
    SNR1 = signalPower1_5Hz / noisePowerAbove2_5Hz;

    fprintf('SNR for Channel 1: %.2f\n', SNR1);
end
