function [Lv, CIE_x, CIE_y, GammaValues, resampledSpectra, desiredwl] = CalibrateColour_Viewpixx_no_goggles

% Load & process the Viewpixx calibration data.
% The processed data are saved in a format that can be loaded for use in LMS2RGB_Vpixx.m,
% the program that converts LMS coordinates into RGB values appropriate for the Viewpixx (or PROpixx).
% It also produces plots of the luminance values + gamma fits and plots the CIE x,y values.
% This version for data without the goggles, Viewpixx in EEG room B114, obtained 9/9/16
% OR the new Viewpixx (Viewpixx2) in Room B116, also no goggles, 9/9/16.
% R Maloney, September 2016
% Modified 7/8/17 for fresh spectral measurements made on the ViewpixxEEG room B114 for Miaomiao

DoGammaFits = 1; %Flag for whether or not you want to fit gamma to the luminance data.

%Compute the Projected source area; ie the size of the screen, ***in m^2***
ProjA = 0.293 * 0.52; % the dimensions, in m, of the Viewpixx screen.
nLevels = 16; %the number of luminance increments per gun: important to know beforehand

%We'll assume we are inside the 'Calibration' directory
%DataDir = fullfile('Ocean Optics spectrometer data','Viewpixx2_09_09_16_no_goggles');   % B116 Viewpixx data location
%DataDir = fullfile('Ocean Optics spectrometer data','ViewpixxEEG_09_09_16_nogoggles'); % B114 Viewpixx data location
DataDir = fullfile('Ocean Optics spectrometer data', 'ViewpixxEEG_07_08_17_no_goggles'); % B114 Viewpixx data location

% File name to save the processed data in:
Proc_fName = 'ViewpixxEEG_Processed_cal_data_7_8_2017.mat';
PlotTitle = 'ViewpixxEEG calibration, no goggles, 7/8/17';

PhotometryName = 'Photometry-2017_8_7-';   % first part of photometry file name.
ColorimetryName = 'ColorSample-2017_8_7-'; % first part of colorimetry file name.

% Each measurement was made in the order R, G, B, K * nLevels, and the suffix on the end of each measurement
% indicates the order in which it was acquired, & hence what test stimulus it belongs to.
% So, since we used 16 nLevels, there should be 4*16 = 64 such files, with suffixes from 1:64.
% These are tallied using the 'TotalMeas' counter.
PlotCols = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5]; % defines the order of the colours, for plotting later.

TotalMeas = 0; %this is a counter for the total number of measurements made across channels

for ThisGunSet = 1:4 %loop across R, G, B & K datasets
    
    for ThisLumLevel = 1:nLevels % Loop across different luminance levels.
        
        TotalMeas = TotalMeas+1; %increment total measurement counter
        
        % ------------------------------ %
        % First, load the photometry data:
        % ------------------------------ %
        
        fName = fullfile(DataDir, [PhotometryName num2str(TotalMeas) '.txt']);
        % Open the text file.
        fileID = fopen(fName,'r');
        % Read columns of data according to format string.
        dataArray = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', '\t',  'ReturnOnError', false);
        % Close the text file.
        fclose(fileID);
        
        % Set aside the luminance.
        % We do that by taking the value in Candela and dividing that by the projected source area (ie screen size).
        % This gives us Luminance in Cd/m^2 (see https://en.wikipedia.org/wiki/Luminance)
        % Candela is the 7th value in the numeric data column (need to convert to double from cell/char first)
        Lv(ThisGunSet,ThisLumLevel) = str2double(cell2mat(dataArray{2}(7))) / ProjA;
        
        % ------------------------------ %
        % Now for the colorimetry data:
        % ------------------------------ %
        fName = fullfile(DataDir, [ColorimetryName num2str(TotalMeas) '.txt']);
        startRow = 5;
        % Open the text file.
        fileID = fopen(fName,'r');
        % Read columns of data according to format string.
        textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
        dataArray = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', '\t', 'ReturnOnError', false);
        % Close the text file.
        fclose(fileID);
        
        % Take the CIE x value (the 4th row in the array)
        CIE_x(ThisGunSet,ThisLumLevel) = str2double(cell2mat(dataArray{2}(4)));
        % And the CIE y value (the 5th row)
        CIE_y(ThisGunSet,ThisLumLevel) = str2double(cell2mat(dataArray{2}(5)));
        % We don't need CIE z because we can compute it thus: 1-x-y
    end
end

% Plot the luminance data
figure(1)
for ThisGunSet = 1:4
    plot(linspace(0,1,nLevels), Lv(ThisGunSet,:) , '.', 'Color', PlotCols(ThisGunSet,:))
    hold on
end
ylabel('Luminance (Cd/m^2)')
xlabel('RGB output level')
title (PlotTitle)

% ---------------------------------- %
% Fit the gamma fits here, if requested.
% ---------------------------------- %
if DoGammaFits
    estimate = [2 1]; % provide estimates for the exponent (gamma) and constant parameters, respectively, in the power function
    % normalise the measurements first before fitting the gamma!
    MaxLums = max(Lv,[],2); % Obtain the maximum luminances for each channel (hopefully this is also the last measured value)
    % Normalise element-by-element, across each row in turn
    Norm_measurements = bsxfun(@rdivide, Lv, MaxLums);
    
    % Fit the power function to each measured channel using 'nlinfit':
    for ii = 1:4
        fittedParameters = nlinfit(linspace(0,1,nLevels), Norm_measurements(ii,:), @GammaPower, estimate);
        GammaValues(ii, :) = fittedParameters;
    end
    
    % Now plot the fitted functions back on top of the raw measurements.
    % They have to be re-scaled by the max. luminances of the raw values because the fits were performed on the normalized values
    contX = [0:1/255:1]'; %continuous values for x
    for ii=1:4
        figure(1)
        plotY = GammaPower(GammaValues(ii,:), contX) * MaxLums(ii);
        hold on
        plot(contX, plotY, '-', 'Color', PlotCols(ii,:));
    end
end

% Plot the CIE x, y values just to check they are roughly in the right place in
% CIE color space
figure(2)
for ThisGunSet = 1:4
    plot(CIE_x(ThisGunSet,:), CIE_y(ThisGunSet,:) , '.', 'Color', PlotCols(ThisGunSet,:))
    hold on
end
axis([0 0.8 0 0.9]); % limit the axes to those of the CIE chart (but note some anomalous values might be outside this range).
title('CIE 1931 chromaticity values')
xlabel('x')
ylabel('y')

% --------------------------------
% Load and plot the spectra.
% --------------------------------

% We will load the spectra obtained for the R, G, B & K channels.
% Remember that they were presented in 1-16 steps in that order, but incremented each time,
% so 16 * 4 = 1:64 spectra.
% we just want those for the max luminance.

% Viewpixx 2: B116
% SpectName{1} = 'Jaz2Spect_09_09_16_R_16.txt';
% SpectName{2} = 'Jaz2Spect_09_09_16_G_32.txt';
% SpectName{3} = 'Jaz2Spect_09_09_16_B_48.txt';
% SpectName{4} = 'Jaz2Spect_09_09_16_K_64.txt';

% Viewpixx EEG: B114
% SpectName{1} = 'JazEEGSpect_09_09_16_R_16.txt';
% SpectName{2} = 'JazEEGSpect_09_09_16_G_32.txt';
% SpectName{3} = 'JazEEGSpect_09_09_16_B_48.txt';
% SpectName{4} = 'JazEEGSpect_09_09_16_K_64.txt';

% For the measurements made 7 Aug 2017, we started each new gun set of spectral measurements with a '1'
% So the max for each gun will be '16'
SpectName{1} = 'JazSpect_07_08_17_R_16.txt';
SpectName{2} = 'JazSpect_07_08_17_G_16.txt';
SpectName{3} = 'JazSpect_07_08_17_B_16.txt';
SpectName{4} = 'JazSpect_07_08_17_K_16.txt';

for ThisGunSet = 1:4
    % Initialize variables.
    fName = fullfile(DataDir, SpectName{ThisGunSet});
    startRow = 18;
    % Open the text file.
    fileID = fopen(fName,'r');
    % Read columns of data according to format string.
    dataArray = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    % Close the text file.
    fclose(fileID);
    
    % Convert the Absolute irradiance values to doubles, & get rid of the last row (contains text)
    % These are found in dataArray{2}
    RawSpectra(:,ThisGunSet) = str2double(dataArray{2}(1:end-1));
end

% Also, grab the wavelengths. These are the same for all spectra.
% They are found in dataArray{1}. Also get rid of final row.
wl = str2double(dataArray{1}(1:end-1));

% *** Resample the spectra. ***
% The spectra are sampled at an irregular stepsize that is less than 1 nm.
% So we want to resample them at a desired range of wavelengths that increase at 1 integer nm.
% We do so using 'interp1'.
% this requires inputs of the current x values (wavelength), the current y values (the spectrum, in abs irradiance),
% and finally, the desired range of x values (desired wavelengths, in nm).
% The output is the interpolated spectra (y values) at the desired range of wavelengths (x values).

desiredwl = ceil(min(wl)):floor(max(wl)); %the desired range of wavelengths in integer 1 nm steps
figure
for jj = 1:4 %across channels
    resampledSpectra(:,jj) = interp1(wl, RawSpectra(:,jj), desiredwl);
    %Plot them to have a look:
    %Plot the raw (original) spectrum:
    plot(wl, RawSpectra(:,jj), '-', 'Color', PlotCols(jj,:));
    hold on
    %Now plot the resampled spectrum on top as circles:
    plot(desiredwl, resampledSpectra(:,jj), 'o',  'Color', PlotCols(jj,:));
end

xlabel('wavelength (nm)')
ylabel('Absolute Irradiance (uW/cm^2/nm)')

% Just plot the resampled spectra, within a sensible range:
figure
for jj = 1:3 %across channels: only R, G & B.
    resampledSpectra(:,jj) = interp1(wl, RawSpectra(:,jj), desiredwl);
    % Plot them to have a look:
    % Plot the raw (original) spectrum:
    plot(wl, RawSpectra(:,jj), '-', 'Color', PlotCols(jj,:), 'LineWidth', 2);
    hold on
end
xlabel('wavelength (nm)')
ylabel('Absolute Irradiance (uW/cm^2/nm)')
xlim([370 730]) % Restrict x-axis to approx range of visible light
%ylim([-0.01 0.045])

% store the wavelengths in the first column of the resampled spectra:
resampledSpectra = [desiredwl', resampledSpectra];
RawSpectra = [wl, RawSpectra];
%Lastly, save the results, including the raw + resampled spectra, and the wl + desired wl.
save(Proc_fName, 'Lv', 'CIE_x', ...
    'CIE_y', 'GammaValues', 'RawSpectra', 'resampledSpectra')





