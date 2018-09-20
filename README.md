# wadelab-calibration
Stimulus code, data and analysis code used for calibrating the Viewpixx

CalibrateUsingJaz_Display_only.m% This displays a full screen colour patch, cycling through Red, Green, Blue phosphors and then 
all channels (monochrome or K), in steps you can specify (eg. 8, 16, etc)
The colour screens can be advanced by pressing the keyboard.
It only displays the stimulus, and does not communicate with the spectrophotometer

CalibrateColour_Viewpixx_no_goggles.m % This will load, plot and process the spectral data for calibration
GammaPower.m % needed for the above; fits the gamma function to the data to obtain the gamma exponents for calibration

ViewpixxEEG_Processed_cal_data_7_8_2017.mat	% The data: Viewpixx EEG (Room B114) obtained 7 Aug 2017 without the 3D goggles
ViewpixxEEG_calibration_7_8_17_no_goggles.pdf % example of the plotted data
