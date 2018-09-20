function CalibrateUsingJaz_Display_only

% This code runs through a Psychtoolbox calibration procedure, presenting
% full screen images that vary the luminance intensity of either the Red,
% Green, or Blue guns (RGB), or all guns together (K).
% They are displayed in that order.
%
% The code displays the same colour to both eyes, and will wait until the
% **spacebar** is pressed before proceeding.
%
% The point is to allow measurements of spectra or luminance to be made
% manually elsewhere on a separate machine - nothing is actually recorded by this
% code, only the stimuli are displayed.
%
% This version differs from CalibrateUsingJaz_Vpixx_Display_only because it does not call the Vpixx
% device at all - it can be used on any display.
%
% Measurements could be taken from behind the polarising goggles, but not that they
% would not be in operation and should be permanently in their 'open' state.
%
% R Maloney, April 2016

%set some parameters for the calibration:
nLevels=3; % The number of luminance steps for each channel

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')
%Define response key to abort the program: 'q'
ResponseQuit = KbName('q'); % q to quit
ResponseBegin = KbName('space'); %space bar to begin calibration.


%%%%-----------------------------------------------------%%%%
%       Set up the Vpixx and PTB display
%%%%-----------------------------------------------------%%%%

try
    
    % Set up the Psychtoolbox screen and initialise the Vpixx device, if using one:
    AssertOpenGL;
    PsychImaging('PrepareConfiguration'); % Prepare pipeline for configuration.
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    %Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    
    % Don't want the debugging stuff for present purposes
    Screen( 'Preference', 'Verbosity', 0 );
    %Choose the screen: it is usually the max screen no. available.
    %Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
    %So if we're on that machine, need to -1 from the screen number:
    [~, CompName] = system('hostname'); %find out the computer name
    if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... %normal strcmp not working here, can't figure out why...
            && (length(Screen('Screens'))>1) %and there is more than 1 display connected...
        WhichScreen = max( Screen( 'Screens' ) )-1;
    else
        WhichScreen = max( Screen( 'Screens' ) ); %should be right for any other machine!
    end
    
    % Open an on screen (grey) window and configure the imaging pipeline
    [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
    
    %Raise priority level
    priorityLevel = MaxPriority( win );
    Priority( priorityLevel );
    
    %open a black window in both eyes
    Screen('FillRect', win, 0);
    Screen('Flip', win);
    
    %Display the welcome screen & wait for user input.
    % Select specific text font, style and size:
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 20);
    userStart = false;
    while ~userStart
        [keyIsDown secs, keyCode ] = KbCheck(); %Check for keyboard response
        if  keyCode(ResponseBegin)
            userStart = true; %If the subject presses spacebar, then begin
        elseif keyCode(ResponseQuit) %if they press 'q', quit, close everything down
            %user has aborted
            Priority( 0 ); %restore priority
            Screen( 'CloseAll' ); %close down screen
            ShowCursor; %reveal mouse cursor
            error('Aborting calibration: You quit the program!')
        end
        
        %print welcome screen/instructions to both eyes.
        fprintf('Welcome. Press SPACE to begin calibration, and to proceed through the different luminance levels. \n \nPress ''q'' to quit at any time.'); %also display in command window
        
        DrawFormattedText(win, ...
            'Welcome. Press SPACE to begin calibration, and to proceed through the different luminance levels. \n \nPress ''q'' to quit at any time.', ...
            'center', 'center', 1);
        Screen('Flip', win);
    end
    
    gunSet=[1 0 0; 0 1 0; 0 0 1; 1 1 1]; %which of the RGB channels (guns) to run, ie R, G, B & R+G+B=K channels
    contLevels=linspace(0,1,nLevels); %the intensities for each channel at each step
    
    %Set aside a table of the RGB(K) values that are run through.
    RGBKvals = zeros(4 * nLevels,3); x=0;
    for ii=1:4
        for jj=1:nLevels
            x = x+1;
            RGBKvals(x,:) = gunSet(ii,:).*contLevels(jj);
        end
    end
    
    missedFrames = 0; %keep a tally of any missed presentation deadlines
    
    %Now run through each luminance intensity, for each of the R, G, B and K channels
    
    HideCursor; %hide the mouse cursor so it doesn't interfere with the measurements
    for thisGunSet=1:4
        
        %Start each channel with a black (0 lum) screen (though not measured)
        Screen('FillRect', win, [0 0 0]);
        [~ , ~ , ~, missed] = Screen('Flip',win); %,[],1);
        
        %keep record of any missed frames:
        if missed > 0
            missedFrames = missedFrames + 1;
        end
        
        for thisContLevelIndex=1:nLevels
            
            %determine the current RGB values to send to the screen (somewhere betw. 0-1)
            fillColour = gunSet(thisGunSet,:) * (contLevels(thisContLevelIndex));
            
            %Print info to the Command Window about current channel and luminance level
            fprintf('\nGun %d level %d\n',thisGunSet,thisContLevelIndex);
            
            % fill the screen for as long as necessary.
            % This should continue looping through the colours until the spacebar is pressed.
            WaitSecs(0.2)
            KbCheck(); % take a quick KbCheck to load it now & flush any stored events
            userProceed = false;
            
            while ~userProceed
                
                [~, ~, keyCode] = KbCheck(); %Check for keyboard response
                if  keyCode(ResponseBegin)
                    userProceed = true; %If the subject presses spacebar, then continue
                elseif keyCode(ResponseQuit) %if they press 'q', quit, close everything down
                    %user has aborted
                    Priority( 0 ); %restore priority
                    Screen( 'CloseAll' ); %close down screen
                    ShowCursor; %reveal mouse cursor
                    error('Aborting calibration: You quit the program!')
                end
                
                Screen('FillRect', win, fillColour );
                Screen('DrawingFinished', win);
                [~ , ~ , ~, missed] = Screen('Flip',win); % ,[],1);
                
                %keep record of any missed frames:
                if missed > 0
                    missedFrames = missedFrames + 1;
                end
                
            end %end of stimulus presentation
            
            %abort if the user presses 'q', and tidy up
            [~, ~, keyCode] = KbCheck();
            if keyCode(ResponseQuit)
                %user has aborted
                Priority( 0 ); %restore priority
                
                Screen( 'CloseAll' ); %close down screen
                ShowCursor; %reveal mouse cursor
                error('Aborting calibration: You quit the program!')
            end
        end     % end of loop across different luminance levels of current channel
    end         % end of loop across channels
    
    %This is the point where we might want to check the calibration but I don't think it's necessary for present purposes,
    %especially since we are already running 4 sets of measurements
    
catch MException
    
    rethrow (MException)
    %Shut things down in the event of an error
    ShowCursor; %reveal the mouse cursor again
    
    %Close down the PTB screen
    Screen('CloseAll');
    % We throw the error again so the user sees the error description.
    psychrethrow(psychlasterror);
    Priority( 0 ); %restore priority
    error('error!')
    
end %end of try/catch loop

%Finished! Tidy up.
ShowCursor; %reveal the mouse cursor again
%Close down the PTB screen
Screen('CloseAll');
Priority( 0 ); %restore priority

missedFrames

%the end.