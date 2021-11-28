%% sEMG DATA PRCOESSING

clc; clear all;

% Loading EMG signal file.
load('EMG_data_MAMA.mat');

% Extracting 'event', 'biceps' and 'triceps' data in different variables.
emg_events_data = EMG_data(1,:);
emg_biceps_data = EMG_data(2,:);
emg_triceps_data = EMG_data(3,:);

% Sampling rate in Hz:
sampling_freq = 1000;

% Initializing the axis of the plot (Time(sec))
x1 = 0 : 1/sampling_freq : (length(EMG_data)/sampling_freq) - (1/sampling_freq);

%% Using direct Bandpass command:

% Designing Bandpassfilter using designfilt and filtfilt:
d = designfilt('bandpassfir','StopbandFrequency1',20,'PassbandFrequency1',30,'PassbandFrequency2',450,'StopbandFrequency2',500,'StopbandAttenuation1',60,'PassbandRipple',1,'StopbandAttenuation2',60,'SampleRate',1000);

% Applying the Bandpass filter on biceps and triceps EMG data.
bp2_2 = filtfilt(d,double(emg_biceps_data));
bp3_2 = filtfilt(d,double(emg_triceps_data));

%% Rectifying the filtered data:

%Rectifying the filtered signals using 'abs' command.
rect_bp2_2 = abs(bp2_2);
rect_bp3_2 = abs(bp3_2);

%% Envelop of the muscle signals (low pass 3-6 hz):

% Using 'filtfilt' and butterworth lowpass filter to envelop the signal:
% Cutoff frequency 3 Hz
LP_1 = 6;
NyqFreq = sampling_freq/2;
Wp_1 = LP_1/NyqFreq;
[F_1,E_1] = butter(2,Wp_1,'low');
enve_sig2_1 = filtfilt(F_1,E_1,rect_bp2_2);
enve_sig3_1 = filtfilt(F_1,E_1,rect_bp3_2);

%% Downsampling:

% Downsampling factor:  resample command
fact = 2;
y1 = downsample (emg_events_data, fact);
y2 = downsample (enve_sig2_1, fact);
y3 = downsample (enve_sig3_1, fact);

%% Extracting Event 'CUE' and 'GO' Index points from sEMG data:

% Generated Time data (sec):
emg_time_data = x1;

% Flag variables used to detect the occurance of an event
a_emg_1 = true;
a_emg_2 = true;

% Counter varibles used to count the occurance of an event
epochs_emg_1 = 0;
epochs_emg_2 = 0;

% Using for loop to extract the 'starting' and 'ending' indexes of the desired
% sets of movements. In our case that would be 1,6,10 and 11.

for k1 = 1:length(emg_events_data)
    
    % For event 1 CUE: 
    if( emg_events_data(k1) == 1 && a_emg_1 == true)
        epochs_emg_1 = epochs_emg_1 + 1;
        a_emg_1 = false;
       
        if(epochs_emg_1 == 1)
            epochs_emg_1_start_set1 = k1;
        end

        if(epochs_emg_1 == 96)
            epochs_emg_1_end_set1 = k1;
        end
     
        if(epochs_emg_1 == 481)
            epochs_emg_1_start_set6 = k1;
        end

        if(epochs_emg_1 == 576)
            epochs_emg_1_end_set6 = k1;
        end
    
        if(epochs_emg_1 == 865)
            epochs_emg_1_start_set10 = k1;
        end

        if(epochs_emg_1 == 960)
            epochs_emg_1_end_set10 = k1;
        end

        if(epochs_emg_1 == 961)
            epochs_emg_1_start_set11 = k1;
        end

        if(epochs_emg_1 == 1056)
            epochs_emg_1_end_set11 = k1;
        end
    end
     
    if(emg_events_data(k1) ~= 1 && a_emg_1 == false) 
        a_emg_1 = true;
    end
   
    % For event 2 GO: 
    if( emg_events_data(k1) == 2 && a_emg_2 == true)
        epochs_emg_2 = epochs_emg_2 + 1;
        a_emg_2 = false;
            
        if(epochs_emg_2 == 1)
            epochs_emg_2_start_set1 = k1;
        end
    
        if(epochs_emg_2 == 96)
            epochs_emg_2_end_set1 = k1;
        end
         
        if(epochs_emg_2 == 481)
            epochs_emg_2_start_set6 = k1;
        end
    
        if(epochs_emg_2 == 576)
            epochs_emg_2_end_set6 = k1;
        end    
        
        if(epochs_emg_2 == 865)
            epochs_emg_2_start_set10 = k1;
        end
    
        if(epochs_emg_2 == 960)
            epochs_emg_2_end_set10 = k1;
        end
    
        if(epochs_emg_2 == 961)
            epochs_emg_2_start_set11 = k1;
        end
    
        if(epochs_emg_2 == 1056)
            epochs_emg_2_end_set11 = k1;
        end
    end
     
    if(emg_events_data(k1) ~= 2 && a_emg_2 == false) 
        a_emg_2 = true;
    end
end

%% Extracting sEMG row signal with on top the filtered signal:

%Extracting set1 emg raw signal:
emg_raw_events_set1 = emg_events_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_raw_bicep_set1 =  emg_biceps_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_raw_tricep_set1 = emg_triceps_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_raw_time_set1 = x1(1:length(emg_raw_events_set1));

%Extracting set1 emg filtered signal:
emg_filt_bicep_set1 =  bp2_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_filt_tricep_set1 = bp3_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);

%Extracting set6 emg signal:
emg_raw_events_set6 = emg_events_data(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_raw_bicep_set6 =  emg_biceps_data(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_raw_tricep_set6 = emg_triceps_data(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_raw_time_set6 = x1(1:length(emg_raw_events_set6));

%Extracting set6 emg filtered signal:
emg_filt_bicep_set6 =  bp2_2(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_filt_tricep_set6 = bp3_2(epochs_emg_1_start_set1:epochs_emg_2_end_set6);

%Extracting set10 emg raw signal:
emg_raw_events_set10 = emg_events_data(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_raw_bicep_set10 =  emg_biceps_data(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_raw_tricep_set10 = emg_triceps_data(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_raw_time_set10 = x1(1:length(emg_raw_events_set10));

%Extracting set10 emg filtered signal:
emg_filt_bicep_set10 =  bp2_2(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_filt_tricep_set10 = bp3_2(epochs_emg_1_start_set10:epochs_emg_2_end_set10);

%Extracting set11 emg raw signal:
emg_raw_events_set11 = emg_events_data(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_raw_bicep_set11 =  emg_biceps_data(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_raw_tricep_set11 = emg_triceps_data(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_raw_time_set11 = x1(1:length(emg_raw_events_set11));

%Extracting set11 emg filtered signal:
emg_filt_bicep_set11 =  bp2_2(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_filt_tricep_set11 = bp3_2(epochs_emg_1_start_set11:epochs_emg_2_end_set11);

%% Plotting sEMG row signal with on top the filtered signal:

% Plotting subplot of EMG Biceps raw and filtered signal for all 4 sets:

figure (10);
subplot(4,1,1);
plot(emg_raw_time_set1, emg_raw_bicep_set1, 'r', emg_raw_time_set1, emg_filt_bicep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set1");
xlim([0 ((length(emg_raw_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,2); 
plot(emg_raw_time_set6, emg_raw_bicep_set6, 'r', emg_raw_time_set6, emg_filt_bicep_set6,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set6");
xlim([0 ((length(emg_raw_time_set6)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,3); 
plot(emg_raw_time_set10, emg_raw_bicep_set10, 'r', emg_raw_time_set10, emg_filt_bicep_set10,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set10");
xlim([0 ((length(emg_raw_time_set10)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,4); 
plot(emg_raw_time_set11, emg_raw_bicep_set11, 'r', emg_raw_time_set11, emg_filt_bicep_set11,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set11");
xlim([0 ((length(emg_raw_time_set11)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');


% Plotting subplot of EMG Triceps raw and filtered signal for all 4 sets:
figure (11);
subplot(4,1,1);
plot(emg_raw_time_set1, emg_raw_tricep_set1, 'r', emg_raw_time_set1, emg_filt_bicep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set1");
xlim([0 ((length(emg_raw_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,2); 
plot(emg_raw_time_set6, emg_raw_tricep_set6, 'r', emg_raw_time_set6, emg_filt_bicep_set6,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set6");
xlim([0 ((length(emg_raw_time_set6)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,3); 
plot(emg_raw_time_set10, emg_raw_tricep_set10, 'r', emg_raw_time_set10, emg_filt_bicep_set10,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set10");
xlim([0 ((length(emg_raw_time_set10)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,4); 
plot(emg_raw_time_set11, emg_raw_tricep_set11, 'r', emg_raw_time_set11, emg_filt_bicep_set11,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set11");
xlim([0 ((length(emg_raw_time_set11)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

%%  Extracting sEMG rectified data with on top the Envelope:

%Extracting set1 emg rectified signal:
emg_rect_events_set1 = emg_events_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_rect_bicep_set1 =  rect_bp2_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_rect_tricep_set1 = rect_bp3_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_rect_time_set1 = x1(1:length(emg_rect_events_set1));

%Extracting set1 emg enveloped signal:
emg_envelop_bicep_set1 =  enve_sig2_1(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_envelop_tricep_set1 = enve_sig3_1(epochs_emg_1_start_set1:epochs_emg_2_end_set1);

%Extracting set6 emg rectified signal:
emg_rect_events_set6 = emg_events_data(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_rect_bicep_set6 =  rect_bp2_2(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_rect_tricep_set6 = rect_bp3_2(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_rect_time_set6 = x1(1:length(emg_rect_events_set6));

%Extracting set6 emg enveloped signal:
emg_envelop_bicep_set6 =  enve_sig2_1(epochs_emg_1_start_set6:epochs_emg_2_end_set6);
emg_envelop_tricep_set6 = enve_sig3_1(epochs_emg_1_start_set6:epochs_emg_2_end_set6);

%Extracting set10 emg rectified signal:
emg_rect_events_set10 = emg_events_data(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_rect_bicep_set10 =  rect_bp2_2(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_rect_tricep_set10 = rect_bp3_2(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_rect_time_set10 = x1(1:length(emg_rect_events_set10));

%Extracting set10 emg enveloped signal:
emg_envelop_bicep_set10 =  enve_sig2_1(epochs_emg_1_start_set10:epochs_emg_2_end_set10);
emg_envelop_tricep_set10 = enve_sig3_1(epochs_emg_1_start_set10:epochs_emg_2_end_set10);

%Extracting set11 emg rectified signal:
emg_rect_events_set11 = emg_events_data(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_rect_bicep_set11 =  rect_bp2_2(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_rect_tricep_set11 = rect_bp3_2(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_rect_time_set11 = x1(1:length(emg_rect_events_set11));

%Extracting set11 emg enveloped signal:
emg_envelop_bicep_set11 =  enve_sig2_1(epochs_emg_1_start_set11:epochs_emg_2_end_set11);
emg_envelop_tricep_set11 = enve_sig3_1(epochs_emg_1_start_set11:epochs_emg_2_end_set11);

%% Plotting sEMG rectified data with on top the Envelope:

% Plotting subplot of EMG Biceps rectified and enveloped signal for all 4 sets:
figure (12);
subplot(4,1,1);
plot(emg_rect_time_set1, emg_rect_bicep_set1, 'r', emg_rect_time_set1, emg_envelop_bicep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set1");
xlim([0 ((length(emg_rect_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,2); 
plot(emg_rect_time_set6, emg_rect_bicep_set6, 'r', emg_rect_time_set6, emg_envelop_bicep_set6,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set6");
xlim([0 ((length(emg_rect_time_set6)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,3); 
plot(emg_rect_time_set10, emg_rect_bicep_set10, 'r', emg_rect_time_set10, emg_envelop_bicep_set10,'b');  ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set10");
xlim([0 ((length(emg_rect_time_set10)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,4); 
plot(emg_rect_time_set11, emg_rect_bicep_set11, 'r', emg_rect_time_set11, emg_envelop_bicep_set11,'b');  ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set11");
xlim([0 ((length(emg_rect_time_set11)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

% Plotting subplot of EMG Triceps rectified and enveloped signal for all 4 sets:
figure (13);
subplot(4,1,1);
plot(emg_rect_time_set1, emg_rect_tricep_set1, 'r', emg_rect_time_set1, emg_envelop_tricep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1");
xlim([0 ((length(emg_rect_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,2); 
plot(emg_rect_time_set6, emg_rect_tricep_set6, 'r', emg_rect_time_set6, emg_envelop_tricep_set6,'b');  ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1 Set6");
xlim([0 ((length(emg_rect_time_set6)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,3); 
plot(emg_rect_time_set10, emg_rect_tricep_set10, 'r', emg_rect_time_set10, emg_envelop_tricep_set10,'b');  ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1 Set10");
xlim([0 ((length(emg_rect_time_set10)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,4); 
plot(emg_rect_time_set11, emg_rect_tricep_set11, 'r', emg_rect_time_set11, emg_envelop_tricep_set11,'b');  ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1 Set11");
xlim([0 ((length(emg_rect_time_set11)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

%% MOVEMENT DATA LOADING 

%Loading movement data
load('kinem_MAMA.mat');

% Extracting 'time points', 'events', 'x-y cursor' and 'x-y target' data into different variables.
time_points = kinem_data(1,:);
raw_events = kinem_data(2,:);
x_cursor = kinem_data(3,:);
y_cursor = kinem_data(4,:);
x_target = kinem_data(5,:);
y_target = kinem_data(6,:);

%Sampling frequency in Hz
sampling_freq2 = 100;

% Generated Time data (sec):
gen_time_data = 0 : 1/sampling_freq2 : (length(time_points)/sampling_freq2) - (1/sampling_freq2);

%% Extracting Event 'CUE' and 'GO' Index points from movement data:

% Flag variable used to detect the occurance of an event
a = true; 

% Counter varibles used to count the occurance of an event
epochs = 0; 

% Again using for loop to extract the 'starting' and 'ending' indexes of the desired
% sets of movements but for movement data.

for k1 = 1:length(raw_events)
    
    % For event '2' CUE:
    if( raw_events(k1) == 2 && a == true)
        epochs = epochs + 1;
        a = false;
        if(epochs == 1)
            epoch_start_set1 = k1;
        end

        if(epochs == 96)
            epoch_end_set1 = k1;
        end
     
        if(epochs == 481)
            epoch_start_set6 = k1;
        end

        if(epochs == 576)
            epoch_end_set6 = k1;
        end
    
        if(epochs == 865)
            epoch_start_set10 = k1;
        end

        if(epochs == 960)
            epoch_end_set10 = k1;
        end

        if(epochs == 961)
            epoch_start_set11 = k1;
        end

        if(epochs == 1056)
            epoch_end_set11 = k1;
        end
    end
     
    if(raw_events(k1) ~= 2 && a == false) 
        a = true;
    end
end

%% Extracting XY MOVEMENT and TARGETS DATA:

%Extracting set1 movement data:
movement_events_set1 = raw_events(203:36200);
movementx_set1 = x_cursor(203:36200);
movementy_set1 = y_cursor(203:36200);
target_x_set1 = x_target(203:36200);
target_y_set1 = y_target(203:36200);
gen_time_data_set1 = gen_time_data(1:length(movement_events_set1));

%Extracting set6 movement data:
movement_events_set6 = raw_events(174002:209200);
movementx_set6 = x_cursor(174002:209200);
movementy_set6 = y_cursor(174002:209200);
target_x_set6 = x_target(174002:209200);
target_y_set6 = y_target(174002:209200);
gen_time_data_set6 = gen_time_data(1:length(movement_events_set6));

%Extracting set10 movement data:
movement_events_set10 = raw_events(314402:349000);
movementx_set10 = x_cursor(314402:349000);
movementy_set10 = y_cursor(314402:349000);
target_x_set10 = x_target(314402:349000);
target_y_set10 = y_target(314402:349000);
gen_time_data_set10 = gen_time_data(1:length(movement_events_set10));

%Extracting set11 movement data:
movement_events_set11 = raw_events(349203:384200);
movementx_set11 = x_cursor(349203:384200);
movementy_set11 = y_cursor(349203:384200);
target_x_set11 = x_target(349203:384200);
target_y_set11 = y_target(349203:384200);
gen_time_data_set11 = gen_time_data(1:length(movement_events_set11));

%% Plotting xy movement data in time(sec):

% Plotting X-Y movement data in time(sec) for all 4 sets:
figure (14);
subplot(4,1,1);
plot(gen_time_data_set1,movementx_set1,'r',gen_time_data_set1, movementy_set1,'b'); ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set1");
xlim([0 ((length(movementx_set1)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');

subplot(4,1,2); 
plot(gen_time_data_set6,movementx_set6,'r',gen_time_data_set6, movementy_set6,'b');  ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set6");
xlim([0 ((length(movementx_set6)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');

subplot(4,1,3); 
plot(gen_time_data_set10,movementx_set10,'r',gen_time_data_set10, movementy_set10,'b'); ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set10");
xlim([0 ((length(movementx_set10)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');

subplot(4,1,4); 
plot(gen_time_data_set11,movementx_set11,'r',gen_time_data_set11, movementy_set11,'b'); ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set11");
xlim([0 ((length(movementx_set11)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');
%% Plotting for XY movement data and targets data:

% Plotting X-Y movement and X-Y Tagret data for all 4 sets:
figure (18);
subplot(4,1,1);
plot(gen_time_data_set1,movementx_set1,'r',gen_time_data_set1, movementy_set1,'b',gen_time_data_set1, target_x_set1,'black',gen_time_data_set1, target_y_set1,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set1");
xlim([0 ((length(movementx_set1)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');

subplot(4,1,2); 
plot(gen_time_data_set6,movementx_set6,'r',gen_time_data_set6, movementy_set6,'b',gen_time_data_set6, target_x_set6,'black',gen_time_data_set6, target_y_set6,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set6");
xlim([0 ((length(movementx_set6)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');

subplot(4,1,3); 
plot(gen_time_data_set10,movementx_set10,'r',gen_time_data_set10, movementy_set10,'b',gen_time_data_set10, target_x_set10,'black',gen_time_data_set10, target_y_set10,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set10");
xlim([0 ((length(movementx_set10)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');

subplot(4,1,4); 
plot(gen_time_data_set11,movementx_set11,'r',gen_time_data_set11, movementy_set11,'b',gen_time_data_set11, target_x_set11,'black',gen_time_data_set11, target_y_set11,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set11");
xlim([0 ((length(movementx_set11)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');
