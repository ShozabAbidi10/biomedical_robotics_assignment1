%% sEMG DATA PRCOESSING

clc; clear all;

% Loading EMG signal file.

load('EMG_data_MAMA.mat');
emg_events_data = EMG_data(1,:);
emg_biceps_data = EMG_data(2,:);
emg_triceps_data = EMG_data(3,:);

% Sampling rate in Hz:
sampling_freq = 1000;

x1 = 0 : 1/sampling_freq : (length(EMG_data)/sampling_freq) - (1/sampling_freq);

%% Using direct Bandpass command:

% Designing Bandpassfilter using designfilt and filtfilt:
d = designfilt('bandpassfir','StopbandFrequency1',20,'PassbandFrequency1',30,'PassbandFrequency2',450,'StopbandFrequency2',500,'StopbandAttenuation1',60,'PassbandRipple',1,'StopbandAttenuation2',60,'SampleRate',1000);
bp2_2 = filtfilt(d,double(emg_biceps_data));
bp3_2 = filtfilt(d,double(emg_triceps_data));

%% Rectifying the filtered data:

rect_bp2_2 = abs(bp2_2);
rect_bp3_2 = abs(bp3_2);

%% Envelop of the muscle signals (low pass 3-6 hz):

% 3 Hz Cut off frequency:
LP_1 = 3;
NyqFreq = sampling_freq/2;
Wp_1 = LP_1/NyqFreq;
[F_1,E_1] = butter(2,Wp_1,'low');
enve_sig2_1 = filtfilt(F_1,E_1,rect_bp2_2);
enve_sig3_1 = filtfilt(F_1,E_1,rect_bp3_2);

% 6 Hz Cut off frequency:
% LP_2 = 6; 
% NyqFreq = sampling_freq/2;
% Wp_2 = LP_2/NyqFreq;
% [F_2,E_2] = butter(2,Wp_2,'low');
% enve_sig2_2 = filtfilt(F_2,E_2,rect_bp2_2);
% enve_sig3_2 = filtfilt(F_2,E_2,rect_bp3_2);

%% Downsampling:

% What should the downsampling factor? 
% Ans: this is something you will have to figure out,knowing the
% sampling frequency of the EMG and the one of the movement data

fact = 1;
y1 = downsample (EMG_data(1,:),fact);
y2 = downsample (enve_sig2_1, fact);
y3 = downsample (enve_sig3_1, fact);

%% Extracting Event 'CUE' and 'GO' Index points from sEMG data:

% Generated Time data (sec):
emg_time_data = downsample (x1,fact);

a_emg_1 = true;
a_emg_2 = true;
epochs_emg_1 = 0;
epochs_emg_2 = 0;

for k1 = 1:length(y1)
    % For event = 1: 
    if( y1(k1) == 1 && a_emg_1 == true)
        epochs_emg_1 = epochs_emg_1 + 1;
        a_emg_1 = false;
        
        if(epochs_emg_1 == 1)
            epochs_emg_1_start_set1 = k1;
        end

        if(epochs_emg_1 == 96)
            epochs_emg_1_end_set1 = k1;
        end
     
        if(epochs_emg_1 == 481)
            epochs_emg_1_start_set2 = k1;
        end

        if(epochs_emg_1 == 576)
            epochs_emg_1_end_set2 = k1;
        end
    
        if(epochs_emg_1 == 865)
            epochs_emg_1_start_set3 = k1;
        end

        if(epochs_emg_1 == 960)
            epochs_emg_1_end_set3 = k1;
        end

        if(epochs_emg_1 == 961)
            epochs_emg_1_start_set4 = k1;
        end

        if(epochs_emg_1 == 1056)
            epochs_emg_1_end_set4 = k1;
        end
    end
     
    if(y1(k1) ~= 1 && a_emg_1 == false) 
        a_emg_1 = true;
    end
   
    % For event = 2: 
    if( y1(k1) == 2 && a_emg_2 == true)
        epochs_emg_2 = epochs_emg_2 + 1;
        a_emg_2 = false;
            
        if(epochs_emg_2 == 1)
            epochs_emg_2_start_set1 = k1;
        end
    
        if(epochs_emg_2 == 96)
            epochs_emg_2_end_set1 = k1;
        end
         
        if(epochs_emg_2 == 481)
            epochs_emg_2_start_set2 = k1;
        end
    
        if(epochs_emg_2 == 576)
            epochs_emg_2_end_set2 = k1;
        end    
        
        if(epochs_emg_2 == 865)
            epochs_emg_2_start_set3 = k1;
        end
    
        if(epochs_emg_2 == 960)
            epochs_emg_2_end_set3 = k1;
        end
    
        if(epochs_emg_2 == 961)
            epochs_emg_2_start_set4 = k1;
        end
    
        if(epochs_emg_2 == 1056)
            epochs_emg_2_end_set4 = k1;
        end
    end
     
    if(y1(k1) ~= 2 && a_emg_2 == false) 
        a_emg_2 = true;
    end
end

%% Extracting sEMG row signal with on top the filtered signal:

% Set_1 Raw Signal:
emg_raw_events_set1 = emg_events_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_raw_bicep_set1 =  emg_biceps_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_raw_tricep_set1 = emg_triceps_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_raw_time_set1 = x1(1:length(emg_raw_events_set1));

% Set_1 Filtered Signal:
emg_filt_bicep_set1 =  bp2_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_filt_tricep_set1 = bp3_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);

% Set_2 Raw Signal:
emg_raw_events_set2 = emg_events_data(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_raw_bicep_set2 =  emg_biceps_data(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_raw_tricep_set2 = emg_triceps_data(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_raw_time_set2 = x1(1:length(emg_raw_events_set2));

% Set_2 Filtered Signal:
emg_filt_bicep_set2 =  bp2_2(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_filt_tricep_set2 = bp3_2(epochs_emg_1_start_set1:epochs_emg_2_end_set2);

% Set_3 Raw Signal:
emg_raw_events_set3 = emg_events_data(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_raw_bicep_set3 =  emg_biceps_data(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_raw_tricep_set3 = emg_triceps_data(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_raw_time_set3 = x1(1:length(emg_raw_events_set3));

% Set_3 Filtered Signal:
emg_filt_bicep_set3 =  bp2_2(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_filt_tricep_set3 = bp3_2(epochs_emg_1_start_set3:epochs_emg_2_end_set3);

% Set_4 Raw Signal:
emg_raw_events_set4 = emg_events_data(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_raw_bicep_set4 =  emg_biceps_data(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_raw_tricep_set4 = emg_triceps_data(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_raw_time_set4 = x1(1:length(emg_raw_events_set4));

% Set_4 Filtered Signal:
emg_filt_bicep_set4 =  bp2_2(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_filt_tricep_set4 = bp3_2(epochs_emg_1_start_set4:epochs_emg_2_end_set4);

%% Plotting sEMG row signal with on top the filtered signal:

% Subplot Biceps:

figure (10);
subplot(4,1,1);
plot(emg_raw_time_set1, emg_raw_bicep_set1, 'r', emg_raw_time_set1, emg_filt_bicep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set1");
xlim([0 ((length(emg_raw_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,2); 
plot(emg_time_set2, emg_raw_bicep_set2, 'r', emg_time_set2, emg_filt_bicep_set2,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set6");
xlim([0 ((length(emg_time_set2)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,3); 
plot(emg_time_set3, emg_raw_bicep_set3, 'r', emg_time_set3, emg_filt_bicep_set3,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set10");
xlim([0 ((length(emg_time_set3)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,4); 
plot(emg_time_set4, emg_raw_bicep_set4, 'r', emg_time_set4, emg_filt_bicep_set4,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Signal with on top the Filtered Signal Set11");
xlim([0 ((length(emg_time_set4)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');


% Subplot Tricepts:

figure (11);
subplot(4,1,1);
plot(emg_time_set1, emg_raw_tricep_set1, 'r', emg_time_set1, emg_filt_bicep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set1");
xlim([0 ((length(emg_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,2); 
plot(emg_time_set2, emg_raw_tricep_set2, 'r', emg_time_set2, emg_filt_bicep_set2,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set6");
xlim([0 ((length(emg_time_set2)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,3); 
plot(emg_time_set3, emg_raw_tricep_set3, 'r', emg_time_set3, emg_filt_bicep_set3,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set10");
xlim([0 ((length(emg_time_set3)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

subplot(4,1,4); 
plot(emg_time_set4, emg_raw_tricep_set4, 'r', emg_time_set4, emg_filt_bicep_set4,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Signal with on top the Filtered Signal Set11");
xlim([0 ((length(emg_time_set4)/sampling_freq)-(1/sampling_freq))]);
legend('Raw Data','Filtered Data');

%%  Extracting sEMG rectified data with on top the Envelope:

% Set_1 rectified Signal:
emg_rect_events_set1 = emg_events_data(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_rect_bicep_set1 =  rect_bp2_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_rect_tricep_set1 = rect_bp3_2(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_rect_time_set1 = x1(1:length(emg_rect_events_set1));

% Set_1 enveloped Signal:
emg_envelop_bicep_set1 =  enve_sig2_1(epochs_emg_1_start_set1:epochs_emg_2_end_set1);
emg_envelop_tricep_set1 = enve_sig3_1(epochs_emg_1_start_set1:epochs_emg_2_end_set1);

% Set_2 rectified Signal:
emg_rect_events_set2 = emg_events_data(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_rect_bicep_set2 =  rect_bp2_2(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_rect_tricep_set2 = rect_bp3_2(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_rect_time_set2 = x1(1:length(emg_rect_events_set2));

% Set_2 enveloped Signal:
emg_envelop_bicep_set2 =  enve_sig2_1(epochs_emg_1_start_set2:epochs_emg_2_end_set2);
emg_envelop_tricep_set2 = enve_sig3_1(epochs_emg_1_start_set2:epochs_emg_2_end_set2);

% Set_3 rectified Signal:
emg_rect_events_set3 = emg_events_data(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_rect_bicep_set3 =  rect_bp2_2(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_rect_tricep_set3 = rect_bp3_2(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_rect_time_set3 = x1(1:length(emg_rect_events_set3));

% Set_3 enveloped Signal:
emg_envelop_bicep_set3 =  enve_sig2_1(epochs_emg_1_start_set3:epochs_emg_2_end_set3);
emg_envelop_tricep_set3 = enve_sig3_1(epochs_emg_1_start_set3:epochs_emg_2_end_set3);

% Set_4 rectified Signal:
emg_rect_events_set4 = emg_events_data(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_rect_bicep_set4 =  rect_bp2_2(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_rect_tricep_set4 = rect_bp3_2(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_rect_time_set4 = x1(1:length(emg_rect_events_set4));

% Set_4 enveloped Signal:
emg_envelop_bicep_set4 =  enve_sig2_1(epochs_emg_1_start_set4:epochs_emg_2_end_set4);
emg_envelop_tricep_set4 = enve_sig3_1(epochs_emg_1_start_set4:epochs_emg_2_end_set4);

%% Plotting sEMG rectified data with on top the Envelope:

% Subplot Biceps:

figure (12);
subplot(4,1,1);
plot(emg_rect_time_set1, emg_rect_bicep_set1, 'r', emg_rect_time_set1, emg_envelop_bicep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set1");
xlim([0 ((length(emg_rect_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,2); 
plot(emg_rect_time_set2, emg_rect_bicep_set2, 'r', emg_rect_time_set2, emg_envelop_bicep_set2,'b'); ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set6");
xlim([0 ((length(emg_rect_time_set2)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,3); 
plot(emg_rect_time_set3, emg_rect_bicep_set3, 'r', emg_rect_time_set3, emg_envelop_bicep_set3,'b');  ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set10");
xlim([0 ((length(emg_rect_time_set3)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,4); 
plot(emg_rect_time_set4, emg_rect_bicep_set4, 'r', emg_rect_time_set4, emg_envelop_bicep_set4,'b');  ylabel('Amplitude (units)'); title("sEMG Biceps Rectified Data vs Envelope Data Set11");
xlim([0 ((length(emg_rect_time_set4)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

% Subplot Tricepts:

figure (13);
subplot(4,1,1);
plot(emg_rect_time_set1, emg_rect_tricep_set1, 'r', emg_rect_time_set1, emg_envelop_tricep_set1,'b'); ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1");
xlim([0 ((length(emg_time_set1)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,2); 
plot(emg_rect_time_set2, emg_rect_tricep_set2, 'r', emg_rect_time_set2, emg_envelop_tricep_set2,'b');  ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1 Set6");
xlim([0 ((length(emg_time_set2)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,3); 
plot(emg_rect_time_set3, emg_rect_tricep_set3, 'r', emg_rect_time_set3, emg_envelop_tricep_set3,'b');  ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1 Set10");
xlim([0 ((length(emg_rect_time_set3)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

subplot(4,1,4); 
plot(emg_rect_time_set4, emg_rect_tricep_set4, 'r', emg_rect_time_set4, emg_envelop_tricep_set4,'b');  ylabel('Amplitude (units)'); title("sEMG Triceps Rectified Data vs Envelope Data Set1 Set11");
xlim([0 ((length(emg_rect_time_set4)/sampling_freq)-(1/sampling_freq))]);
legend('Rectified Data','Enveloped Data');

%% MOVEMENT DATA LOADING 

load('kinem_MAMA.mat');

time_points = kinem_data(1,:);
raw_events = kinem_data(2,:);
x_cursor = kinem_data(3,:);
y_cursor = kinem_data(4,:);
x_target = kinem_data(5,:);
y_target = kinem_data(6,:);
sampling_freq2 = 100;
tt = time_points-142.0620; 

% Generated Time data (sec):
gen_time_data = 0 : 1/sampling_freq2 : (length(time_points)/sampling_freq2) - (1/sampling_freq2);

%% Extracting Event 'CUE' and 'GO' Index points from movement data:

a = true; 
epochs = 0; 

for k1 = 1:length(raw_events)
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
            epoch_start_set2 = k1;
        end

        if(epochs == 576)
            epoch_end_set2 = k1;
        end
    
        if(epochs == 865)
            epoch_start_set3 = k1;
        end

        if(epochs == 960)
            epoch_end_set3 = k1;
        end

        if(epochs == 961)
            epoch_start_set4 = k1;
        end

        if(epochs == 1056)
            epoch_end_set4 = k1;
        end
    end
     
    if(raw_events(k1) ~= 2 && a == false) 
        a = true;
    end
end

%% Extracting XY MOVEMENT and TARGETS DATA:

% Set_1
movement_events_set1 = raw_events(203:36200);
movementx_set1 = x_cursor(203:36200);
movementy_set1 = y_cursor(203:36200);
target_x_set1 = x_target(203:36200);
target_y_set1 = y_target(203:36200);
gen_time_data_set1 = gen_time_data(1:length(movement_events_set1));

% Set_2
movement_events_set2 = raw_events(174002:209200);
movementx_set2 = x_cursor(174002:209200);
movementy_set2 = y_cursor(174002:209200);
target_x_set2 = x_target(174002:209200);
target_y_set2 = y_target(174002:209200);
gen_time_data_set2 = gen_time_data(1:length(movement_events_set2));

% Set_3
movement_events_set3 = raw_events(314402:349000);
movementx_set3 = x_cursor(314402:349000);
movementy_set3 = y_cursor(314402:349000);
target_x_set3 = x_target(314402:349000);
target_y_set3 = y_target(314402:349000);
gen_time_data_set3 = gen_time_data(1:length(movement_events_set3));

% Set_4
movement_events_set4 = raw_events(349203:384200);
movementx_set4 = x_cursor(349203:384200);
movementy_set4 = y_cursor(349203:384200);
target_x_set4 = x_target(349203:384200);
target_y_set4 = y_target(349203:384200);
gen_time_data_set4 = gen_time_data(1:length(movement_events_set4));

%% Plotting xy movement data in time(sec):
% Subplot 

figure (14);
subplot(4,1,1);
plot(gen_time_data_set1,movementx_set1,'r',gen_time_data_set1, movementy_set1,'b'); ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set1");
xlim([0 ((length(movementx_set1)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');

subplot(4,1,2); 
plot(gen_time_data_set2,movementx_set2,'r',gen_time_data_set2, movementy_set2,'b');  ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set6");
xlim([0 ((length(movementx_set2)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');

subplot(4,1,3); 
plot(gen_time_data_set3,movementx_set3,'r',gen_time_data_set3, movementy_set3,'b'); ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set10");
xlim([0 ((length(movementx_set3)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');

subplot(4,1,4); 
plot(gen_time_data_set4,movementx_set4,'r',gen_time_data_set4, movementy_set4,'b'); ylabel('Amplitude (units)');title("movement signal X and Y in time (sec) Set11");
xlim([0 ((length(movementx_set4)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY');
%% Plotting for XY movement data and targets data:
% Subplot 

figure (18);
subplot(4,1,1);
plot(gen_time_data_set1,movementx_set1,'r',gen_time_data_set1, movementy_set1,'b',gen_time_data_set1, target_x_set1,'black',gen_time_data_set1, target_y_set1,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set1");
xlim([0 ((length(movementx_set1)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');

subplot(4,1,2); 
plot(gen_time_data_set2,movementx_set2,'r',gen_time_data_set2, movementy_set2,'b',gen_time_data_set2, target_x_set2,'black',gen_time_data_set2, target_y_set2,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set6");
xlim([0 ((length(movementx_set2)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');

subplot(4,1,3); 
plot(gen_time_data_set3,movementx_set3,'r',gen_time_data_set3, movementy_set3,'b',gen_time_data_set3, target_x_set3,'black',gen_time_data_set3, target_y_set3,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set10");
xlim([0 ((length(movementx_set3)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');

subplot(4,1,4); 
plot(gen_time_data_set4,movementx_set4,'r',gen_time_data_set4, movementy_set4,'b',gen_time_data_set4, target_x_set4,'black',gen_time_data_set4, target_y_set4,'green'); xlabel('Time (unit seconds)'); ylabel('Amplitude (units)');title("movement signals X and Y with Targets Set11");
xlim([0 ((length(movementx_set4)/sampling_freq2)-(1/sampling_freq2))]);
legend('MovementX','MovementY','TargetX','TargetY');
