%% Shimmer Gait Phases
imu_pie5= load('Sujeto7_Nico_Session6_imu_pie_Calibrated_PC');
imup5.timestamp= (imu_pie5.imu_pie_TimestampSync_Unix_CAL - imu_pie5.imu_pie_TimestampSync_Unix_CAL(1))/1000;
imup5.gyroY = imu_pie5.imu_pie_Gyro_Y_CAL*-1;
imup5.gyroZ = imu_pie5.imu_pie_Gyro_Z_CAL;
windowSize = 30; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
imup5.gyroY = filter(b,a,imup5.gyroY);
%figure; plot(imup5.gyroY(1000:20000));
%
clc; close all; clearvars -except imup5

aux_t = imup5.timestamp(5000:15000);
aux_y = imup5.gyroY(5000:15000);

clear imup5

imup5.gyroY = aux_y;
imup5.timestamp = aux_t;

clear aux*

plot(imup5.gyroY)

sw_gain = 3;
j = 1;
k = 1;

max_peak_gain = 50;
distance_toeoff = 70;
[max_negative_peaks,pos_peaks]= findpeaks(-imup5.gyroY,'MinPeakDistance',distance_toeoff,'MinPeakHeight',max_peak_gain);

local_peak_gain = 50;
distance_heelstrike = 70;
[total_negative_peaks,pos_total_peaks]= findpeaks(-imup5.gyroY,'MinPeakDistance',distance_heelstrike,'MinPeakProminence',10);
j=1;
for i=2:length(total_negative_peaks)
    if total_negative_peaks(i) <= total_negative_peaks(i-1)*0.6
        local_negative_peaks(j) = total_negative_peaks(i);
        pos_local_peaks(j) = pos_total_peaks(i);
        j = j+1;
    end
end

%hold on ; plot(imup5.timestamp(pos_local_peaks),-local_negative_peaks,'*')

for i=1:length(imup5.gyroY)
    if imup5.gyroY(i) > sw_gain
        swing_phase.angle(j) = imup5.gyroY(i); 
        swing_phase.time(j) = imup5.timestamp(i);
        phase.id(i) = 3; %%% SWING PHASE
        phase.time(i) = imup5.timestamp(i);
        j = j+1;
    elseif imup5.gyroY(i) <= sw_gain*1.5 && imup5.gyroY(i) >= -sw_gain*4
        stance_phase.angle(k) = imup5.gyroY(i); 
        stance_phase.time(k) = imup5.timestamp(i);
        phase.id(i) = 1; %%% STANCE PHASE
        phase.time(i) = imup5.timestamp(i);
        k = k+1;    
    else
        phase.id(i) = -10;
        phase.time(i) = imup5.timestamp(i);
        for z=1:length(max_negative_peaks)
            if imup5.timestamp(pos_peaks(z)) == imup5.timestamp(i) 
                phase.id(i) = 2; %%% TOE OFF
                phase.time(i) = imup5.timestamp(i);
            end
        end
        for w=1:length(local_negative_peaks)
            if imup5.timestamp(pos_local_peaks(w)) == imup5.timestamp(i) 
                phase.id(i) = 0; %%% HEEL STRIKE
                phase.time(i) = imup5.timestamp(i);
            end
        end
    end
end


phase.id_filtered = fillme(phase.id,1,2,3);
phase.id_filtered = fillme(phase.id_filtered,3,0,1);
count = 0;
try
for i = 2:length(phase.id_filtered)-10
    if phase.id_filtered(i) == -10 
        prev_id = phase.id_filtered(i-1);
        prev_i = i;
        exit = 0;
   
        while exit == 0
           i = i + 1;

           if phase.id_filtered(i) ~= -10
                prev_id = phase.id_filtered(i-1);
                post_i = i;
           
                exit = 1;
           end
        end
        if prev_id == 0 && post_id == 2
            phase.id_filtered(prev_i:post_i-1) = 1;
        elseif prev_id == 1 && post_id == 3
            phase.id_filtered(prev_i:post_i-1) = 2;
        elseif prev_id == 2 && post_id == 0
            phase.id_filtered(prev_i:post_i-1) = 3;
        elseif prev_id == 3 && post_id == 1
            phase.id_filtered(prev_i:post_i-1) = 0;
        else
            phase.id_filtered(prev_i:post_i-1) = 0;
        end
    end
end,
catch
phase.id_filtered(phase.id_filtered == -10) = 0;
end
 yyaxis left
%plot(phase.time,phase.id,'*'); hold on;
plot(phase.time,phase.id_filtered,'ok');
 yyaxis right
% %plot(swing_phase.time,swing_phase.angle)
% %plot(stance_phase.time,stance_phase.angle)
plot(imup5.timestamp,imup5.gyroY)

t0.vector = getPhaseTime(phase.time, phase.id_filtered, 0);
t1.vector = getPhaseTime(phase.time, phase.id_filtered, 1);
t2.vector = getPhaseTime(phase.time, phase.id_filtered, 2);
t3.vector = getPhaseTime(phase.time, phase.id_filtered, 3);

t0.mean = mean(t0.vector);
t1.mean = mean(t1.vector);
t2.mean = mean(t2.vector);
t3.mean = mean(t3.vector);
t0.std = std(t0.vector);
t1.std = std(t1.vector);
t2.std = std(t2.vector);
t3.std = std(t3.vector);
t0.name = 'heel strike';
t1.name = 'flat foot';
t2.name = 'toe off';
t3.name = 'swing phase';
display('Heel strike: ' + string(t0.mean) + ' +- ' + string(t0.std))
display('Flat foot: ' + string(t1.mean) + ' +- ' + string(t1.std))
display('Toe off: ' + string(t2.mean) + ' +- ' + string(t2.std))
display('Swing phase: ' + string(t3.mean) + ' +- ' + string(t3.std))
