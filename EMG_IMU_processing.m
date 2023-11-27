clear all; close all; clc;

path = "C:\Users\Sophia\OneDrive - ESCUELA COLOMBIANA DE INGENIERIA JULIO GARAVITO\Estudio_Tesis\";
prefix = ["Sujet1_Andy2_", "Sujet4_Mate2_", "Sujet5_Juan2_", "Sujet6_Serg2_", "Sujet7_Nico3_", "Sujet8_Gerar_", "Sujet9_Carlo_", "Sujet10_Seb_", "Sujet11_Juam_", "Sujet12_Serg_"];
folders = ["WOE", "WE", "WET"];
%folders= ["WS1","WS2","WS3"];
muscles = ["Session2_biceps_vasR", "Session4_tibial_gasR", "Session6_tibial_gasR", "Session8_biceps_vasR" ];
muscles2= ["Session1_biceps_vasL", "Session3_tibial_gasL", "Session5_tibial_gasL", "Session7_biceps_vasL" ];
channels = ["biceps_vasR_EMG_CH2_BandPass_Filter_CAL", "tibial_gasR_EMG_CH2_BandPass_Filter_CAL", "tibial_gasR_EMG_CH1_BandPass_Filter_CAL", "biceps_vasR_EMG_CH1_BandPass_Filter_CAL"];
channels2 = ["biceps_vasL_EMG_CH1_BandPass_Filter_CAL", "tibial_gasL_EMG_CH1_BandPass_Filter_CAL", "tibial_gasL_EMG_CH2_BandPass_Filter_CAL", "biceps_vasL_EMG_CH2_BandPass_Filter_CAL"];
timestamps = ["biceps_vasR_TimestampSync_Unix_CAL", "tibial_gasR_TimestampSync_Unix_CAL", "tibial_gasR_TimestampSync_Unix_CAL", "biceps_vasR_TimestampSync_Unix_CAL" ];
timestamps2 = ["biceps_vasL_TimestampSync_Unix_CAL", "tibial_gasL_TimestampSync_Unix_CAL", "tibial_gasL_TimestampSync_Unix_CAL", "biceps_vasL_TimestampSync_Unix_CAL" ];
mvc_maxs = [];
data = [];
n_subjects = 10;
Fs = 1024;
cutoff= 15;
[p,q]=butter(4,cutoff/(Fs/2),'low');
%n_subjects
for id = 1:30
    %Find MVC
    final_path = path + "Sujeto" + int2str(id) + "\MVC\" + prefix(id);
    for m_id = 1:length(muscles)
        Rmuscle = load(final_path+muscles(m_id)+"_Calibrated_PC.mat");
        Rmuscle_trials = abs(Rmuscle.(channels(m_id)));
        Rmuscle_trials= filtfilt(p,q, Rmuscle_trials);
        Rmuscle.timestamp= (Rmuscle.(timestamps(m_id)) - Rmuscle.(timestamps(m_id))(1))/1000;
        %figure(1); plot(Rmuscle.timestamp, Rmuscle_trials)
        c_fig= figure(1);
        [xi,~] = getpts(c_fig);
        %calculate RMS from three contractions of the MVC and save it in variable MAX
        MAX= [];
        for i=1:(length(xi)/2)-1
            MAX(i)= rms(Rmuscle_trials(find(Rmuscle.timestamp>xi((2*i)-1),1):find(Rmuscle.timestamp>xi(2*i),1)));
        end
        MAX= mean(MAX);
        mvc_maxs(id,m_id) = MAX;
    end
  %%  
    pre_path = path + "Sujeto" + int2str(id) + "\TEST\";
    
    for m_id = 1:length(folders)
        aux_path = pre_path + folders(m_id);     
        names = dir(aux_path);
        files = {names.name};
        file = char(files(4));
        num = file(length(char(prefix(id)))+8);
        
        final_path = aux_path + "\" +  prefix(id) + "Session" + num;
        %IMU data
        imu_pie5= load(final_path+"_imu_der_Calibrated_PC");
        %Time vector
        imup5.timestamp= (imu_pie5.imu_der_TimestampSync_Unix_CAL - imu_pie5.imu_der_TimestampSync_Unix_CAL(1))/1000;
        imup5.gyroY = imu_pie5.imu_der_Gyro_Y_CAL*-1;
        %Moving average filter
        windowSize = 30; 
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
        imup5.gyroY = filter(b,a,imup5.gyroY);
        %Find signal peaks 
        [val1,loc1_imu] = findpeaks(imup5.gyroY, 'MinPeakDistance',100,'MinPeakHeight',30);
        figure; plot(imup5.timestamp,imup5.gyroY);
        findpeaks(imup5.gyroY, 'MinPeakDistance',100,'MinPeakHeight',30);
        [val2,loc2] = findpeaks(-imup5.gyroY, 'MinPeakDistance',50,'MinPeakHeight',30);
        %hold on; plot(-imup5.gyroY);
        %Invert signal and find positive and negative peaks
        %findpeaks(-imup5.gyroY, 'MinPeakDistance',50,'MinPeakHeight',30);
        %zero-crossings in the gyro signal
        zx=  imup5.gyroY - (-imup5.gyroY);
        zeros1 = transpose(find(zx(1:end-1)>0 & zx(2:end) < 0));
        zeros2 = transpose(find(zx(1:end-1)<0 & zx(2:end) > 0));
        zerosn= [zeros1, zeros2];
        zerosn= sort(zerosn);  
        %Indices related to zero-crossings
        novo=[];
        for i=1:length(loc1_imu)-1
            if (length(find(zerosn<loc1_imu(i)))) == 0
               i=i+1; 
            end
            novo(end+1)= zerosn(find(zerosn>loc1_imu(i),1));
            novo(end+1)= zerosn(length(find(zerosn<loc1_imu(i))));
            
            
            %end
        end
        novo=sort(novo);
        y=[];
        for i=1:length(novo)
            y(i)= imup5.timestamp(novo(i));
        end
        hold on; plot(novo, y,'r*')
        %Extracting swing phase and stance phase times
        sw=[];
        st=[];
        for pt=1:(length(y)/2)
            sw(pt)= y(2*pt)- y(2*pt-1);
            if pt>1 
                st(pt-1)= y(2*pt-1)-y(2*pt-2);
            end
        end
        sw= mean(sw);
        st = mean(st);
        dist_imu=[];
        for i=1:(length(loc1_imu)-1)
            dist_imu(end+1)= loc1_imu(i+1)-loc1_imu(i);
        end
        %Resampling vectors of gait cycles to create single matrix and
        %extract features
        min_dis_imu= min(dist_imu);
        ciclos_imu= length(dist_imu);
        matriz_aj_imu= zeros(ciclos_imu, min_dis_imu);
        for i=1:size(matriz_aj_imu,1) %cantidad de datos en primera dimensión
            matriz_aj_imu(i,:)= resample(imup5.gyroY(loc1_imu(i):loc1_imu(i+1)-1), min_dis_imu, dist_imu(i)); %seleccionar señal de primer pico a siguiente 
        end
        %Plot mean signal with standard deviation of all gait cycles
        means_sm_imu= mean(matriz_aj_imu,1);
        std_sm_imu= std(matriz_aj_imu,1);    
        ejex_imu= linspace(0,100, min_dis_imu);
        a_imu= smoothdata(means_sm_imu+std_sm_imu, "gaussian", 10);
        b_imu= smoothdata(means_sm_imu-std_sm_imu, "gaussian", 10);
        
        %EMG signal for Biceps, Vastus, Tibialis and Gastrocnemius  
        emgBV = load(final_path+"_biceps_vasR_Calibrated_PC.mat");
        emgTG = load(final_path+"_tibial_gasR_Calibrated_PC.mat");
        mode_Rbiceps= emgBV.biceps_vasR_EMG_CH2_BandPass_Filter_CAL;
        mode_Rgastro= emgTG.tibial_gasR_EMG_CH2_BandPass_Filter_CAL;
        mode_Rtibial= emgTG.tibial_gasR_EMG_CH1_BandPass_Filter_CAL;
        mode_Rvasto= emgBV.biceps_vasR_EMG_CH1_BandPass_Filter_CAL;
        
        %Time vectors
        mode_timeBV = (emgBV.biceps_vasR_TimestampSync_Unix_CAL - emgBV.biceps_vasR_TimestampSync_Unix_CAL(1))/1000;
        mode_timeTG = (emgTG.tibial_gasR_TimestampSync_Unix_CAL - emgTG.tibial_gasR_TimestampSync_Unix_CAL(1))/1000;
        %Structure EMG signals
        muscle_signals = struct;
        muscle_signals.biceps = mode_Rbiceps;
        muscle_signals.gastro = mode_Rgastro;
        muscle_signals.tibial = mode_Rtibial;
        muscle_signals.vasto = mode_Rvasto;
        muscle_signals.biceps_t = mode_timeBV;
        muscle_signals.gastro_t = mode_timeTG;
        muscle_signals.tibial_t = mode_timeTG;
        muscle_signals.vasto_t = mode_timeBV;
        muscle_signals.names = ["biceps","gastro","tibial","vasto"];
        %Processing, filtering and normalization (envelope)
        for name_id = 1:length(muscle_signals.names)
            signal = muscle_signals.(muscle_signals.names(name_id));
            signal = filtfilt(p,q, signal);
            signal = (abs(signal)./mvc_maxs(id,m_id))*100;
            w=200;
            envelope15 = sqrt(movmean((signal.^2), w));
            %Iterating through the peaks detected in the IMU data (loc1_imu) and aligning them with the EMG signals.
            mode=[];
            for i=1:length(loc1_imu)-1
                mode(i) = find(muscle_signals.(muscle_signals.names(name_id)+"_t")>imup5.timestamp(loc1_imu(i)),1);
            end
            mode = unique(mode);
            %Determining distances between consecutive peaks in the IMU signal and using these to segment the EMG signal
            dist=[];
            for i=1:(length(mode)-1)
                ayuda = mode(i+1)- mode(i);
                if ayuda ~= 0
                    dist(end+1)= mode(i+1)- mode(i);
                end
            end
            %
            min_dist= min(dist);
            ciclos= length(dist);
            matriz_aj= zeros(ciclos, min_dist);
            for i=1:size(matriz_aj,1) %cantidad de datos en primera dimensión
                matriz_aj(i,:)= resample(envelope15(mode(i):mode(i+1)-1), min_dist, dist(i)); %seleccionar señal de primer pico a siguiente 
            end
            %Extracting features (RMS, mean, std) per cycle
            rms_mode= rms(matriz_aj,2);
            means_sm = mean(matriz_aj,1);
            std_sm = std(matriz_aj,1);
            a= smoothdata(means_sm+std_sm, "gaussian", 10);
            b= smoothdata(means_sm-std_sm, "gaussian", 10);
            ejex= linspace(0,100, min_dist);
            
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_mean_rms") = mean(rms_mode);
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_std_rms") = std(rms_mode);
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_ejex") = ejex;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_a") = a;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_b") = b;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_mean") = means_sm;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_std") = std_sm;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_imu_ejey") =means_sm_imu;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_imu_tiempo") =ejex_imu;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_a") = a_imu;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_b") = b_imu;
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_sw") = sw;   
            data.("s"+int2str(id)).(folders(m_id)+"_"+muscle_signals.names(name_id)+"_st") = st;
                
        end
    end 
end

close all;
clearvars -except data
