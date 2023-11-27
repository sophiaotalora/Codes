clear all; close all; clc;
path = "C:\Users\NTA\MATLAB Drive\ICORR\";
%path= "C:\Users\sophi\MATLAB Drive\ICORR\";
windowSize = 30;
b = (1/windowSize)*ones(1,windowSize);
a = 1;ex
fs_imu= 52;
for id = 1:10
    %% Borg scale
    %Segement data in three states of fatigue
    f_path = path + "S" + int2str(id);
    cd(f_path);  
    Borg= load('BORG.txt');
    escala_borg_10s = Borg(:,2)';
    t2= Borg(:,1)';
    esc1= []; temp1= [];
    esc2= []; temp2= [];
    esc3= []; temp3= [];
    for i=1:length(escala_borg_10s)       
        if escala_borg_10s(i)<= 3 
            esc1= [esc1, round(escala_borg_10s(i))];
            temp1= [temp1, t2(i)];
        end
        if escala_borg_10s(i)> 3 && escala_borg_10s(i)<= 7  
            esc2= [esc2, round(escala_borg_10s(i))];
            temp2= [temp2, t2(i)];
        end
        if escala_borg_10s(i)> 7 && escala_borg_10s(i)<= 11 
            esc3= [esc3, round(escala_borg_10s(i))];
            temp3= [temp3, t2(i)];
        end
    end
    if esc1(:, end) ~= 3
        esc1= [esc1, round((esc1(:, end)+esc2(1))/2)];
        temp1= [temp1, ((temp1(:, end)+temp2(1))/2)];
    end
    
    if esc2(:, end) ~= 7
        esc2= [esc2, round((esc2(:, end)+esc3(1))/2)];
        temp2= [temp2, ((temp2(:, end)+temp3(1))/2)];
    end
    temp2(1)= temp1(:, end)+1;
    temp3(1)= temp2(:, end)+1;
    
    figure;
    plot(t2, escala_borg_10s,'-x')
    xlabel('Tiempo (segundos)');
    ylabel('Escala de Borg');
    legend('Valores Originales (20 segundos)');
    
    %% Optical fiber data (elbow angle sensor)
    %
    ANG= importdata('EYF.txt');
    ang = ANG(:,2);
    Fs= 52;
    ang=lowpass(ang,0.5,Fs);
    %figure; plot(ang);
    %Find peaks to divide the signal into elbow flexion-extension cycles
    [r2, d2]= findpeaks(ang, 'MinPeakDistance', 60);
    len_w2=[];
    for s=1:length(d2)-1
            len_w2= [len_w2, (d2(s)-d2(s+1))*(-1)];
    end
    min_gen= min(len_w2);
    FIB.time= (0:length(ang)-1)*(1/(length(ang)/temp3(:, end)));
    %figure; plot(ang);
    %Segmentation based on Borg Scale Ranges
    FIB.Borg1_3 = ang(1:find(FIB.time>temp1(:, end),1));
    FIB.Borg4_7= ang(find(FIB.time>temp2(1),1):find(FIB.time>temp2(:, end),1));
    FIB.Borg8_10= ang(find(FIB.time>temp3(1),1):find(FIB.time>(temp3(:, end)-1),1)); 
    
    FIB.t_IMU1_3=FIB.time(1:find(FIB.time>temp1(:, end),1));
    FIB.t_IMU4_7= FIB.time(find(FIB.time>temp2(1),1):find(FIB.time>temp2(:, end),1));
    FIB.t_IMU8_10= FIB.time(find(FIB.time>temp3(1),1):find(FIB.time>(temp3(:, end)-1),1)); 
    st=[1480,1578,1229,1475,1316,1206,1531,1188,0,1449];
    ed=[1875,2059,1803,1951,1695,1667,2073,1936,686,2442];
    FIB.Borg1_3 = FIB.Borg1_3([1:st(id), ed(id):end]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1);yyaxis left; plot((FIB.Borg1_3)-40); hold on; findpeaks((FIB.Borg1_3)-40, 'MinPeakDistance', 60); hold on;%subplot(3,1,2);plot(FIB.Borg4_7);subplot(3,1,3);plot(FIB.Borg8_10);
    figure(2);yyaxis left; plot((FIB.Borg4_7)-40); hold on; findpeaks((FIB.Borg4_7)-40, 'MinPeakDistance', 60); hold on;
    figure(3);yyaxis left; plot((FIB.Borg8_10)-40); hold on; findpeaks((FIB.Borg8_10)-40, 'MinPeakDistance', 60); hold on;
    borg_sep=["Borg1_3", "Borg4_7", "Borg8_10"];
    [r3, d3]= findpeaks(FIB.(borg_sep(1)), 'MinPeakDistance', 60);
    [r7, d7]= findpeaks(FIB.(borg_sep(2)), 'MinPeakDistance', 60);
    [r10, d10]= findpeaks(FIB.(borg_sep(3)), 'MinPeakDistance', 60);
    %Feature extraction
    for z=1:length(borg_sep)
        sig= FIB.(borg_sep(z));
        [r, d]= findpeaks(sig, 'MinPeakDistance', 60);
        %figure;plot(sig); hold on; plot(d,r,'x');  title("S"+(int2str(id))+"borg"+(int2str(z)));  
        fib_norm=[];fib_mean=[];fib_rms=[];fib_ROM=[];fib_std=[];
        for s=1:length(d)-1
            len_w= d(s)-d(s+1);
            fib_norm= [fib_norm, (len_w*100)/min_gen*(-1)];
            fib_mean= [fib_mean, mean(sig(d(s):d(s+1)))];
            fib_std= [fib_std, std(sig(d(s):d(s+1)))];
            fib_rms= [fib_rms, rms(sig(d(s):d(s+1)))];
            fib_ROM= [fib_ROM, max(sig(d(s):d(s+1))) - min(sig(d(s):d(s+1)))];
        end
        FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_norm= fib_norm'; FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_mean= fib_mean';
        FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_rms= fib_rms'; FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_std= fib_std';
        FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_ROM= fib_ROM';
        FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)+"fib_FIN")=[FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_norm,FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_std, FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_mean, FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_rms, FORMA2.("S"+int2str(id)).FIB.(borg_sep(z)).fib_ROM];
        %
        w = 3;    
        w_samples = w * Fs;
        num_w = floor(length(sig)/w_samples);  
        %struct1(z)
        ws = zeros(w_samples, num_w);   
        for k = 1:num_w
            inicio = (k - 1) * w_samples + 1;
            fin = inicio + w_samples - 1;
            ws(:, k) = sig(inicio:fin);    
        end
        FIB.(borg_sep(z)+"ws")= ws;  
    end  
    FORMA2.("S"+int2str(id)).FIB.fib_FIN= [FORMA2.("S"+int2str(id)).FIB.(borg_sep(1)+"fib_FIN");FORMA2.("S"+int2str(id)).FIB.(borg_sep(2)+"fib_FIN");FORMA2.("S"+int2str(id)).FIB.(borg_sep(3)+"fib_FIN")];
    Az= [size(FORMA2.("S"+int2str(id)).FIB.(borg_sep(1)+"fib_FIN")); size(FORMA2.("S"+int2str(id)).FIB.(borg_sep(2)+"fib_FIN")); size(FORMA2.("S"+int2str(id)).FIB.(borg_sep(3)+"fib_FIN"))];
    a=repelem(1,Az(1)); b= repelem(2,Az(2)); c= repelem(3,Az(3));
    FORMA2.("S"+int2str(id)).FIB.fib_teste=[a'; b';c'];

    %% IMU
    imu1=load("IMU1.txt");
    imu2= load('IMU2.txt');
    % Extract gyro and acceleration from IMU.
    IMU.IMU1.gyroX = imu1(:,7); IMU.IMU1.gyroY = imu1(:,8); IMU.IMU1.gyroZ = imu1(:,9);
    IMU.IMU1.accX = imu1(:,1); IMU.IMU1.accY = imu1(:,2); IMU.IMU1.accZ = imu1(:,3);
    IMU.IMU2.gyroX = imu2(:,7); IMU.IMU2.gyroY = imu2(:,8); IMU.IMU2.gyroZ = imu2(:,9);
    IMU.IMU2.accX = imu2(:,1); IMU.IMU2.accY = imu2(:,2); IMU.IMU2.accZ = imu2(:,3);
  
    %Filter
    filt=["gyroX", "gyroY", "gyroZ", "accX","accY", "accZ"];
    imux= ["IMU1","IMU2"];
    for i=1:length(imux)
        for j=1:length(filt)
            IMU.(imux(i)).(filt(j))= filter(b,a,IMU.(imux(i)).(filt(j)));
            IMU.(imux(i)).time= (0:length(IMU.IMU1.gyroX)-1)*(1/(length(IMU.IMU1.gyroX)/temp3(:, end)));
        end
    end
    figure;
    for j=1:length(filt)
        subplot(3,3,j); plot(IMU.(imux(1)).(filt(j)));
    end
    sgtitle("S" + int2str(id)+ " IMU1")
    %
    figure;
    for j=1:length(filt)
        subplot(3,3,j); plot(IMU.(imux(2)).(filt(j)));
    end
    sgtitle("S" + int2str(id) + " IMU2")  
    %Segments data into subranges of Borg scale
    for i=1:length(imux)
        for j=1:length(filt)
            IMU.("IMU"+int2str(i)).("Borg1_3"+filt(j)) = IMU.("IMU"+int2str(i)).(filt(j))(1:find(IMU.("IMU"+int2str(i)).time>temp1(:, end),1));
            IMU.("IMU"+int2str(i)).("Borg4_7"+filt(j))= IMU.("IMU"+int2str(i)).(filt(j))(find(IMU.("IMU"+int2str(i)).time>temp2(1),1):find(IMU.("IMU"+int2str(i)).time>temp2(:, end),1));
            IMU.("IMU"+int2str(i)).("Borg8_10"+filt(j))= IMU.("IMU"+int2str(i)).(filt(j))(find(IMU.("IMU"+int2str(i)).time>temp3(1),1):find(IMU.("IMU"+int2str(i)).time>(temp3(:, end)-1),1));           
        
        end 
    end
    %
    for i=1:length(imux)
        IMU.("IMU"+int2str(i)).t_IMU1_3=IMU.("IMU"+int2str(i)).time(1:find(IMU.("IMU"+int2str(i)).time>temp1(:, end),1));
        IMU.("IMU"+int2str(i)).t_IMU4_7= IMU.("IMU"+int2str(i)).time(find(IMU.("IMU"+int2str(i)).time>temp2(1),1):find(IMU.("IMU"+int2str(i)).time>temp2(:, end),1));
        IMU.("IMU"+int2str(i)).t_IMU8_10= IMU.("IMU"+int2str(i)).time(find(IMU.("IMU"+int2str(i)).time>temp3(1),1):find(IMU.("IMU"+int2str(i)).time>(temp3(:, end)-1),1));
    end
    
   figure; subplot(3,1,1);plot(IMU.("IMU"+int2str(2)).t_IMU1_3,IMU.("IMU"+int2str(2)).("Borg1_3"+filt(3)))
   subplot(3,1,2);plot(IMU.("IMU"+int2str(2)).t_IMU4_7,IMU.("IMU"+int2str(2)).("Borg4_7"+filt(3)))
   subplot(3,1,3);plot(IMU.("IMU"+int2str(2)).t_IMU8_10,IMU.("IMU"+int2str(2)).("Borg8_10"+filt(3)))
   
    st=[1480,1578,1229,1475,1316,1206,1531,1188,0,1449];
    ed=[1875,2059,1803,1951,1695,1667,2073,1936,686,2442];
    for i=1:length(imux)
        for j=1:length(filt)
            IMU.("IMU"+int2str(i)).("Borg1_3"+filt(j)) = IMU.("IMU"+int2str(i)).("Borg1_3"+filt(j))([1:st(id), ed(id):end]);
        end 
    end
    
%%
    borg_sep=["Borg1_3", "Borg4_7", "Borg8_10"];

    for i=1:length(imux)
        for z=1:length(borg_sep)
            for j=1:length(filt)
                sig= IMU.("IMU"+int2str(i)).(borg_sep(z)+filt(j)); 
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean")=[]; FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std")=[];
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms")=[]; FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")=[];
                if z==1
                    for s=1:length(d3)-1
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean"),mean(sig(d3(s):d3(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std"),std(sig(d3(s):d3(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms"),rms(sig(d3(s):d3(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM"), max(sig(d3(s):d3(s+1))) - min(sig(d3(s):d3(s+1)))];
                    end
                elseif z==2
                    for s=1:length(d7)-1
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean"),mean(sig(d7(s):d7(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std"),std(sig(d7(s):d7(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms"),rms(sig(d7(s):d7(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM"), max(sig(d7(s):d7(s+1))) - min(sig(d7(s):d7(s+1)))];                    
                    end
                elseif z==3
                    for s=1:length(d10)-1
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean"),mean(sig(d10(s):d10(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std"),std(sig(d10(s):d10(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms"),rms(sig(d10(s):d10(s+1)))];
                        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM"), max(sig(d10(s):d10(s+1))) - min(sig(d10(s):d10(s+1)))];
                    end
                end
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean")= FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean")';
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std")= FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std")';
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms")= FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms")';
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")=FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")';
                FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+"imu_FIN")(:,(j*4)-3:j*4)=[FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_mean"),FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_std"),FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_rms"),FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"imu_ROM")];
                
                w = 3;    
                w_samples = w * fs_imu;
                num_w = floor(length(sig)/w_samples);   
                ws = zeros(w_samples, num_w);   
                for k = 1:num_w
                    inicio = (k - 1) * w_samples + 1;
                    fin = inicio + w_samples - 1;
                    ws(:, k) = sig(inicio:fin);    
                end
                %
                IMU.("IMU"+int2str(i)).(borg_sep(z)+filt(j)+"ws")= ws;  
                          
            end
                
        end 
        
        %Organize processed data into structure     
        Ai= [size(FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(1)+"imu_FIN")); size(FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(2)+"imu_FIN")); size(FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(3)+"imu_FIN"))];
        a=repelem(1,Ai(1)); b= repelem(2,Ai(2)); c= repelem(3,Ai(3));
        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).IMU_teste=[a'; b';c'];
        FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).imu_FIN= [FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(1)+"imu_FIN"); 
                                                              FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(2)+"imu_FIN");
                                                              FORMA2.("S"+int2str(id)).("IMU"+int2str(i)).(borg_sep(3)+"imu_FIN")];
    end
                %

    %}
    close all      
    FORMA2.("S"+int2str(id)).final= [FORMA2.("S"+int2str(id)).("IMU"+int2str(1)).imu_FIN,FORMA2.("S"+int2str(id)).("IMU"+int2str(2)).imu_FIN FORMA2.("S"+int2str(id)).FIB.fib_FIN, FORMA2.("S"+int2str(id)).FIB.fib_teste];
    FORMA2.("S"+int2str(id)).finalBORG= [FORMA2.("S"+int2str(id)).FIB.fib_teste];
    
    %}
end
clearvars -except FORMA2

