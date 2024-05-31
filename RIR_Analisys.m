function [RIRsq,DC,t,AcousticParams] = RIR_Analisys(method, input, room_dimensions,printFlag)
%RIR_Analisys Calculates acoustic parameters as in ISO3382. use method 1
%for wav file (input is name of wavfile), use method 2 for direct input,
%(structure like with .y and .fs)
    if method == 1
        [y,fs] = audioread(input);
    else
        y = input.y;
        fs = input.fs;
    end

    if size(y,2) > 1
        y = mean(y,2);
    end
    sound_duration = length(y)/fs;
    t = ((0:length(y)-1)/fs)';
    y = bandpass(y,[88.39 176.8],fs);
    RIRsq = y.^2;
    schroeder_cumsum = cumsum(flipud(RIRsq));
    schroeder_normalized = schroeder_cumsum / max(schroeder_cumsum);
    DC = flipud(10*log10(schroeder_normalized));
    
    [LF_EDT,LF_T20,LF_T30,LF_T60] = GetLinearFits(DC,t);
    
    [~, EDT_Idx1] = min(abs(LF_EDT+10));
    [~, EDT_Idx2] = min(abs(LF_EDT));
    T60delEDT = 6*(t(EDT_Idx1)-t(EDT_Idx2));
    
    [~, EDT_Idx1] = min(abs(LF_T20+20));
    [~, EDT_Idx2] = min(abs(LF_T20));
    T60delT20 = 3*(t(EDT_Idx1)-t(EDT_Idx2));
    
    [~, EDT_Idx1] = min(abs(LF_T30+30));
    [~, EDT_Idx2] = min(abs(LF_T30));
    T60delT30 = 2*(t(EDT_Idx1)-t(EDT_Idx2));
     
    [~,maxIdx] = max(RIRsq);
    RIRsq_trunc = RIRsq(maxIdx:end);
    Energy0_50 = trapz(RIRsq_trunc(1:round(0.05*fs)));
    Energy0_80 = trapz(RIRsq_trunc(1:round(0.08*fs)));
    Energy0_end = trapz(RIRsq_trunc(1:end));
    
    D50 = Energy0_50/Energy0_end;
    D80 = Energy0_80/Energy0_end;
    C50 = 10*log(D50/(1-D50));
    C80 = 10*log(D80/(1-D80));
    G = 10*log10(T60delT20/prod(room_dimensions))+45;
    
    AcousticParams = struct('T60delEDT',T60delEDT,'T60delT20',T60delT20,'T60delT30',T60delT30 ...
        ,'D50',D50,'D80',D80,'C50',C50,'C80',C80,'G',G);


    if printFlag == true
        figure
        plot(t,y)
        title('Impulse Response')
        xlabel('Time (s)')
        ylabel('Amplitude')
        xlim([0 sound_duration])

        figure
        plot(t,DC,"LineWidth",2)
        title('Schroeder Integration')
        xlabel('Time (s)')
        ylabel('Decay (Db)')
        ylim([-80 0.1])
        xlim([0 sound_duration])

        hold on
        plot(t,LF_EDT)
        plot(t,LF_T20)
        plot(t,LF_T30)
        legend('Decay Curve','Linear fit for EDT','Linear fit for T20','Linear fit for T30','Location','southwest')
        hold off

        fprintf("\nT60 from EDT: %f",T60delEDT)
        fprintf("\nT60 from T20: %f",T60delT20)
        fprintf("\nT60 from T30: %f",T60delT30)

        fprintf("\nD50: %f",D50)
        fprintf("\nC50: %f",C50)
        fprintf("\nC80: %f",C80)
        fprintf("\nG: %f",G)
        fprintf("\n")
    end 
end

