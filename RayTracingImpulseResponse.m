clear; close all; clc;

%% SetUp
room_dim = [6, 4, 3.4];
source_pos = [0.3, 0.3, 0.5];
mic_pos = [3.3, 2.3, 1.3];
mic_radius = 0.0875;
impResTime = 10;

plotRoom(room_dim,mic_pos,source_pos,1)

%% Generate Rays
N = 5000;
rng(0)
rays = RandSampleSphere(N);

%% Reflections and Scattering Coefficients
FVect = [125 250 500 1000 2000 4000 8000];

abs_coeffs = [];
abs_coeffs(:,1) = [0.02,0.02,0.03,0.03,0.04,0.05,0.05]; %Concrete
abs_coeffs(:,2) = [0.70,0.45,0.65,0.60,0.75,0.65,0.65]; %AbsPanels
abs_coeffs(:,3) = [0.14,0.10,0.06,0.08,0.10,0.10,0.10]; %Door

scatt_coeffs = [];
scatt_coeffs(:,1) = [0.30,0.50,0.60,0.60,0.70,0.70,0.70];
scatt_coeffs(:,2) = [0.20,0.40,0.50,0.50,0.60,0.60,0.60];
scatt_coeffs(:,3) =  [0.30,0.50,0.60,0.60,0.70,0.70,0.70];

num_FBands = length(FVect);

abs_map = {};
abs_map{1} = zeros(round(room_dim(2)*100),round(room_dim(3)*100),num_FBands);
abs_map{2} = zeros(round(room_dim(2)*100),round(room_dim(3)*100),num_FBands);
abs_map{3} = zeros(round(room_dim(1)*100),round(room_dim(3)*100),num_FBands);
abs_map{4} = zeros(round(room_dim(1)*100),round(room_dim(3)*100),num_FBands);
abs_map{5} = zeros(round(room_dim(1)*100),round(room_dim(2)*100),num_FBands);
abs_map{6} = zeros(round(room_dim(1)*100),round(room_dim(2)*100),num_FBands);

scatt_map = {};
scatt_map{1} = zeros(round(room_dim(2)*100),round(room_dim(3)*100),num_FBands);
scatt_map{2} = zeros(round(room_dim(2)*100),round(room_dim(3)*100),num_FBands);
scatt_map{3} = zeros(round(room_dim(1)*100),round(room_dim(3)*100),num_FBands);
scatt_map{4} = zeros(round(room_dim(1)*100),round(room_dim(3)*100),num_FBands);
scatt_map{5} = zeros(round(room_dim(1)*100),round(room_dim(2)*100),num_FBands);
scatt_map{6} = zeros(round(room_dim(1)*100),round(room_dim(2)*100),num_FBands);

%% Define wall materials
[abs_map, scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,1,abs_coeffs(:,1),...
    scatt_coeffs(:,2),[0,0],[room_dim(2),room_dim(3)]); %SmallWall
[abs_map, scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,2,abs_coeffs(:,1),...
    scatt_coeffs(:,2),[0,0],[room_dim(2),room_dim(3)]); %OpSmallWall
[abs_map, scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,3,abs_coeffs(:,1),...
    scatt_coeffs(:,2),[0,0],[room_dim(1),room_dim(3)]); %LargeWall
[abs_map, scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,4,abs_coeffs(:,1),...
    scatt_coeffs(:,2),[0,0],[room_dim(1),room_dim(3)]); %OpLargeWall
[abs_map, scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,5,abs_coeffs(:,1),...
    scatt_coeffs(:,2),[0,0],[room_dim(1),room_dim(2)]); %Floor
[abs_map, scatt_map] = updateAbsScattCoeffs(abs_map,scatt_map,6,abs_coeffs(:,1),...
    scatt_coeffs(:,2),[0,0],[room_dim(1),room_dim(2)]); %Ceiling

%% Reflections map
ref_map = cell(1,num_FBands);
for i = 1:numel(abs_map)
    ref_map{i} = sqrt(1-abs_map{i});
end

%% Energy Histogram
histTimeStep = 0.0010;
nTBins = round(impResTime/histTimeStep);
nFBins = length(FVect);
TFHist = zeros(nTBins,nFBins);

%% Ray Tracing
for iBand = 1:nFBins
    fprintf("Calculating rays for band %d\n",iBand)
    % Perform ray tracing independently for each frequency band.
    for iRay = 1:size(rays,1)
        % Select ray direction
        ray = rays(iRay,:);
        % All rays start at the source/transmitter
        ray_xyz = source_pos;
        % Set initial ray direction. This direction changes as the ray is
        % reflected off surfaces.
        ray_dxyz = ray;
        % Initialize ray travel time. Ray tracing is terminated when the
        % travel time exceeds the impulse response length.
        ray_time = 0;
        % Initialize the ray energy to a normalized value of 1.     Energy
        % decreases when the ray hits a surface.
        ray_energy = 1;

        while (ray_time <= impResTime)

            % Determine the surface that the ray encounters
            [surfaceofimpact,displacement] = getImpactWall(ray_xyz,...
                                             ray_dxyz,room_dim);
            
            % Determine the distance traveled by the ray
            distance = sqrt(sum(displacement.^2));

            % Determine the coordinates of the impact point
            impactCoord = ray_xyz+displacement;

            if surfaceofimpact > 4
                pointOfImpact = [impactCoord(1),impactCoord(2)];
            elseif surfaceofimpact < 3
                pointOfImpact = [impactCoord(2),impactCoord(3)];
            else
                pointOfImpact = [impactCoord(1),impactCoord(3)];
            end

            % Update ray location/source
            ray_xyz = impactCoord;

            % Update cumulative ray travel time
            c = 343; % speed of light (m/s)
            ray_time = ray_time+distance/c;

            % Apply surface reflection to ray's energy
            % This is the amount of energy that is not lost through
            % absorption.

            ReflectionAtPoint = ref_map{surfaceofimpact}(ceil(100*pointOfImpact(1)),ceil(100*pointOfImpact(2)),iBand);
            ray_energy = ray_energy*ReflectionAtPoint;

            % Apply diffuse reflection to ray energy
            % This is the fraction of energy used to determine what is
            % detected at the receiver
            DifussionAtPoint = scatt_map{surfaceofimpact}(ceil(100*pointOfImpact(1)),ceil(100*pointOfImpact(2)),iBand);
            rayrecv_energy = ray_energy*DifussionAtPoint;

            % Determine impact point-to-receiver direction.
            rayrecvvector = mic_pos-impactCoord;

            % Determine the ray's time of arrival at receiver.
            distance = sqrt(sum(rayrecvvector.*rayrecvvector));
            recv_timeofarrival = ray_time+distance/c;

            if recv_timeofarrival>impResTime
                break
            end

            if ray_energy < 0.000001
                break
            end

            % Determine amount of diffuse energy that reaches the receiver.
            % See (5.20) in [2].

            % Compute received energy
            N = getWallNormalVector(surfaceofimpact);
            cosTheta = sum(rayrecvvector.*N)/(sqrt(sum(rayrecvvector.^2)));
            cosAlpha = sqrt(sum(rayrecvvector.^2)-mic_radius^2)/sum(rayrecvvector.^2);
            E = (1-cosAlpha)*2*cosTheta*rayrecv_energy;

            % Update energy histogram
            tbin = floor(recv_timeofarrival/histTimeStep + 0.5);
            TFHist(tbin,iBand) = TFHist(tbin,iBand) + E;

            % Compute a new direction for the ray.
            % Pick a random direction that is in the hemisphere of the
            % normal to the impact surface.
            d = rand(1,3);
            d = d/norm(d);
            if sum(d.*N)<0
                d = -d;
            end

            % Derive the specular reflection with respect to the incident
            % wall
            ref = ray_dxyz-2*(sum(ray_dxyz.*N))*N;

            % Combine the specular and random components
            d = d/norm(d);
            ref = ref/norm(ref);
            ray_dxyz = DifussionAtPoint*d+(1-DifussionAtPoint)*ref;
            ray_dxyz = ray_dxyz/norm(ray_dxyz);
        end
    end
end

%% View
figure(1)
bar(histTimeStep*(0:size(TFHist)-1),TFHist)
grid on
xlabel("Time (s)")
legend(["125 Hz","250 Hz","500 Hz","1000 Hz","2000 Hz","4000 Hz"])

%% Generate Impulse Response
fs = 44100;
V = prod(room_dim);
t0 = ((2*V*log(2))/(4*pi*c^3))^(1/3); % eq 5.45 in [2]
poissonProcess = [];
timeValues = [];
t = t0;
while (t<impResTime)
    timeValues = [timeValues t]; %#ok
    % Determine polarity.
    if (round(t*fs)-t*fs) < 0 
        poissonProcess = [poissonProcess 1]; %#ok
    else
        poissonProcess = [poissonProcess -1];%#ok
    end
    % Determine the mean event occurence (eq 5.44 in [2])
    mu = min(1e4,4*pi*c^3*t^2/V); 
    % Determine the interval size (eq. 5.44 in [2])
    deltaTA = (1/mu)*log(1/rand); % eq. 5.43 in [2])
    t = t+deltaTA;
end
randSeq = zeros(ceil(impResTime*fs),1);
for index=1:length(timeValues)
    randSeq(round(timeValues(index)*fs)) = poissonProcess(index);
end
flow = [115 225 450 900 1800 3600 7200];
fhigh = [135 275 550 1100 2200 4400 8800];
NFFT = 8192;
win = hann(882,"symmetric");
sfft = dsp.STFT(Window = win,OverlapLength=441,FFTLength=NFFT,FrequencyRange="onesided");
isfft = dsp.ISTFT(Window=win,OverlapLength=441,FrequencyRange="onesided");
F = sfft.getFrequencyVector(fs);
RCF = zeros(length(FVect),length(F));
for index0 = 1:length(FVect)
    for index=1:length(F)
        f = F(index);
        if f<FVect(index0) && f>=flow(index0)
            RCF(index0,index) = .5*(1+cos(2*pi*f/FVect(index0)));
        end
        if f<fhigh(index0) && f>=FVect(index0)
            RCF(index0,index) = .5*(1-cos(2*pi*f/(FVect(index0)+1)));
        end
    end
end
frameLength = 441;
numFrames = length(randSeq)/frameLength;
y = zeros(length(randSeq),numel(FVect));
for index=1:numFrames
    x = randSeq((index-1)*frameLength+1:index*frameLength);
    X = sfft(x);    
    X = X.*RCF.';
    y((index-1)*frameLength+1:index*frameLength,:) = isfft(X);
end
impTimes = (1/fs)*(0:size(y,1)-1);
hisTimes = histTimeStep/2 + histTimeStep*(0:nTBins);
W = zeros(size(impTimes,2),numel(FVect));
BW = fhigh-flow;
for k=1:size(TFHist,1)
    gk0 = floor((k-1)*fs*histTimeStep)+1;
    gk1 = floor(k*fs*histTimeStep);
    yy = y(gk0:gk1,:).^2;
    val = sqrt(TFHist(k,:)./sum(yy,1)).*sqrt(BW/(fs/2));
    for iRay=gk0:gk1
        W(iRay,:)= val;
    end
end
%% Create Impulse Response
y_2 = y.*W;
ip = sum(y_2,2);
ip = ip./max(abs(ip));
ip = ip + rand(length(ip),1)/1000;
vectorTiempo = (1/fs)*(0:numel(ip)-1);
figure
plot(vectorTiempo,ip.^2)
grid on
xlabel("Time (s)")
ylabel("Impulse Response")

schroeder_cumsum = cumsum(flipud(ip.^2));
schroeder_normalized = schroeder_cumsum / max(schroeder_cumsum);
L = flipud(10*log10(schroeder_normalized));
[SignalSize, ~] = size(ip); 

vectorTiempo = (1:SignalSize)/fs;
[LF_EDT,LF_T20,LF_T30,LF_T60] = GetLinearFits(L,vectorTiempo);

[~, EDT_Idx] = min(abs(LF_EDT+60));
T60delEDT = vectorTiempo(EDT_Idx)
[~, T20_Idx] = min(abs(LF_T20+60));
T60delT20 = vectorTiempo(T20_Idx)
[~, T30_Idx] = min(abs(LF_T30+60));
T60delT30 = vectorTiempo(T30_Idx)
[~, T60_Idx] = min(abs(LF_T60+60));
T60delT60 = vectorTiempo(T60_Idx)

figure
plot(vectorTiempo,L,'LineWidth',2)
title('Schroeder Integration')
xlabel('Time (s)')
ylabel('Decay (Db)')
hold on
plot(vectorTiempo,LF_EDT)
plot(vectorTiempo,LF_T20)
plot(vectorTiempo,LF_T30)
plot(vectorTiempo,LF_T60)
legend('Energy Decay Cruve','Linear fit for EDT','Linear fit for T20','Linear fit for T30', ...
     'Linear fit for T60','Location','southwest')
title('Linear fits')
xlabel('Time (s)')
ylabel('Decay (Db)')
ylim([-80 0])
hold off





