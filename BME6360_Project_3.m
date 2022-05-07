%BME6360 Project 3
%Jacey Davies

clear;
close all;
clc;

%Loading the data for subject 1
[signal, states, parameters]=load_data; %Then select the 9 .dat files for S001

%Collecting 800ms of signal samples after the start of each intensification
a = 1;
b = 1;
c = 1; %These are for tracking indices
samples = zeros(193,414720); %To collect all of the 800ms samples
samples_stand = zeros(193,345600); %For when the desired character is not intensified
samples_odd = zeros(193,69120); %For when the desired character is intensified
stim_type = zeros(1,414720); %To keep track of the stimulation type in one vector

for i = 1:64 %To go through all channels
    for j = 1:315954 %To go through all data points
    if (states.Flashing(j)==1) && (states.Flashing(j-1)==0) %To find the transition from 0 to 1
        samples(:,a) = signal(j:j+192,i); %Record the channel's samples in a column
        if (states.StimulusType(j) == 0)
            samples_stand(:,b) = signal(j:j+192,i); %For no desired character
            %I use j+192 to include the start data point + 800ms (or 192 data
            %points)
            
            stim_type(a) = 0; %Record the stimulation type
            b = b + 1; %To keep track of samples_stand index
        else
            samples_odd(:,c) = signal(j:j+192,i); %For desired character
            
            stim_type(a) = 1; %Record the stimulation type
            c = c + 1; %To keep track of samples_odd index
        end
        a = a + 1; %To keep track of samples index
    else
        %Skip non-transition points
    end
    end
end

%I know this complicated for loop is not the most efficient way to do this
%but it works (I think) and I'm too scared to change it at this point

%Finding the evoked potentials/means
mean_stand = mean(samples_stand,2); %Find mean of each row
mean_odd = mean(samples_odd,2); %Find mean of each row

time = 0:(1/240)*1000:800; %Time from transition point (0ms) to 800ms after

figure(1)
plot(time,mean_stand','r')
hold on
plot(time,mean_odd','k')
xlim([0 800])
xlabel('Time After Stimulus (ms)')
ylabel('Signal Amplitude (mV)') %Units??? Why is it so small??
hold off

%Computing the correlation between the state.StimulusType and the response 
%amplitude for each time sample and channel
correl = zeros(1,size(samples,1)); %I actually know the size of this so I preset the vector

for i = 1:size(samples,1) %To go through each xi or time sample
    m = corrcoef(stim_type,samples(i,:)'); %Stim type is a vector of 0 and 1 for stimulation type
    %samples(i,:)' is the "column" for each xi/time sample (see the explanation for correlation document that I based this on)
    
    correl(i) = m(2); %Corrcoef() makes a matrix where the second element is the actual correlation
end

figure(2)
plot(time,correl,'k')
xlim([0 800])
xlabel('Time After Stimulus (ms)')
ylabel('r for Standard vs. Oddball')

figure(3)
plot(time,correl.^2,'k')
xlim([0 800])
xlabel('Time After Stimulus (ms)')
ylabel('r^2 for Standard vs. Oddball')

%We redoing this because I'm so tired and I need that damn topograph
sample_time = .8; %s
sample_rate = 240; % Hz
num_channels = 64;
num_samples = sample_rate * sample_time;

% For all channels, collect a 800ms of signal samples after the start of each
% intensification, i.e., whenever state.Flashing changes from 0 to 1 (note: each
% character epoch of the data set starts at the first flash, i.e. state.Flashing=1 for the
% first data sample in each epoch).

flashing = states.Flashing(2:end);
flashing2 = states.Flashing(1:end-1);
changes =[0; flashing - flashing2];

index = find(changes); %location of sdifferent intensifications 

data = zeros(num_samples, num_channels, length(index));
        
    for j = 1:length(index)
data (:,:,j) = signal(index(j):index(j) + num_samples-1, :);
    end
    
% Compute the correlation between the state.StimulusType and the response
% amplitude for each time sample and channel.
data2 = permute(data,[3 1 2]); 

correlation = zeros(num_channels, num_samples);

for i = 1:num_channels
    for j = 1:num_samples
        correlation(i, j) = corr(squeeze(data2(:, j, i)), double(states.StimulusType(index)));
    end
end

datavector_100 = correlation(:,24).^2;
datavector_200 = correlation(:,48).^2;
datavector_300 = correlation(:,72).^2;
datavector_400 = correlation(:,96).^2;
datavector_500 = correlation(:,120).^2;

figure(4)
subplot(1,5,1)
topoplot(datavector_100,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,2)
topoplot(datavector_200,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,3)
topoplot(datavector_300,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,4)
topoplot(datavector_400,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,5)
topoplot(datavector_500,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

%Plotting average temporal response waveforms & temporal correlations for
%channels 11, 51, 56, & 60
samples_11 = zeros(193,6480);
samples_11_stand = zeros(193,5400);
samples_11_odd = zeros(193,1080);

samples_51 = zeros(193,6480);
samples_51_stand = zeros(193,5400);
samples_51_odd = zeros(193,1080);

samples_56 = zeros(193,6480);
samples_56_stand = zeros(193,5400);
samples_56_odd = zeros(193,1080);

samples_60 = zeros(193,6480);
samples_60_stand = zeros(193,5400);
samples_60_odd = zeros(193,1080);

b = 1;
c = 1;

    for i = 1:length(index)
samples_11(:,i) = signal(index(i):index(i)+192,11);
samples_51(:,i) = signal(index(i):index(i)+192,51);
samples_56(:,i) = signal(index(i):index(i)+192,56);
samples_60(:,i) = signal(index(i):index(i)+192,60);
        if (states.StimulusType(index(i)) == 0)
            samples_11_stand(:,b) = signal(index(i):index(i)+192,11); %For no desired character
            samples_51_stand(:,b) = signal(index(i):index(i)+192,51);
            samples_56_stand(:,b) = signal(index(i):index(i)+192,56);
            samples_60_stand(:,b) = signal(index(i):index(i)+192,60);
            
            b = b + 1; %To keep track of samples_stand index
        else
            samples_11_odd(:,c) = signal(index(i):index(i)+192,11); %For no desired character
            samples_51_odd(:,c) = signal(index(i):index(i)+192,51);
            samples_56_odd(:,c) = signal(index(i):index(i)+192,56);
            samples_60_odd(:,c) = signal(index(i):index(i)+192,60);
            
            c = c + 1; %To keep track of samples_odd index
        end
    end
    
correl_11 = zeros(1,193); %I actually know the size of this so I preset the vector
correl_51 = zeros(1,193);
correl_56 = zeros(1,193);
correl_60 = zeros(1,193);

for i = 1:193 %To go through each xi or time sample
    m = corrcoef(double(states.StimulusType(index)),samples_11(i,:)');
    correl_11(i) = m(2); %Corrcoef() makes a matrix where the second element is the actual correlation
    
     m = corrcoef(double(states.StimulusType(index)),samples_51(i,:)');
    correl_51(i) = m(2);
    
     m = corrcoef(double(states.StimulusType(index)),samples_56(i,:)');
    correl_56(i) = m(2);
    
     m = corrcoef(double(states.StimulusType(index)),samples_60(i,:)');
    correl_60(i) = m(2);
    
end

mean_11_stand = mean(samples_11_stand,2); %Find mean of each row
mean_11_odd = mean(samples_11_odd,2); %Find mean of each row

mean_51_stand = mean(samples_51_stand,2); %Find mean of each row
mean_51_odd = mean(samples_51_odd,2); %Find mean of each row

mean_56_stand = mean(samples_56_stand,2); %Find mean of each row
mean_56_odd = mean(samples_56_odd,2); %Find mean of each row

mean_60_stand = mean(samples_60_stand,2); %Find mean of each row
mean_60_odd = mean(samples_60_odd,2); %Find mean of each row

figure(5)
subplot(1,4,1)
plot(time,mean_11_stand','r')
hold on
plot(time,mean_11_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
ylabel('Signal Amplitude (mV)')
hold off

subplot(1,4,2)
plot(time,mean_51_stand','r')
hold on
plot(time,mean_51_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
hold off

subplot(1,4,3)
plot(time,mean_56_stand','r')
hold on
plot(time,mean_56_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
hold off

subplot(1,4,4)
plot(time,mean_60_stand','r')
hold on
plot(time,mean_60_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
hold off

han=axes(figure(5),'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Time After Stimulus (ms)');

figure(6)
subplot(1,4,1)
plot(time,correl_11.^2,'k')
xlim([0 800])
ylim([0 0.03])
ylabel('r^2')

subplot(1,4,2)
plot(time,correl_51.^2,'k')
xlim([0 800])
ylim([0 0.03])

subplot(1,4,3)
plot(time,correl_56.^2,'k')
xlim([0 800])
ylim([0 0.03])

subplot(1,4,4)
plot(time,correl_60.^2,'k')
xlim([0 800])
ylim([0 0.03])

han=axes(figure(6),'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Time After Stimulus (ms)');



%Loading the data for subject 2
[signal, states, parameters]=load_data; %Then select the 9 .dat files for S002

%Collecting 800ms of signal samples after the start of each intensification
a = 1;
b = 1;
c = 1; %These are for tracking indices
samples = zeros(193,414720); %To collect all of the 800ms samples
samples_stand = zeros(193,345600); %For when the desired character is not intensified
samples_odd = zeros(193,69120); %For when the desired character is intensified
stim_type = zeros(1,414720); %To keep track of the stimulation type in one vector

for i = 1:64 %To go through all channels
    for j = 1:315954 %To go through all data points
    if (states.Flashing(j)==1) && (states.Flashing(j-1)==0) %To find the transition from 0 to 1
        samples(:,a) = signal(j:j+192,i); %Record the channel's samples in a column
        if (states.StimulusType(j) == 0)
            samples_stand(:,b) = signal(j:j+192,i); %For no desired character
            %I use j+192 to include the start data point + 800ms (or 192 data
            %points)
            
            stim_type(a) = 0; %Record the stimulation type
            b = b + 1; %To keep track of samples_stand index
        else
            samples_odd(:,c) = signal(j:j+192,i); %For desired character
            
            stim_type(a) = 1; %Record the stimulation type
            c = c + 1; %To keep track of samples_odd index
        end
        a = a + 1; %To keep track of samples index
    else
        %Skip non-transition points
    end
    end
end

%I know this complicated for loop is not the most efficient way to do this
%but it works (I think) and I'm too scared to change it at this point

%Finding the evoked potentials/means
mean_stand = mean(samples_stand,2); %Find mean of each row
mean_odd = mean(samples_odd,2); %Find mean of each row

figure(7)
plot(time,mean_stand','r')
hold on
plot(time,mean_odd','k')
xlim([0 800])
xlabel('Time After Stimulus (ms)')
ylabel('Signal Amplitude (mV)') %Units??? Why is it so small??
hold off

%Computing the correlation between the state.StimulusType and the response 
%amplitude for each time sample and channel
correl = zeros(1,size(samples,1)); %I actually know the size of this so I preset the vector

for i = 1:size(samples,1) %To go through each xi or time sample
    m = corrcoef(stim_type,samples(i,:)'); %Stim type is a vector of 0 and 1 for stimulation type
    %samples(i,:)' is the "column" for each xi/time sample (see the explanation for correlation document that I based this on)
    
    correl(i) = m(2); %Corrcoef() makes a matrix where the second element is the actual correlation
end

figure(8)
plot(time,correl,'k')
xlim([0 800])
xlabel('Time After Stimulus (ms)')
ylabel('r for Standard vs. Oddball')

figure(9)
plot(time,correl.^2,'k')
xlim([0 800])
xlabel('Time After Stimulus (ms)')
ylabel('r^2 for Standard vs. Oddball')

%We redoing this because I'm so tired and I need that damn topograph
sample_time = .8; %s
sample_rate = 240; % Hz
num_channels = 64;
num_samples = sample_rate * sample_time;

% For all channels, collect a 800ms of signal samples after the start of each
% intensification, i.e., whenever state.Flashing changes from 0 to 1 (note: each
% character epoch of the data set starts at the first flash, i.e. state.Flashing=1 for the
% first data sample in each epoch).

flashing = states.Flashing(2:end);
flashing2 = states.Flashing(1:end-1);
changes =[0; flashing - flashing2];

index = find(changes); %location of sdifferent intensifications 

data = zeros(num_samples, num_channels, length(index));
        
    for j = 1:length(index)
data (:,:,j) = signal(index(j):index(j) + num_samples-1, :);
    end
    
% Compute the correlation between the state.StimulusType and the response
% amplitude for each time sample and channel.
data2 = permute(data,[3 1 2]); 

correlation = zeros(num_channels, num_samples);

for i = 1:num_channels
    for j = 1:num_samples
        correlation(i, j) = corr(squeeze(data2(:, j, i)), double(states.StimulusType(index)));
    end
end

datavector_100 = correlation(:,24).^2;
datavector_200 = correlation(:,48).^2;
datavector_300 = correlation(:,72).^2;
datavector_400 = correlation(:,96).^2;
datavector_500 = correlation(:,120).^2;

figure(10)
subplot(1,5,1)
topoplot(datavector_100,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,2)
topoplot(datavector_200,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,3)
topoplot(datavector_300,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,4)
topoplot(datavector_400,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

subplot(1,5,5)
topoplot(datavector_500,'eloc64.txt','EEG');
caxis([0 0.02])
colorbar

%Plotting average temporal response waveforms & temporal correlations for
%channels 11, 51, 56, & 60
samples_11 = zeros(193,6480);
samples_11_stand = zeros(193,5400);
samples_11_odd = zeros(193,1080);

samples_51 = zeros(193,6480);
samples_51_stand = zeros(193,5400);
samples_51_odd = zeros(193,1080);

samples_56 = zeros(193,6480);
samples_56_stand = zeros(193,5400);
samples_56_odd = zeros(193,1080);

samples_60 = zeros(193,6480);
samples_60_stand = zeros(193,5400);
samples_60_odd = zeros(193,1080);

b = 1;
c = 1;

    for i = 1:length(index)
samples_11(:,i) = signal(index(i):index(i)+192,11);
samples_51(:,i) = signal(index(i):index(i)+192,51);
samples_56(:,i) = signal(index(i):index(i)+192,56);
samples_60(:,i) = signal(index(i):index(i)+192,60);
        if (states.StimulusType(index(i)) == 0)
            samples_11_stand(:,b) = signal(index(i):index(i)+192,11); %For no desired character
            samples_51_stand(:,b) = signal(index(i):index(i)+192,51);
            samples_56_stand(:,b) = signal(index(i):index(i)+192,56);
            samples_60_stand(:,b) = signal(index(i):index(i)+192,60);
            
            b = b + 1; %To keep track of samples_stand index
        else
            samples_11_odd(:,c) = signal(index(i):index(i)+192,11); %For no desired character
            samples_51_odd(:,c) = signal(index(i):index(i)+192,51);
            samples_56_odd(:,c) = signal(index(i):index(i)+192,56);
            samples_60_odd(:,c) = signal(index(i):index(i)+192,60);
            
            c = c + 1; %To keep track of samples_odd index
        end
    end
    
correl_11 = zeros(1,193); %I actually know the size of this so I preset the vector
correl_51 = zeros(1,193);
correl_56 = zeros(1,193);
correl_60 = zeros(1,193);

for i = 1:193 %To go through each xi or time sample
    m = corrcoef(double(states.StimulusType(index)),samples_11(i,:)');
    correl_11(i) = m(2); %Corrcoef() makes a matrix where the second element is the actual correlation
    
     m = corrcoef(double(states.StimulusType(index)),samples_51(i,:)');
    correl_51(i) = m(2);
    
     m = corrcoef(double(states.StimulusType(index)),samples_56(i,:)');
    correl_56(i) = m(2);
    
     m = corrcoef(double(states.StimulusType(index)),samples_60(i,:)');
    correl_60(i) = m(2);
    
end

mean_11_stand = mean(samples_11_stand,2); %Find mean of each row
mean_11_odd = mean(samples_11_odd,2); %Find mean of each row

mean_51_stand = mean(samples_51_stand,2); %Find mean of each row
mean_51_odd = mean(samples_51_odd,2); %Find mean of each row

mean_56_stand = mean(samples_56_stand,2); %Find mean of each row
mean_56_odd = mean(samples_56_odd,2); %Find mean of each row

mean_60_stand = mean(samples_60_stand,2); %Find mean of each row
mean_60_odd = mean(samples_60_odd,2); %Find mean of each row

figure(11)
subplot(1,4,1)
plot(time,mean_11_stand','r')
hold on
plot(time,mean_11_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
ylabel('Signal Amplitude (mV)')
hold off

subplot(1,4,2)
plot(time,mean_51_stand','r')
hold on
plot(time,mean_51_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
hold off

subplot(1,4,3)
plot(time,mean_56_stand','r')
hold on
plot(time,mean_56_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
hold off

subplot(1,4,4)
plot(time,mean_60_stand','r')
hold on
plot(time,mean_60_odd','k')
xlim([0 800])
ylim([-2.5 4.5])
hold off

han=axes(figure(11),'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Time After Stimulus (ms)');

figure(12)
subplot(1,4,1)
plot(time,correl_11.^2,'k')
xlim([0 800])
ylim([0 0.03])
ylabel('r^2')

subplot(1,4,2)
plot(time,correl_51.^2,'k')
xlim([0 800])
ylim([0 0.03])

subplot(1,4,3)
plot(time,correl_56.^2,'k')
xlim([0 800])
ylim([0 0.03])

subplot(1,4,4)
plot(time,correl_60.^2,'k')
xlim([0 800])
ylim([0 0.03])

han=axes(figure(12),'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Time After Stimulus (ms)');

