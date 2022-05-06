%BME6360 Project 2
%Jacey Davies

clear;
close all;
clc;

%Loading the Project 2 data
load Project2.mat;

%Putting the neuronal spike times into separate vectors
neuron1 = neuron(1).times;
neuron2 = neuron(2).times;

%Plotting the rasterplot for neuron 1 
figure(1)
graph = [6 3 2 1 4 7 8 9];
for i = 1:8 %for all directions
indDir = find(direction==i);
numTrials = length(indDir);
subplot(3,3,graph(i));hold on;
for j = 1:numTrials
%Here you must find the spikeTimes 1sec before and after for each trial given direction i
spikeTimes = (neuron1(find(neuron1>go(indDir(j))-1 & neuron1<go(indDir(j))+1))-go(indDir(j)))';
n = length(spikeTimes);
plot([spikeTimes; spikeTimes], [ones(1,n)*j-1; ones(1,n)*j],'k-');
end
xlim([-1 1]);
ylim([0 numTrials]);
hold off;
end

han=axes(figure(1),'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Trial Number');
xlabel(han,'Time(sec)');

%Plotting the rasterplot for neuron 2 
figure(2)
graph = [6 3 2 1 4 7 8 9];
for i = 1:8 %for all directions
indDir = find(direction==i);
numTrials = length(indDir);
subplot(3,3,graph(i));hold on;
for j = 1:numTrials
%Here you must find the spikeTimes 1sec before and after for each trial given direction i
spikeTimes = (neuron2(find(neuron2>go(indDir(j))-1 & neuron2<go(indDir(j))+1))-go(indDir(j)))';
n = length(spikeTimes);
plot([spikeTimes; spikeTimes], [ones(1,n)*j-1; ones(1,n)*j],'k-');
end
xlim([-1 1]);
ylim([0 numTrials]);
hold off;
end

han=axes(figure(2),'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Trial Number');
xlabel(han,'Time(sec)');

%Creating a PETH for neuron 1
bin = .02; % 20 ms bins
win = 1; % 1s window
edgesPeri = -win:bin:win;

figure(3)
for i=1:8
 indDir=find(direction==i); %find trials in a given direction
 numTrials=length(indDir);
 psth = 0;
 subplot(3,3,graph(i));
 hold on;
 for j =1:numTrials
spikeTimes = (neuron1(find(neuron1>go(indDir(j))-1 & neuron1<go(indDir(j))+1))-go(indDir(j)))';
psth=psth+histc(spikeTimes,edgesPeri)'/numTrials/bin; %/numTrials/bin normalizes histogram
 end
 subplot(3,3,graph(i))
 bar(edgesPeri,psth)
 ylim([0 30])
 hold off;
end

han=axes(figure(3),'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Number of Spikes');
xlabel(han,'Time(sec)');

%Creating a PETH for neuron 2
figure(4)
for i=1:8
 indDir=find(direction==i); %find trials in a given direction
 numTrials=length(indDir);
 psth = 0;
 subplot(3,3,graph(i));
 hold on;
 for j =1:numTrials
spikeTimes = (neuron2(find(neuron2>go(indDir(j))-1 & neuron2<go(indDir(j))+1))-go(indDir(j)))';
psth=psth+histc(spikeTimes,edgesPeri)'/numTrials/bin; %/numTrials/bin normalizes histogram
 end
 subplot(3,3,graph(i))
 bar(edgesPeri,psth)
 ylim([0 30])
 hold off;
end

han=axes(figure(4),'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Number of Spikes');
xlabel(han,'Time(sec)');

%Finding the firing rate for 1 sec after go cue for neuron 1
spikeCount=zeros(8,1);

for i=1:8
 indDir=find(direction==i); %find trials in a given direction
 numTrials(i)=length(indDir);
    for j =1:numTrials(i)
      centerTime=go(indDir(j)); %to center on start of movement time
      allTimes=neuron1-centerTime; %center spike times
      %spikeCount(i)=spikeCount(i)+sum(allTimes>-1&allTimes<1); %for 2s around go time
      spikeCount(i)=spikeCount(i)+sum(allTimes>=0&allTimes<1); %for 1s after go time
    end  
%divide by the number of trials & bin size (1 s) for a mean firing rate
  spikeCount(i)=spikeCount(i)/numTrials(i);

end

%Graphing the tuning curve for neuron 1
figure(5)
dir_deg = [0 45 90 135 180 225 270 315];

%Line of best fit
coefficients = polyfit(dir_deg, spikeCount, 2); %binomial
xFit = linspace(min(dir_deg), max(dir_deg), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'r-');

%Actual plot
hold on;
scatter(dir_deg,spikeCount,'filled','k')
xlabel(['Direction of motion (' char(176) ')'])
ylabel('Mean firing rate [Hz]')
ylim([0 12])
xticks([0 45 90 135 180 225 270 315])
xlim([0 315])
hold off;

%Finding the firing rate for 1 sec after go cue for neuron 2
spikeCount=zeros(8,1);

for i=1:8
 indDir=find(direction==i); %find trials in a given direction
 numTrials(i)=length(indDir);
    for j =1:numTrials(i)
      centerTime=go(indDir(j)); %to center on start of movement time
      allTimes=neuron2-centerTime; %center spike times
      %spikeCount(i)=spikeCount(i)+sum(allTimes>-1&allTimes<1); %for 2s around go time
      spikeCount(i)=spikeCount(i)+sum(allTimes>=0&allTimes<1); %for 1s after go time
    end  
%divide by the number of trials & bin size (1 s) for a mean firing rate
  spikeCount(i)=spikeCount(i)/numTrials(i);

end

%Graphing the tuning curve for neuron 2
figure(6)
%Line of best fit
coefficients = polyfit(dir_deg, spikeCount, 2); %binomial
xFit = linspace(min(dir_deg), max(dir_deg), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'r-');

%Actual plot
hold on;
scatter(dir_deg,spikeCount,'filled','k')
xlabel(['Direction of motion (' char(176) ')'])
ylabel('Mean firing rate [Hz]')
ylim([0 12])
xticks([0 45 90 135 180 225 270 315])
xlim([0 315])
hold off;

