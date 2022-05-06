%BME6360 Project 1
%Jacey Davies

clear;
close all;
clc;

%Loading the data
load Project1.mat;

%Creating the row vector t
t = [0:1/Fs:31/Fs]';

%Plotting the data
figure(1)
plot(t,X,'k-');
xlabel('t(ms)')
ylabel('Amplitude (mV)')

%Plotting the mean of the data with standard deviation error bars
figure(2)
shadedErrorBar(t,X,{@mean,@std});
xlabel('t(ms)')
ylabel('Amplitude (mV)')

%Applying PCA
[E,A,L]=pca(X);

%Plotting first 3 principle components
figure(3)
subplot(3,1,1);
plot(t,E(:,1),'k-');
ylabel('PC_1')

subplot(3,1,2);
plot(t,E(:,2),'k-');
ylabel('PC_2')

subplot(3,1,3);
plot(t,E(:,3),'k-');
ylabel('PC_3')
xlabel('t(ms)')

%Plotting the projection of the data onto the m-dimensional space (m=2)
figure(4)
scatter(A(:,1),A(:,2),100,'k','.');
xlabel('PC_1')
ylabel('PC_2')

%Plotting one observation and its reconstruction using this space
figure(5)
subplot(2,1,1)
p = 5; %Chose the 5th observation arbitrarily
plot(t,X(p,:),'k-') %Plotting the pth observation
ylabel('x_k(t)')

subplot(2,1,2)
M = mean(X,1);
X_est = M + A(p,1:2)*E(p,1:2)'; %Estimate using m = 2
plot(t,X_est,'k-')
ylabel('\^x_k(t)')
xlabel('t(ms)')

%Spike sorting the data using k-means clustering
A_2D = [A(:,1),A(:,2)];
idx = kmeans(A_2D,2); %k = 2 clusters

figure(6)
gscatter(A(:,1),A(:,2),idx,'rg','.',[],'off') %Scatterplot for k = 2
xlabel('PC_1')
ylabel('PC_2')

figure(7)
idx = kmeans(A_2D,3); %k = 3 clusters
gscatter(A(:,1),A(:,2),idx,'rgb','.',[],'off') %Scatterplot for k = 3
xlabel('PC_1')
ylabel('PC_2')

%Plotting the average spike for each cluster, with standard deviation shown
figure(8)
subplot(1,3,1)
shadedErrorBar(t,X(idx==1,:),{@mean,@std},'lineprops','r');
title(['n = ' num2str(length(idx(idx==1)))])
ylabel('Amplitude (mV)')
xlabel('t(ms)')

subplot(1,3,2)
shadedErrorBar(t,X(idx==2,:),{@mean,@std},'lineprops','g');
title(['n = ' num2str(length(idx(idx==2)))])
xlabel('t(ms)')

subplot(1,3,3)
shadedErrorBar(t,X(idx==3,:),{@mean,@std},'lineprops','b');
title(['n = ' num2str(length(idx(idx==3)))])
xlabel('t(ms)')

%Repeating the process with m = 3
%Plotting the projection of the data onto the m-dimensional space (m=3)
figure(9)
scatter3(A(:,1),A(:,2),A(:,3),100,'k','.');
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
grid on

figure(10)
subplot(2,1,1)
p = 5; %Chose the 5th observation arbitrarily
plot(t,X(p,:),'k') %Plotting the pth observation
ylabel('x_k(t)')

subplot(2,1,2)
M = mean(X,1);
X_est = M + A(p,1:3)*E(p,1:3)'; %Estimate using m = 3
plot(t,X_est,'k')
ylabel('\^x_k(t)')
xlabel('t(ms)')

%Spike sorting with k mean clustering
A_3D = [A(:,1),A(:,2),A(:,3)];
idx = kmeans(A_3D,3); %k = 3 clusters

figure(11)
scatter3(A(idx==1,1),A(idx==1,2),A(idx==1,3),100,'r','.') %Scatterplot for k = 3
hold on
scatter3(A(idx==2,1),A(idx==2,2),A(idx==2,3),100,'g','.')
scatter3(A(idx==3,1),A(idx==3,2),A(idx==3,3),100,'b','.')
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
grid on
hold off

figure(12)
idx = kmeans(A_3D,4); %k = 4 clusters
scatter3(A(idx==1,1),A(idx==1,2),A(idx==1,3),100,'r','.') %Scatterplot for k = 4
hold on
scatter3(A(idx==2,1),A(idx==2,2),A(idx==2,3),100,'g','.')
scatter3(A(idx==3,1),A(idx==3,2),A(idx==3,3),100,'b','.')
scatter3(A(idx==4,1),A(idx==4,2),A(idx==4,3),100,'m','.')
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
grid on
hold off

%Plot average spike for ea cluster
figure(13)
subplot(2,2,1)
shadedErrorBar(t,X(idx==1,:),{@mean,@std},'lineprops','r');
title(['n = ' num2str(length(idx(idx==1)))])
ylabel('Amplitude (mV)')

subplot(2,2,2)
shadedErrorBar(t,X(idx==2,:),{@mean,@std},'lineprops','g');
title(['n = ' num2str(length(idx(idx==2)))])

subplot(2,2,3)
shadedErrorBar(t,X(idx==3,:),{@mean,@std},'lineprops','b');
title(['n = ' num2str(length(idx(idx==3)))])
ylabel('Amplitude (mV)')
xlabel('t(ms)')

subplot(2,2,4)
shadedErrorBar(t,X(idx==4,:),{@mean,@std},'lineprops','m');
title(['n = ' num2str(length(idx(idx==4)))])
xlabel('t(ms)')


