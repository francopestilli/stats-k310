%% Creating a variables and visualizing data
%
% - This tutorial will show how to plot a histogram with areas and bar-heights.
% - How compute the mean and the median and how to use MatLab corresponding
%   Functions.
% - How to compute the r.m.s. (root-mean-square) value of a sample.
% - How to compute the standard deviation and use the corresponding matlab
%   function.
% 
% Pestilli, Franco K310 Spring 2016 Indiana University Bloomington

%% (0) Simulating data.
% We create a vector of random data points. The values will be
% randomly drawn from a gaussian distribution with a mean of 0 and a
% standard deviation of 1.  To do so we will use a built in MatLab
% function:
rnd_data = 2+1*randn(100,1); % Type: 'doc randn' for more information 

%% (1) Plotting a histogram with equal bins.
% The height of each bar will indicate the number of elements in the sample
% that are part of the interval indicated by the bar horizontal width.
figure('name','test','color','w')
subplot(3,1,1); % Subplot plots multiple axis in a single matlab figure
xbins1 = -6:6; 
hist(rnd_data,xbins1);
hold on
plot(rnd_data, 20*ones(size(rnd_data)),'go','markersize',14)
set(gca,'tickdir','out', 'box','off','xlim',[-8 8], ...
        'xtick',[-8:2:8],'fontsize',25)
ylabel('Frequency','fontsize',25)

%% (2) Plotting a histogram with unequal bins.
subplot(3,1,2);  % Subplot plots multiple axis in a single matlab figure
xbins2 = [-4 -2 0 .5 3];
hist(rnd_data,xbins2)
hold on
plot(rnd_data,zeros(size(rnd_data)),'r+')
set(gca,'tickdir','out', 'box','off','xlim',[-8 8], ...
        'xtick',[-8:2:8],'fontsize',25)
ylabel('Frequency','fontsize',25)

%% (3) Plotting a histogram by computing probabilities:
[y,x] = hist(rnd_data,xbins1);
% Divide by the sum: 
% this does not scale to probability
subplot(3,1,3);  % Subplot plots multiple axis in a single matlab figure
% Divide by the area:
% This method scales the y-axis to probability
bar(x,y/trapz(x,y));
set(gca,'tickdir','out', 'box','off','xlim',[-8 8],'xtick',[-8:8])
ylabel('Probability')
xlabel('Value of sample data')
set(gca,'tickdir','out', 'box','off','xlim',[-8 8], ...
        'xtick',[-8:2:8],'fontsize',25)

%% (3) Build two data sets sampling from different distributions.
% Here after we will plot two distributions with different shape.
% One shape is symmetric (same number of values to the left and right of
% the mean) the other shape is not symmetric, it is skewed. The skewed
% distribution is a Chi^2 distribution. When we summarize both
% distributions with mean, median and STD the significance of each measure
% changes. Both mean and median describe well the symmetric data drawn from
% a normal (Gaussian) distirbution but not so well the dat adrawn from the
% Chi^2 distribution.
%
chirnd_data = chi2rnd(5,6000,1); % Please explore the help of this new function
rnd_data    = 6+2*randn(6000,1); % Type: 'doc randn' for more information 

% We plto the first sample, not symmetric.
figure('name','Two histograms','color','w')
subplot(2,1,1); % Subplot plots multiple axis in a single matlab figure
[y,x] = hist(chirnd_data,60);
bar(x,y,'k')
hold on
set(gca,'tickdir','out', 'box','off','xlim',[0 20],'xtick',[0:5:20],'fontsize',25)
ylabel('Frequency','fontsize',25)

% Add the mean
m  = mean(chirnd_data);
sd = std(chirnd_data);
plot([m-sd m+sd],[20 20],'r-','linewidth',2)
plot([m m],[20 20],'ro','markerfacecolor','r','markersize',10)

% Add the median
m = median(chirnd_data);
plot([m m],[15 15],'b^','markerfacecolor','b','markersize',10)

% We plot the second sample, symmetric.
subplot(2,1,2); % Subplot plots multiple axis in a single matlab figure
[y,x] = hist(rnd_data,60);
bar(x,y,'k')
hold on
set(gca,'tickdir','out', 'box','off','xlim',[0 20],'xtick',[0:5:20],'fontsize',25)
ylabel('Frequency')
xlabel('Value of sample data')

% Add the mean
m  = mean(rnd_data);
sd = std(rnd_data);
plot([m-sd m+sd],[20 20],'r-','linewidth',2)
plot([m m],[20 20],'ro','markerfacecolor','r','markersize',10)

% Add the median
m = median(rnd_data);
plot([m m],[15 15],'b^','markerfacecolor','b','markersize',10)


%% (4) Plot a bimodal dataset. Use mean, median and STD to summarize the sample.
% Here after we will create two samples from different distributions,
% meaning sampling from distributions with different means and standard
% deviations. We will then combine the samples into a single distrbution of
% and compute mean and standard deviation of this new distribution.
% We show that the mean and the standard deviation do not represent well
% the data in this 'bimodally' distributed set.
% 
rnd_data1 = 0+1*randn(2000,1);
rnd_data2 = 5+2*randn(2000,1); % Type: 'doc randn' for more information 
rnd_data = [rnd_data1; rnd_data2];
[y, x] = hist(rnd_data,200);
y = y/trapz(x,y); % Trasform in probabilities

figure('name','Bimodal dataset','color','w')
bar(x,y); hold on
set(gca,'tickdir','out', 'box','off','xlim',[-10 10],'xtick',[-10:5:10],'fontsize',25)
ylabel('Frequency','fontsize',25)
xlabel('Value of sample data','fontsize',25)

% Add the mean
m  = mean(rnd_data);
sd = std(rnd_data);
plot([m-sd m+sd],[.1 .1],'r-','linewidth',2)
plot([m m],[.1 .1],'ro','markerfacecolor','r','markersize',30)

% Add the median
m = median(rnd_data);
plot([m m],[.12 .12],'b^','markerfacecolor','b','markersize',30)


%% (5) Effects of outliers on the mean, median and STD of symmetric distribution.
% Notice that when we add outliers to the data the mean does not accurately
% describe the central tendency of the data and the standard deviation does
% not accurately reflect the spread in the data.

% Create a new dataset with the same values as rnd_data and a few outliers
% horzcat simply adds the new values at the end of the vector d0sorted
rnd_data   = randn(1,120); % Sample without OUTLIERS 
rnd_data_o = horzcat(rnd_data, 60, 50, 60);  % Add outliers, 60, 50 and 60 

% Compute mean, median and STD, for the two samples.
m1   = mean(rnd_data_o);
sd1  = std(rnd_data_o);
med1 = median(rnd_data_o);
m   = mean(rnd_data);
sd  = std(rnd_data);
med = median(rnd_data);

% Now plot the results in a new figure
figure('name','Adding outliers','color','w')
hold on;
[y,x]   = hist(rnd_data,14);
[y1, x] = hist(rnd_data_o,x);
binmax  = max([y,y1]); % Let's find the maximum height of the graph, 
                      % here instead of horzcat we simply use [] to
                      % concatenate.
plot(x,y, 'k-', 'linewidth',6);
plot(x,y1,'r-','markerfacecolor','k', 'linewidth',2, 'markersize',12);
legend({'Distribution without','Distribution with outliers'},'Location','SouthWest')

% Means are circles
plot(m,   binmax, 'ko', 'markersize',20,'markerfacecolor','k');
plot(m1,  binmax, 'ro', 'markersize',20,'markerfacecolor','r');

% Medians are triangles.
plot(med, binmax-3,'k^','markerfacecolor','k', 'markersize',20);
plot(med1,binmax-3,'r^','markerfacecolor','r', 'markersize',10);
ylabel('Frequency')
xlabel('Value of data in sample')

% Take a look at all the many things we can set on a graph to make it
% appealing.
set(gca,'xlim',[-3 3], 'tickdir','out', ...
    'ytick',[0 10 20 25], ...
    'yticklabel',{'0', '10', '20', ''}, ...
    'xtick',[-2, 0, 2], ...
    'fontsize',20)
    
% End