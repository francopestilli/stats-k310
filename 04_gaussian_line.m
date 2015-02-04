%% (6) The Gaussian distribution.
% Next we will compare the gaussin distribution computed numerically using
% randn.m with the gaussian distribution computed analytically. Using an
% explicit equation fo the gaussian. 
%
% We will introduce the 'for' loop.
%
% The equation for a gaussian is
%     
%    g(x) = (1 / (s*sqrt(2*pi))) *exp(-((x-m).^2)/(2*s^2))
% 
% - m - is the mean and 
% - s - is the standard deviation.  
% 
% We can plot the Gaussian with the appropriate mean
% and standard deviation on both histograms and see how well they fit the
% data.
nSamples = 200;

% Create a new dataset with the same values as rnd_data and a few outliers
% horzcat simply adds the new values at the end of the vector d0sorted
rnd_data   = randn(1,nSamples); % Type: 'doc randn' for more information 
m   = mean(rnd_data);
sd  = std(rnd_data);
med = median(rnd_data);

% To plot the Gaussins we need to create a vector of values over which we
% will evaluate the function. These values must span the minimum and max
% value in the data we want to describe. We will (arbitraily) space the
% values over whcih we evaluate the gaussian by an interval of 0.01.
dx = 0.01;
g.x = [min(rnd_data):dx:max(rnd_data)]; % Here we create a vector 
                                        % of linearly increasing values.
% Note that here we are introducing a structure. doc it!

% Next we will use a 'for loop'. 
% FOR loops are fundamental in programming. They allow the computer to
% repeate an operation many times. A decided number of times. They also
% allow reading inside long vectors at different locations. FOr exampel
% herefater we will Read x at each entry and evaluate the Gaussian function
% at each x.
%
% To do so we loop over the x values and calculate the p(x) value returned
% by the gaussian using the equation above.  This is done by plugging each
% x value into the equation one at a time.
% 
% Note that we are using sd and m from the data set above called
% rnd_data.
for ii = 1:length(g.x) % For each element in x perform the operations before 'end'
    g.y(ii) = (1 / (sd*sqrt(2*pi))) * exp(-((g.x(ii)-m)^2)/(2*sd^2));
end

% In matLab many operations can be performed avoiding FOR loops. This is
% why we have not encountered a for loop yet. The following is the way you
% woudl perform the operation without a FOR loop.
% 
% The dot operators perform arithmetic element-wise and they are useful for
% writing vectorized code.  So remember .* is different from *
g.y = (1 / (sd*sqrt(2*pi))) * exp(-((g.x-m).^2)/(2*sd^2));

% Now we will open a new figure window and plot the gaussian over the
% histogram
figure('name','Gaussian distribution', 'color','w'); 
hold on;

% First we plot the histogram, but the y-axis is now in units of relative
% fractions instead of absolute counts.  by doing so, we can now plot the
% actual Gaussian probability distribution on top of the histogram and it
% will comparable in scale.
[h.y,h.x] = hist(rnd_data,round(nSamples*.1)); % We are introducing 'round.m' doc it! 

% The Gaussian is a probability distribution. To plor the data and the
% Gaussin on top of eachothers we need to normalize the histogram so that
% the heigth does nto represent 'frequecny' but 'probability.' To do so, we
% compute the total area of the histogram by multiplying the bin width by
% the height of the bins.

% We divide by the area of the bins (black in the figure)
% This will compute the total area of the histogram by multiplying the bin
% width by the height of the bins.
area    = diff(h.x(1:2)) * sum(h.y);
h.yprob = h.y/area;

% Now when we plot the histogram bins we scale them by the total area
% of the histogram to rescale the histogram to have an area of 1.
% This is an essential feature of a probability distribution. The total
% probability for all the events the disgribution describes is '1' or '100'
% probability, this means that the area under the curve is 1.
bar(h.x,h.yprob,'k');

% Now we plot a line with the calculated height of the gaussian, y, for
% each value, x.
plot(g.x,g.y,'-r','linewidth',4);
ylabel('Frequency')
xlabel('Value of data in sample')
set(gca,'xlim',[-3 3], 'tickdir','out', ...
    'ytick',[0 .2 .4 .6], ...
    'xtick',[-2, -1, 0, 1, 2], ...
    'fontsize',20)

%% (7) Now we will plot the gaussian for a sample with outliers:
% How does the gaussian represent this sample with outliers?
nSamples = 200;

% Generate gaussian random data
rnd_data   = randn(1,nSamples); % Type: 'doc randn' for more information 
rnd_data_o = horzcat(rnd_data, 60, 50, 60); 
m   = mean(rnd_data_o);
sd  = std(rnd_data_o);
med = median(rnd_data_o);

% Generate the Gaussian
dx = 0.01;
g.x = [min(rnd_data_o):dx:max(rnd_data_o)]; 
g.y = (1 / (sd*sqrt(2*pi))) * exp(-((g.x-m).^2)/(2*sd^2));
[h.y,h.x] = hist(rnd_data_o,round(nSamples*.1)); % We are introducing 'round.m' doc it! 
h.area    = diff(h.x(1:2)) * sum(h.y);
h.yprob = h.y/h.area;

% Plot
figure('name','Gaussian distribution outliers', 'color','w'); 
hold on;
bar(h.x,h.yprob,'k');
plot(g.x,g.y,'-r','linewidth',4);
ylabel('Frequency')
xlabel('Value of data in sample')
set(gca,'xlim',[-10 60], 'tickdir','out', ...
    'ytick',[0 .1 .2], ...
    'xtick',[-2, 0, 2:10:60], ...
    'fontsize',20);

%% Plotting a line
% Here ater we will plot a line and make a few examples of their slopes

% Define the equation of a line
x = linspace(-5,5,10);% Please study the help of the new function 'linspace'
slopes = [2 1 1.5];
intercepts = [2 1 2];

% We loop over each slope and intercept and compute the values for the
% line. We store them into a variabile and plot them later on
for ii = 1:length(slopes)
    y(ii,:) = slopes(ii).*x + intercepts(ii);
end

% We now prepare colors and line properties for each line
% First we use a structure to store in a compact for all the formatting
% information. Color and lineWdth.
%
% A structure is "a variables that contains variables" like 'property' below
% contains 'color' and 'linewidth.' The structure I create below contains
% two 'fields' color and linewidth. Color is a cell array, which means it
% can contains multiple types of variables, numbers or characters, for
% example. The filed linewdth is a numeric vector. How do you recognize a
% cell array from a numeric array? The curly parenthesis {} vs. [].
property.color     = {'r-','b-','g-'}; 
property.linewidth = [2 1 4];

% We open a window where we will plot the lines
figure('name','Lines!','color','w');
hold on

% We plot one line at the time
for ii = 1:length(slopes)
    plot(x,y(ii,:),property.color{ii},'linewidth',property.linewidth(ii))
end

% Lets add axis going through the center of the coordiante frame (0,0)
plot([0 0],[-10 10],'k--')
plot([-10 10],[0 0],'k--')
set(gca,'tickdir','out','box','off','fontsize',14)
legend({'Red line','Blue line','Green line'},'box', 'off','location','southeast');

