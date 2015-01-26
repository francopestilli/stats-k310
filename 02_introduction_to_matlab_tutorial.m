%% Creating a variables and visualizing data
%
% This tutorial will show how to define variables, create a (simulated)
% data set, plot it.
% 
% Pestilli, Franco K310 Spring 2014 Indiana University Bloomington

%% (1) Defining a variable in the work space.
% Please copy and paste the following lines inside the matlab command
% window. One at the time. Explore what changes int he workspace after each
% operation.
var = 2;         % We define a simple variable and assign an arbitrary value to it. 
var2 = [];       % We define a variable with no value.
var2 = var+2;    % We perform a mathematical operation on the variable.
prod = var2*var; % We perform a second operation on the variable. We multiply it by the first variable. 

% So Variables can be defined, they can have values and once a value is
% assigned to the variable the variable can be used in place of that value.
% You do not have to type 2 again jus use the variable. This simple
% property turns out to be very helpful when the number of lines and the
% complexity of your code increase. Imaging tracking a series of '2' and
% '3' in a 1,000,000 lines script.
%
% The final variable is called 'prod' this turns out to be a mistake. First
% check the value of the variable type:
whos % This will return all the variables in the workspace with associated information. 
% Then type the variable at command line:
prod
% MatLab will return its value

% Next let's delete the variable:
clear prod 
% clear.m is 'builtin' MatLab function that deletes variables type 'doc
% clear' for more information Next type:
help prod
% This will return the help for a function called 'prod.m' MatLab by
% default overwrites functions with variables. So if a variables takes a
% name of a function that exists in the MatLab path access to the function
% is temporarily lost. Let's try again.
prod = 2;
help prod
clear prod
help prod
prod(2,2)

%% (2) Variable types
%  Matlab has many types of variables
strings = 'this is a string of characters';
vectors = [1 2 3];
arrays = [1 2 3; 4 5 6; 7 8 9];

%% (3) MatLab commands:
% MatLab has many basic commands. Most of them perform mathematical
% operations. So far we only used a couple, later on we will encounter many
% more. For example:
2+1 % Addition
2-1 % Subtraction
2*3 % Multiplication
4/2 % Division
2^2 % Power
% Type doc * form more information about the basic MatLab operators. Some
% operators are specific for vectors and arrays. For example:
arrays.^2

%% (4) A first example: Simulating data.
%  Next we will step a little bit ahead and use a combination of operations
%  and builtin functions to generate 'data.' Some of this might feel
%  confusing at first but it will help to break the ice.

% (4.1) We create a vector of random data points. The values will be
% randomly drawn from a gaussian distribution with a mean of 0 and a
% standard deviation of 1.  To do so we will use a built in MatLab
% function:
rnd_data = randn(1,200); % Type: 'doc randn' for more information 

%% (4.2) We will now show visualize the data. We will plot a histogram of
% the data by summing within 20 equidistant bins.
figure(1); hist(rnd_data, 20);

% (4.3) We will calculate the mean and standard deviation of the data
m   = mean(rnd_data);
sd  = std(rnd_data);

%% (4.4) Now we will plot the location of the mean +/- 1 standard deviation
% of the data.
hold on; % This function tells MatLab to keep focussing all the plots in the current figure, instead of opening a new one.

% (4.5) Plot the location of the mean
plot(m,20,'ko','markerfacecolor','r', 'markersize',10);

% (4.6) Plot the location of the two standard standard deviations
sd_minus = m-sd;
sd_plus = m+sd;
plot([sd_minus, sd_plus],[20, 20], 'r-', 'linewidth', 2);

%% (4.7) Now let's create two different second data sets. 
% The first one will have mean 2 and standard deviation 3. 
% The second one will have mean 0 and standard deviation 2.
% Please note how the mean 'adds' to the values and
% the standard deviation 'multiplies' the values
rnd_data1 = 2 + 3*randn(1,2000); 
rnd_data2 = 0 + 2*randn(1,2000); 

%% (4.8) Plot data set 1 in a new figure with mean and standard deviation. 
% Please not that this time we are using short cuts and we are never saving out
% variables for the mean and the standard deviation. We are computing them
% on the fly as I build the plotting command.
figure(2); 
hist(rnd_data1,20)
hold on
% Plot the mean
plot(mean(rnd_data1),20,'ko','markerfacecolor','r', 'markersize',10);
% Plot the standard deviation
plot([mean(rnd_data1)-std(rnd_data1), mean(rnd_data1)+std(rnd_data1)],[20, 20], 'r-', 'linewidth', 2);

% Try doc plot to read more about how to use the function. This is one of
% the fundamental functions in MatLab we will use it many times.

%% (5) More advanced plotting.
% We would like to compare the two data sets. For this, we would like to
% plot the two histograms in the same figure. Here is a way to do that. We
% are going to make a new figure.
% We 'prepare' the histograms of the data sets.
[y,x] = hist(rnd_data1,20); 
% We use hist.m in a slightly different way. When called with outputs
% hist.m changes behavior and it does not plot any longer. It computes a
% histograms returns it into variables but does not plot. This can be a
% little bit annoying in MatLab things behave differently depending on how
% you act with them.
[y2,x2] = hist(rnd_data2,x); 
% We used histogram in an even new way. Instead of requesting 20 bins as we
% did before, we passed the bins directly (x). This will allow us to
% measure the values in rnd_data and rnd_data2 at the same positions. Now
% we are going to plot the two histograms in the same figure with two
% different colors:
figure(3);
set(gcf,'name','Two data sets'); % Here I change the name of the figure to something informative.
plot(x,y,'r-','color',[.8 .5 .2],'linewidth',2);
hold on
plot(x2,y2,'b-','color',[.2 .4 .9],'linewidth',5)
% We plot their respective means
plot([mean(rnd_data1),mean(rnd_data1)],[10,10],'ro-', ...
    'color',[.8 .5 .2],'markerfacecolor',[.8 .5 .2]);
plot([mean(rnd_data2),mean(rnd_data2)],[20,20],'bo-', ...
    'color',[.2 .4 .9],'markerfacecolor',[.2 .4 .9]); 
% Please notice that here we are calling plot in a new way.

% Plot the standard deviations
plot([mean(rnd_data1)-std(rnd_data1), mean(rnd_data1)+std(rnd_data1)],[10, 10], ...
    'r-', 'linewidth', 2, 'color',[.8 .5 .2]);
plot([mean(rnd_data2)-std(rnd_data2), mean(rnd_data2)+std(rnd_data2)],[20, 20], ...
    'r-', 'linewidth', 2, 'color',[.2 .4 .9]);

% Next we will define the axis. This is very important and now on we will
% *always* define labels for the axis.
ylabel('Number of occurrencies')
xlabel('Value of random number generated')
% We have learned how to create vectors of random numbers to simulate
% distributions of data. We have learned also how to create histograms and
% plot them. Finally we have learned that many MatLab functions can be
% called with different inputs and depending on the inputs they can behave
% slightly differently.
