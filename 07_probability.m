%% 07 Probability
%
% This tutorial covers:
% (1) ways to compute probabilities form counts, how to 
% (2) sampling from a small sample
% (3) sampling from a large population
% (4) Conditional probability
% (5) Addtition and multiplication rule
%
% Copyright Franco Pestilli Indiana University Spring 2015 K310

%% Probability from a sample
% Let's take a look at what happens when we draw from a data smaple
S = [1 2 2 3 3 3 4 4 4 4 5 5 5 6 6 7]; % This is our sample
n  = numel(S);

% From counts to probabilities:
%
% Let's compute the probability of each number in the sample.
% We choose the sampling spacing between the bins.
binWidth = 1;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(S):binWidth:max(S); 

% We use histogram to compute the shape of the distribution empirically from P
[y,x] = hist(S,bins);

% First we show this distribution which counts the occurrences of data
% in each bin
figure('name','Sample distribution','color','w')
subplot(1,2,1);
bar(x,y,'k')
title('Sample frequency distribution')
set(gca,'tickdir','out', 'box','off')
ylabel('Number of occurrences')
xlabel('Data value (a.u.)')
axis square

% Next, we trasform the distribution in probabilities. This is done by
% diving the counts by the total number of data points (N) and by
% accounting for the type of bins we used for computing the counts.
prob = y / n * binWidth;

subplot(1,2,2);
bar(x,prob,'k')
title('Sample probability distribution')
set(gca,'tickdir','out', 'box','off')
ylabel('P(D)')
xlabel('Data value (a.u.)')
axis square

%% Drawing from a sample
% (1) Drawing with replacement
% We will first draw a sample of 2 data points from the sample S the
% draw will be perfomed with replacement.
n = 2;
replacement = true; % Read more about true: "help true"
sample1 = randsample(S,n,replacement);

% We plot the sample on top of the distribution.
figure('name','Sample')
subplot(1,2,1)
plot(S,zeros(size(S)),'ro','linewidth',2)
hold on 
y = 0.5*linspace(1,size(sample1,2),size(sample1,2)); % We prepare a y that changes incremetally 
plot(sample1,y,'b+','linewidth',2)
ylabel('Arbitrary')
xlabel('Data value')
axis square
legend({'Data available in sample','Sampled values'},'box','off','location','northeastoutside')

% (2) Drawing without replacement
% Next we will draw 10 data points without replacement
replacement = false; % Read more about false: "help false"
sample2 = randsample(S,n,replacement);

% We make a plot.
subplot(1,2,2)
plot(S,zeros(size(S)),'ro','linewidth',2)
hold on 
plot(sample2,zeros(size(sample2)),'ko','linewidth',2)
y = 0.5*linspace(1,size(sample2,2),size(sample2,2)); % We prepare a y that changes incremetally 
plot(sample2,y,'b+','linewidth',2)
ylabel('Arbitrary')
xlabel('Data value')
axis square
legend({'Data available in sample','Removed Data','Sampled values'},'box','off','location','northeastoutside')

%% Conditional probability (Dependence)
% Dependence: What is the probability of drawing a 2 after drawing a 1?
S = [1 2 3 4 5];

% Let's compute the probabilities of drawing each number
[y,x] = hist(S,S);
p = y/length(S);

% Now we draw twice.
% First we draw '1' then we compute the probability of drawing 2.
S = S(2:end);

% Let's compute the probabilities of drawing each number
[y,x] = hist(S,S);
p = y/length(S);
P(2)

%% Dependence AND: What is the probability of drawing a 2 AND 1?
S = [1 2 3 4 5];

% Let's compute the probabilities of drawing each number
[y,x] = hist(S,S);
p = y/length(S);
p2 = p(2); % Probability of 2

% First we draw '2' then we compute the probability of drawing 1.
S = S([1 3:5]);
[y,x] = hist(S,S);
p = y/length(S);

prob = p2 * p(1)

%% Dependence XOR: What is the probability of drawing a 2 OR 1?
S = [1 2 3 4 5];

% Let's compute the probabilities of drawing each number
[y,x] = hist(S,S);
p = y/length(S);

% We use the summation rule, the two events are mutually exclusive
prob = p(2) + p(1)

%% Independent OR: What is the probability of drawing a 2 AND 1?
S = [1 2 3 4 5];

% Let's compute the probabilities of drawing each number
[y,x] = hist(S,S);
p = y/length(S);

% We use the multiplication rule, the two events are independet
prob = p(2) * p(1)

%% Dealing with a large population.
% Here after we will generate a population of data to use later to cumpute
% probabilities.

N  = 10000;
mu = 100; % Mean fo the population
sd = 10; % Standard deviation of the population
P  = mu + sd * randn(N,1);

% Let's compute the histogram of the population and trasform the population
% into a probability distribution.

% We choose the sampling spacing between the bins.
binWidth = .5;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(P):binWidth:max(P); 

% We use histogram to compute the shape of the distribution empirically from P
[Y,X] = hist(P,bins);

% First we show this distribution which counts the occurrences of data
% in each bin
figure('name','Population distribution','color','w')
subplot(1,2,1);
bar(X,Y)
title('Population frequency distribution')
set(gca,'tickdir','out', 'box','off')
ylabel('Number of occurrences')
xlabel('Data value (a.u.)')
axis square

% Next, we trasform the distribution in probabilities. This is done by
% diving the counts by the total number of data points (N) and by
% accounting for the type of bins we used for computing the counts.
PROB = Y / (N * binWidth);

subplot(1,2,2);
bar(X,PROB,'k')
title('Population probability distribution')
set(gca,'tickdir','out', 'box','off')
ylabel('P(D)')
xlabel('Data value (a.u.)')
axis square

% Drawing from a population
% Finally we will now draw from a full population using the same methods
% introduced above. This will be helpful later on to test hypotheses
% statistically.
%
% (1) Drawing with replacement
% We will first draw a sample of 10 data points from the population P the
% draw will be perfomed with replacement.
n = 10;
replacement = true; % Read more about true: "help true"
sample1 = randsample(P,n,replacement);

% We plot the sample on top of the distribution.
subplot(1,2,1)
hold on  % Please not one 'hold on' per plot, is this necessary
plot(sample1,zeros(size(sample1)),'ro')
subplot(1,2,2)
hold on  % Please not one 'hold on' per plot, is this necessary
plot(sample1,zeros(size(sample1)),'ro')

% (2) Drawing without replacement
% Next we will draw 10 data points without replacement
replacement = false; % Read more about false: "help false"
sample2 = randsample(P,n,replacement);

% We make a plot.
subplot(1,2,1)
hold on  % Please not one 'hold on' per plot, is this necessary
plot(sample2,zeros(size(sample2)),'b+');
subplot(1,2,2)
hold on  % Please not one 'hold on' per plot, is this necessary
plot(sample2,zeros(size(sample2)),'b+');
