function zero_eight_chance_error
%% 08 Chance error
%
% This tutorial covers two topics:
% (A) Sampling from a population, and the relation between sample size and error in estimating
% parameters such as the mean.
%
% (A.1) We build a large distribution with given mean and STD.
% (A.2) We sample with different sample sizes and compute the mean.
%       The mean will be slightly wrong each time. How much wrong depends on the sample size.
% (A.3) We will explore the relation between sample size and chance error in estimating the mean.
%       This is a montecarlo approach.
%
% (B) A box model of gambling. This is a simple version of the model in
%     Freeman et al., 2007. For this topic we will learn how to make a
%     function in MatLab to re-use the function to test different gambling
%     schemes and different outcomes form gambling.
%
% (B.1) Reproduce Kerrich's experiment (10,000 tosses of a coin).
% (B.2) Expected value and standard error
% (B.3) Show how to implement the central limit theorem.
%
% Copyright Franco Pestilli Indiana University Spring 2015 K310 

%% (A) Drawing from a population
% Here after we will generate a population of data to use later to compute
% probabilities.
N  = 1000;
mu = 100; % Mean fo the population
sd = 10;  % Standard deviation of the population
P  = mu + sd * randn(N,1);

% Let's compute the histogram of the population and trasform the population
% into a probability distribution.
%
% We choose the sampling spacing between the bins.
binWidth = .5;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(P):binWidth:max(P); 

% We use histogram to compute the shape of the distribution empirically from P
[Y,X] = hist(P,bins);

% Next, we trasform the distribution in probabilities. This is done by
% diving the counts by the total number of data points (N) and by
% accounting for the type of bins we used for computing the counts.
PROB = Y / (N * binWidth);

figure('name','Population distribution','color','w')
bar(X,PROB,'k')
title('Population probability distribution')
set(gca,'tickdir','out', 'box','off')
ylabel('P(D)')
xlabel('Data value (a.u.)')
axis square

%% (A.1) Drawing from a population
% Next we will draw from a full population using with replacement. We want
% to measure how well the sample represents the population. To do that we
% will use the mean. How well does the mean fo the sample represents that
% of the population?
%
% We will draw a sample of 3 data points from the population P the
% draw will be perfomed with replacement.
n = 3;
replacement = true; % Read more about true: "help true"
sample1 = randsample(P,n,replacement);
% Compute the mean this mean will be compared to the mean of the population
% from which the sample comes (mu):
m_sample1 = mean(sample1);

% We plot the sample on top of the distribution.
figure('name','Population distribution','color','w')

% We plot first the full-population
bar(X,PROB,'k')
title('Population probability distribution')
set(gca,'tickdir','out', 'box','off')
ylabel('P(D)')
xlabel('Data value (a.u.)')
axis square
hold on
% We marke the mean fo the population:
plot(mu,0,'r^','markersize',20,'markerfacecolor','r')
plot([mu mu],[0 max(PROB(:))],'r--')
% Then we mark the sampled values
plot(sample1,zeros(size(sample1)),'bo','markersize',12)
% Finally the mean of the smaple
plot(m_sample1,zeros(size(sample1)),'b^','markerfacecolor','b','markersize',12)

% (A.2) Now let's see how well we estimate the mean if we sample more
n = 300;
sample2   = randsample(P,n,replacement);
m_sample2 = mean(sample1);
% Then we mark the sampled values
% We omit the sample so not to crowd the plot, but if you want to see the
% sample use the following:
% plot(sample2,zeros(size(sample2)),'co','markersize',12)
% Finally the mean of the smaple
plot(m_sample2,zeros(size(sample2)),'c^','markerfacecolor','c','markersize',12)

% (A.3) To really evaluate the accuracy of the estimate of the mean from samples
% of different size we would want to:
% (*)  Sample with a certain numerosity
% (**) Compute the mean
% (***)Repeat the process many times to evaluate how reliable the estimates
%      are.
n = [3, 10, 30, 100, 200, 300]; % We will evaluate different sample sizes
m = 100; % To evaluate the 'goodness' of each sampel size in estimating 
         % the true mean we want to repeat the process many times.
for im = 1:m
    for in = 1:length(n)
        sample_tmp = randsample(P,n(in),replacement);
        m_sample(in,im) = mean(sample_tmp);
    end
end

% Now let's visualize the results of our experiment:
figure('name','Reliability of mean with different samples','color','w')
% We compare to the true mean
plot([0 max(n)],[mu mu],'r--')
hold on
plot(n,m_sample,'bo')
ylabel('Estimated mean of the population')
xlabel('Sample size used for estimates')
set(gca,'tickdir','out', 'box','off')
legend({'True mean of the population'},'box','off')

%% (B.1) Box model.
% This file is a function. You can see this by looking at the very first
% line of the file. It is defined there as a matlab 'function.'
% Functions allows to:
% - A convenient way to be called over and over.
% - They have inputs and outputs.
% - They do not share computations or variables with the rest fo the MatLab
%   environemnt.
%
% Below, at the end of this file we have difined another function. This
% function called boxGamblingModel allows creating a model similar to
% the ones introduced in Chapter 16 of the book (Freeman, Pisani, Purves,
% Statistics 2007, 4th Edition).

% The funcetion is defined below. It requires some INPUTS:
% valuesInBox - these are all the available values in the box we are
%               building, for example, [1 2 3].
% numberOfValues - These are the number of repeated values. For example:
%                  [3,5,2] woudl mean [1 1 1 2 2 2 2 2 3 3].
% numberOfDraws - This indicates the number of times that we will draw out
%                 of the box.
%
% The function returns some OUTPUTS:
% drawnValues - This variables contains all the values drawn from the box.
% box - The second output is the full 'box' created in the process.
keyboard
valuesInBox    = [1 2 3 4 5];
numberOfValues = [2 2 2 2 2];
numberOfDraws  = 3;
[drawnValues, box] = boxGamblingModel(valuesInBox,numberOfValues,numberOfDraws);

% Now that we have a model (the function boxGamblingModel) let's repeat
% the experiment that John Kerrich performed. We will toss a coin many
% times from 0 to 10,000 and count the number of heads:
valuesInBox    = [1 0]; %'1'=head, '0'=tail 
numberOfValues = [1 1]; % A coins has 1 head and one tail.
for iToss = 1:10000
    draws{iToss} = boxGamblingModel(valuesInBox,numberOfValues,iToss);
end

% Now lets compute the number of heads, this can be done by counting the
% heads, meanign the '1' or more simply by summing the values inside draws.
for iToss = 1:10000
    nHeads(iToss) = sum(draws{iToss});
end

% We will Plot Figure 1 in Chapter 16:
figure('name','Kerrich''s coin tossing experiment','color','w')
subplot(2,1,1)
x = 1:10000; % THe number of tosses
y = nHeads - (x./2); % The number of observed heads minus half the number of tosses
semilogx(x,y,'r-')
hold on
semilogx([1 10000],[0 0],'k--')
title('Chance error')
xlabel('Number of tosses')
ylabel(sprintf('Number of heads minus\nhalf the number of tosses'))
set(gca,'tickdir','out', 'box','off')

subplot(2,1,2)
title('Chance error')
x = 1:10000; % The number of tosses
y = ((nHeads./x) - 0.5)*100; % The number of observed heads minus half the number of tosses
semilogx(x,y,'r-')
hold on
semilogx([1 10000],[0 0],'k--')
title('Chance error as percentage')
xlabel('Number of tosses')
ylabel('Number of heads minus 50%')
set(gca,'tickdir','out', 'box','off')

%% (B.2) The expected value
% The expected value for a sum of draws made at random with replacement
% from a box equals: (number of draws) x (average of box)
%
% We first set up a model for computing the expected value:
valuesInBox    = [1 2 3 4 5 6]; % for example a 6-face dice 
numberOfValues = [1 1 1 1 1 1]; % A coins has 1 head and one tail.
numberOfTosses = 1:100;         % We will repeat the experiment with different tosses
for iToss = numberOfTosses
    [d{iToss}, box] = boxGamblingModel(valuesInBox,numberOfValues,numberOfTosses(iToss));
    
    % Now lets compute the expected value for each experiment.
    ev(iToss) = (numberOfTosses(iToss)) * (mean(box));
    
    % We also compute the observed value, the sum of the elements shoudl be
    % at around
    ov(iToss) = sum(d{iToss});
end

% We plot expected value and observed value:
figure('name','Expected and observed values','color','w')
semilogx(numberOfTosses,ev,'r-')
hold on
semilogx(numberOfTosses,ov,'b-')
xlabel('Number of tosses')
ylabel('Observed and expected value')
legend({'Expected value','Observed value'},'box','off')
set(gca,'tickdir','out', 'box','off')

%% (B.2) The expected value and standard error
% The expected value for a sum of draws made at random with replacement
% from a box equals: (number of draws) x (average of box)
%
% We first set up a model for computing the expected value:
valuesInBox    = [1 2 3 4 5 6]; % for example a 6-face dice 
numberOfValues = [1 1 1 1 1 1]; 
numberOfTosses = 1:100;         % We will repeat the experiment with different tosses
for iToss = numberOfTosses
    [d{iToss}, box] = boxGamblingModel(valuesInBox,numberOfValues,numberOfTosses(iToss));
    
    % Now lets compute the expected value for each experiment.
    ev(iToss) = (numberOfTosses(iToss)) * (mean(box));
    
    % We also compute the observed value, the sum of the elements shoudl be
    % at around
    ov(iToss) = sum(d{iToss});
end

% We plot expected value and observed value:
figure('name','Expected and observed values','color','w')
semilogx(numberOfTosses,ev,'r-')
hold on
semilogx(numberOfTosses,ov,'b-')
xlabel('Number of tosses')
ylabel('Observed and expected value')
legend({'Expected value','Observed value'},'box','off')
set(gca,'tickdir','out', 'box','off')

% We clear some variables before we move on
clear d ov ev 

% The standard error is defined as:
% (sqrt(number of draws)) * (std(box))
%
% We first set up a model for computing the standard error:
valuesInBox    = [1 2 3 4 5 6]; % for example a 6-face dice 
numberOfValues = [1 1 1 1 1 1]; 
numberOfTosses = [1 10 50 100];    % We will repeat the experiment with 
                                % different tosses.
nRepeats    = 50;              % We will repeat the process several times 
                                % to compute a distirbution of observed
                                % values.
for iToss = 1:length(numberOfTosses)
    % We repeat each toss many times to create a distribution
    for ib = 1:nRepeats
        [d{iToss,ib}, box] = boxGamblingModel(valuesInBox,numberOfValues,numberOfTosses(iToss));
        
        % We also compute the observed value, the sum of the elements shoudl be
        % at around
        ov(iToss,ib) = sum(d{iToss,ib});
    end
    
    % Now lets compute the expected value for each experiment.
    ev(iToss) = (numberOfTosses(iToss)) * (mean(box));

    % Now lets compute the expected value for each experiment.
    se(iToss) = sqrt(numberOfTosses(iToss)) * (std(box));    
end

% We plot expected value and observed value:
figure('name','Standard error, observed and expected values','color','w')
semilogx(numberOfTosses,ev,'r-','linewidth',2)
hold on
semilogx(numberOfTosses',ov,'bo')
semilogx([numberOfTosses;numberOfTosses], ...
         [ev-2*se;       ev+2*se],'r-','linewidth',2)
xlabel('Number of tosses')
ylabel('Observed and expected value')
set(gca,'tickdir','out', 'box','off')


end % Main function

function [drawnValues, box] = boxGamblingModel(valuesInBox,numberOfValues,numberOfDraws)
%
% This is a general model for building gamblings and draws out of boxes.
%
% INPUTS:
% valuesInBox - these are all the available values in the box we are
%               building, for example, [1 2 3].
% numberOfValues - These are the number of repeated values. For example:
%                  [3,5,2] would mean [1 1 1 2 2 2 2 2 3 3].
% numberOfDraws - This indicates the number of times that we will draw out
%                 of the box.
%
% OUTPUTS:
% drawnValues - This variables contains all the values drawn from the box.
% box - The second output is the full 'box' created in the process.
%
% Copyright Franco Pestilli Indiana University Spring 2015 K310 

% We build the 'box' byt reading each value in 'valuesInBox' and repeating 
% it as many times requested in 'numberOfValues'
box = [];
for iv = 1:length(valuesInBox)
    box = [box, repmat(valuesInBox(iv),1,numberOfValues(iv))]; 
end

% We now draw (with replacement) out of the box as many samples as 
% requested in numberOfDraws (on number per draw).
drawnValues = randsample(box,numberOfDraws,true);

end
