%% 10 Hypotheses testing
% 
% This matlab tutorial is intended to complement lecture #2.  
%
% It will show how to:
%
% 1. Explore datasets with one variable and two conditions
% 2. Implement nonparametric alternatives to the t-test
% 3. Explore datasets with two variables and one condition
% 4. Compute the Pearson correlation coefficients
% 5. Establish the error of the estimate, using bootstrap.
% 6. Test hypotheses on the correlation between two variables, using
%    bootstrap.
% 7. Use Montecarlo Methods to test hypotheses.
% 
%
% Copyright Franco Pestilli Indiana University Spring 2015 K310 
% 
% This tutorial is modified from a tutorial Franco Pestilli wrote for
% Class Psych 216A Spring 2012 Stanford University

%% 1. Exploring a complex dataset: one variable, two conditions
%
%  Motivation: 
%    Samples from a population will show variation in their means. So, when
%    drawing two samples from a single population, there is some chance
%    that the samples will have different means, even though the underlying
%    population is the same.
%
% Suppose two samples are drawn from the same distribution. There were always
% be some non-zero likelihood of obtaining a difference in the means of these
% samples, although no difference should actually be there.
%
% The following simulation shows that given a population and finite
% sample size there will always be some probability of getting a difference
% between two random samples.
%    
% Furthermore, we will show that the probability of getting large spurious
% differences increases as the sample size decreases.

% Generate a population with mean 0, and sd of 1
N = 10000;     % Choose the population size
mu = 0; % Mean fo the population
sd = 10;  % Standard deviation of the population
P  = mu + sd * randn(N,1); % Create the population, we use randn.m to sample from 
                           % a gaussian distrbution with mean mu and std sd

% Next, we will simulate three samples of different size drawn from
% the distribution P.
n = [8 64 512]; % We choose these three sample sizes.

for ii = 1:length(n)
 sa = [];sb = [];
 % We use randsample to draw samples of size n. Randsample.m is a useful
 % matlab function that accepts a vector (population) and draws samples out
 % of it. Samples can be of any size that is smaller the than population
 % size.
 sa = randsample(P,n(ii)); 
 sb = randsample(P,n(ii));
 
 % Now that we have two samples. We can compute the difference of the means
 % of these two samples. We call these samples 'empiricl' because we are
 % simulating a possible data set.
 d_empirical(ii) = mean(sa) - mean(sb);
end

% Because the three samples were drawn from the same distribution their
% differences should be 0. But with samples of small size some spurious
% difference can be observed due to chance error.
%
% Hereafter, we compute the probabilty of getting a spurious non-zero
% difference between three samples. We test how the size of the sample
% changes this probability. We will show that the smaller the sample size
% the higher the probabilty of a spurios non-zero difference. 

% We want to make a pretty plot of the results. So, we set up as structure
% of strings identifying the color to use for each distribution.
color = {'b','r','g','k'};
figure('name','Testing effects of samples size','color','w')
% We set hold to on will make sure the three plots will be plotted o the same
% axes
hold on

% Now we compute the distribution of differences between the means of two
% random samples. We compute this distribution for each sample size and
% plot the distributions one on top of eachother (with different color).
for jj = 1:length(n)
 % We drawn many pairs of two samples 
 % with our sample sizes: 8, 64 and 512
 sa = zeros(N,n(jj)); % zeros.m allows to preallocate memory
 sb = sa;
 for ii = 1:N
  sa(ii,:) = randsample(P,n(jj));
  sb(ii,:) = randsample(P,n(jj));
 end
 
 % We compute the distribution of differences between the means of sample
 % pairs. Note this is a vector operation, we subtract two vectors
 % element-wise.
 d(jj,:) = mean(sa,2) - mean(sb,2);

 % Let's compute the histogram of the distribution and trasform the
 % distribution into a probabilities.
 %
 % We choose the sampling spacing between the bins.
 binWidth = .1;
 % We make bins that cover the full range fo values in P and are sampled at
 % the resolution we decided using binWidth.
 bins = min(P):binWidth:max(P);
 
 % We now make a histogram of the distribution of these differences.
 % We scale the dstributions properly to transform into probabilites.
 [nn, xx] = hist(d(jj,:),bins);

 bar(xx,nn / (N * binWidth),color{jj});
end

% Here we do some formatting of the figure and axes. We set the max and min
% for the axis.
axis([-10 10 0 1]) 
% We set a descriptive title and labels. 
title(sprintf('Probability of getting an spurious\ndifference between two random samples\nas function of sample size.'),'fontsize',14)
ylabel('Probability of occurrence','fontsize',14)
xlabel('Difference between two random samples','fontsize',14)
% We add a legend to describe the data.
legend({'Sample size: 8','Sample size: 64', 'Sample size: 512'},'fontsize',14)
set(gca,'tickdir','out', 'box','off','fontsize',14)

% Notice that the blue distribution, the one generated by resampling with
% the smallest sample size (8) has the largest standard deviation. This
% means that with such small sample size there is a higher chance of
% observing a spurious difference. For example the probability of
% occurrence of a difference of nearly +/-0.5 should be around 0.1 (10%,
% blue bar). The probability goes down to less then 0.01 with a larger
% sample ize (65, red distribution). 

% Now lets take a look at where the original (d_empirical) differences we
% had simulated above lay in respect to each distribution.
% 
% We will plot a vertical line indicating the empirical difference. The
% height of the distributions of the matching color (i.e., metching sample
% size) at the location of the empirical difference (verticla colored line)
% is the probability that the empirical difference was otained by chance
% given that sample size.
for ii = 1:length(n)
 plot([d_empirical(ii), d_empirical(ii)],[0 1],'-','Color',color{ii},'LineWidth',5);
end

% Questions:
% If a scientific article were to report a *statistically* significant difference of 0.5 on this scale.
% How would you be able to tell whether their analyses and statistical tests were performed correctly?
% How could we relate this to a power analysis?

%% 2. Nonparametric alternatives to the t-test
%
% A t-test is a statistical method to establish the likelihood of the
% difference between the means of two samples drawn from an unknown
% distribution. The shape of this distribution is assumed to be gaussian.
%
% The null-hypothesis for this test is: The two sample come from the same
% distribution. 
%
% A t-test finds the probability that the two means were drawn from the
% same distribution. If that probability is small we can reject the null
% hypothesis and support the alternative hypothesis that the two samples
% did not come from the same distribution. 
%
% Here after we show how to use bootstrap to generate the equivalent of a
% t-test.

% We first generate the distribution for a population that is not gaussian.
%
% Generate the population using mixture of a Chi-Square distribution with 6
% degrees of freedom and a gaussin distribution centered at 0 and 1
% standard diviation. 
%
% Because this population is not Gaussian, a parametric t-test would not be
% appropriate.
N = 10000;     % population size
P = [random('chisquare',6,[N/2,1]);random('norm',0,1,[N/2,1])];

% We compute three descriptive statistics for the population, mean, median
% and standard deviation:
mu = mean(P);
me = median(P);
sd = std(P);

% Let's take a look at the population, we just generated.
% 
% We make a histogram of population and appreciate how well mean, median and
% standard deviation might describe the distribution.
figure('name','Mean, median and std for a non gaussian population','color','w')
hist(P,100);
hold on % keep the current axes on for plotting we will add more to this plot.
set(gca,'tickdir','out', 'box','off','fontsize',14)

% Here we use get.m to get the limits of the y-axis from the plot.
% We will use these values for the next plots.
y = get(gca,'yLim'); % y(2) is the maximum value on the y-axis.

% We now add mean (in red), median (green) and +/- 1 standard deviation
% around the mean.
plot([mu, mu],y*.99,'r-', 'lineWidth',2)
plot([me, me],y*.99,'g-', 'lineWidth',2)
plot([mu - sd,mu + sd],[y(2)*.5,y(2)*.5],'r-', 'lineWidth',2); % We plot the line up to the maximum on the y-axis
title('Simulated true distribution (Chisquare, 6 DOF + Gaussian)','fontsize',14)   % we add a title.

% Now let's set up a bootstrap test of the difference between the means of
% two samples that does not rely on the assumption that the popluation is
% gaussian. Our population is not.
%
% We simulate two empirical samples of size n and m drawn from our
% distribution.  As if we were to run two conditions of the same
% experiment.
%
% By drawing the two sampels from the same distribution we are effectively
% making the null hypothesis true.
% 
% We use randsample.m generate the samples.
n = 4; m = 8; 
sa = randsample(P,n);
sb = randsample(P,m);

% We compute the empirical difference between the means of the two samples.
d_empirical = mean(sa) - mean(sb);

% The first test we show is called a Randomization test (Test of the
% difference of means of two samples).
%
% So, we want to test the null hypothesis (H0) that the means of the two
% samples come from the same distribution.
%
% To do so we will first aggregate the two samples. Under the null
% hypothesis the two samples are interchanganble because they come from the
% same distribution.
sh0 = [sa; sb];

% We then resample them 10,000 times. Effectively building 10,000 new
% samples from the aggregated sample, sh0. 
N = 10000; % number of times we will resample
sa_rand = zeros(n,N);
sb_rand = zeros(m,N);
for ii = 1:N 
 sa_rand(:,ii) = randsample(sh0,n,true);
 sb_rand(:,ii) = randsample(sh0,m,true);
end

% We then compute the differences between the means of these resampled
% samples. Please notice this is a vector operation.
d = mean(sa_rand,1) - mean(sb_rand,1);

% Because our samples were drawn from the same distribution (i.e. the H0
% is true), this distribution is centered at 0. Let's take a look at it.
figure('name','Randomization test','color','w')
% Let's compute the histogram of the distribution and trasform the
% distribution into a probabilities.
%
% We choose the sampling spacing between the bins.
binWidth = .1;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(d):binWidth:max(d);

% We now make a histogram of the distribution of these differences.
% We scale the dstributions properly to transform into probabilites.
[nn, xx] = hist(d,bins);
bar(xx,nn / (N * binWidth),'k');

title('Distribution of the differences between the means of samples Sa and Sb','FontSize',14)
ylabel('Probability of occurrence','FontSize',14)
xlabel('Difference (sa - sb)','FontSize',14)
set(gca,'tickdir','out', 'box','off','FontSize',14)
hold on

% Here we use get.m to get the limits of the y-axis from the plot. We will
% use these values for the next plots.
y = get(gca,'yLim'); % y(2) is the maximum value on the y-axis.

% We can now show the empirical difference we simulated (d_empirical). The
% height of the distribution at the value of d_empirical is the probability
% that that difference comes from samples drawn from the same distribution.
plot([d_empirical,d_empirical],y*.99,'r-','lineWidth',2)

% We can obtain the probability of H0 being true by counting the number of
% randomly obtained values that are larger than the actual observed value
% and divide this by the total number of simulations that were run.
p = sum(abs(d) > abs(d_empirical))/N;
text(d_empirical,y(2)*.95,sprintf('H0 is true with %2.2f probability.',p),'FontSize',14)

% The second test we show is called a Bootstrap test (Test of the
% difference of mean of two samples).
%
% This example is identical to the previous one, with one important
% difference. During the bootstrap process samples are generated 'with
% replacement.'
%
% This means that on each new draw there will be equal likelihood to sample
% any element of the aggregated sample.
% 
% We aggregate the two samples.
sh0 = [sa; sb];

% We resample them 10,000 times with replacement.
N = 10000;
sa_rand = zeros(n,N);
sb_rand = zeros(m,N);
for ii = 1:N
 sa_rand(:,ii) = randsample(sh0,n,true);
 sb_rand(:,ii) = randsample(sh0,m,true);
end

% Please notice that matlab allows us to use the same code for these two
% tests with a single difference we set the 'replace' option in
% randomsample.m to 'true.' (True is a logical variable that matlab
% recognizes as 1).

% Now we compute the differences between the means of these resampled
% samples.
d = median(sa_rand) - median(sb_rand);

% Because our samples were drawn from the same distribution (the H0 is
% true), this distribution is centered at 0.
figure('name','Bootstrap test of the difference between means','color','w');
% Let's compute the histogram of the distribution and trasform the
% distribution into a probabilities.
%
% We choose the sampling spacing between the bins.
binWidth = .1;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(d):binWidth:max(d);

% We now make a histogram of the distribution of these differences.
% We scale the dstributions properly to transform into probabilites.
[nn, xx] = hist(d,bins);
bar(xx,nn / (N * binWidth),'k');
title('Distribution of the differences between the means of samples sa and sb','FontSize',14)
ylabel('Probability of occurrence','FontSize',14)
xlabel('Difference (sa - sb)','FontSize',14)
set(gca,'tickdir','out', 'box','off','FontSize',14)
hold on

% Here we uset get.m to get the limits of the y-axis from the plot. We will
% use these values for the next plots.
y = get(gca,'yLim'); % y(2) is the maximum value on the y-axis.

% We can now show the empirical difference we simulated (d_empirical). The
% height of the distribution at the value of d_empirical is the
% probability that that difference comes from samples drawn from the same
% distribution. 
plot([d_empirical,d_empirical],y,'r-','lineWidth',2)

% We can obtain the probability of H0 being true by counting the number of
% randomly obtained values that are larger (and smaller) than the actual observed value
% and divide this by the total number of simulations that were run.
p = sum(abs(d) > abs(d_empirical))/N;

% Here we use the fucntion sprintf.m to generate a string containing the
% result of the text. sprintf.m and fprintf.m are very useful matlab
% functions tat allow to generate and display strings of text by taking
% matlab variables as inputs.
txt = sprintf('H0 is true with %2.2f probability.',p);

% The function text.m allows to add text to the current axis.
% We plot the text at 5% lower than the maximum value on the y-axis
text(d_empirical,y(2)*0.95,txt,'FontSize',14)

% The following lines show how to set up a parametric t-test using matlab.
% This is out of scope for the present tutorial. 
[h,pp] = ttest2(sa,sb);
txt2 = sprintf('\nH0 is true with %2.2f probability (parametric t-test)\n',pp);
text(d_empirical,y(2)*0.9,txt2,'FontSize',14)

%% 3. Exploring a more complex dataset: two variables, one condition
% Hereafter, we show an example of a slightly more complex data set.
%
% So suppose we have one condition and measure not just one quantity (which
% was the subject of Lecture 1) but measure two distinct quantities. For
% example, suppose we measure both the heights and weights of male adults.
% 
% What can we do with the data? 

% Let's generate two samples rapresenting two variables drawn from gaussian
% distributions. The two variables are correlated with eachother. 

n = 6; % This is the size of our samples
% We add correaltion between the two samples.
w = .8;  
% We generate the two samples.
s1 = randn(n,1);
s2 = w .* s1 + (1 - w).*randn(n,1);


%% 4. Correlation
%
% We can quantify how correlated two variables are by using the metric
% 'correlation.'
%
% Correlation values lie in the range -1 to 1, where 1 indicates a perfect
% negative linear relationship, 0 indicates no relationship, and 1
% indicates a perfect positive linear relationship. (Note that we
% are Pearson's product-moment correlation, but there are other variants of
% correlation.) 
temp = corrcoef(s1,s2);
r = temp(1,2);

% In words, we z-score each variable (subtract off the mean, divide by the
% standard deviation) and then compute the average product of the
% variables. 
%
% Technical note: in the above formula, std should be computed
% using a version of standard deviation where we normalize by n instead of
% n - 1. Tht is why we write std(v1,1) instea of simply std(v1).

% Here we show the results in a plot
figure('name','Simulated correlation plot','color','w')
% scatter.m is a ueful matlab function, in this case it would be equivalent
% to plot(s1,s2,'bo')
scatter(s1,s2); 
ylabel('Values in variable 2','FontSize',14)
xlabel('Values in variable 1','FontSize',14)
set(gca,'tickdir','out', 'box','off','FontSize',14)

% We display the simulted and computed correlations at 80% the max y-value.
y = get(gca,'yLim');
x = get(gca,'xLim');
text(x(1)+.15,y(2)*.8,sprintf('Correlations:\nSimulated %2.2f\nEstimated %2.2f',w,r),'FontSize',14)


%% 5. Error bars
% - How do we obtain error bars on the correlation observed in a set of
%   data?
% - How do we test whether the correlation observed in a set of data is
%   significantly different from zero? 
%
% Hereafter, we compute errorbars on the correlation coefficient using
% bootstrap.
%
% We compute N bootstrap (resampling with replacement) with the same size
% of the original samples. 
%
% We then compute the correlation coefficient with each sample.
N = 1000; % number of randomizations
r_dist = zeros(N,1);
for ii = 1:N
 % Notice that we resample index and we use the same indexes to address
 % both samples. otherwise we lose the correaltion value. 
 index = randsample(1:n,n,1);
 temp = corrcoef(s1(index), s2(index));
 r_dist(ii)    = temp(2,1);
   %sum((s1(index) - mean())/std(s1(index),1) .* ...
   %                   ( - mean(s2(index))/std(s2(index),1)))/n;
end


% We copute the two-tailed 95% confidence intervals:
ci = prctile(r_dist,[2.5,97.5]);

% Now we plot the correlation we estimated with the the 95% confidence
% intervals.
figure('name','Correlation Confidence Intervals','color','w')
bar(r); % plot a bar to indiecate the etimated correlation
hold on
set(gca,'tickdir','out', 'box','off','FontSize',14)
plot([1 1],ci,'r-','lineWidth',3); % Plot the error on the estimate
title('Estimated correlation coefficient with error')
ylabel('Correlation s1, s2')
axis([0 2 0 1]); % set the limits to the current axis


%% 6. Test hypothesis
%
% Now we want to test the hypothesis that the correlation r, was not
% obtained by chance. We want to find the probability for the null
% hypothesis to be true.
%
% If the null hypothesis were true, we would be able to shuffle the
% order of each variable and obtain datasets that are equivalent to the
% original dataset. So, to obtain a p-value, we shuffle each variable,
% calculate a correlation value for each shuffled sample. We repeat this
% process a large number of times and count the number of times that
% randomly obtained correlation values are more extreme than the actual
% observed correlation value.

% We build a bootstrap distribution of samples under the null hypothess
% that s1 and s2 were not correlated.
N = 10000; % number of randomizations
r_dist = zeros(N,1);
s1_rnd = zeros(size(s1)); s2_rnd = s1_rnd;
for ii = 1:N
 % Notice the difference here. We resample s1 and s2 indepedently. This is
 % because under the null hypothesis the two sampels have no correlation.
 s1_rnd = randsample(s1,n,1);
 s2_rnd = randsample(s2,n,1);
 temp = corrcoef(s1_rnd, s2_rnd);
 r_dist(ii)    = temp(2,1);
end

% Important to remember the difference between resampling to compute error
% bars and to test hypothesis. 
%
% (A) When we want to compute error bars we resample from the two samples
% using the same indices (paired indices) during the resmapling process.
% This maintains the correlation level.
%
% (B) When we want to test the null hypothesis that two samples are not
% correlated (r=0) We resample frm the two samples using different indices.
% This effectively instatiates the null hypothesis that the two samles have
% no correlation and returns a series of correlations values that occurr
% only due to chance error.
%
% We show the distribution of correlations obtained by chance given the two
% samples.
figure('name','Test for Correlation (bootstrap)','color','w');
% Let's compute the histogram of the distribution and trasform the
% distribution into a probabilities.
%
% We choose the sampling spacing between the bins.
binWidth = .01;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(r_dist):binWidth:max(r_dist);

% We now make a histogram of the distribution of these differences.
% We scale the dstributions properly to transform into probabilites.
[nn, xx] = hist(r_dist,bins);
bar(xx,nn / (N * binWidth),'k');
title('Distribution of correlations between s1 and s2','FontSize',14)
ylabel('Probability of occurrence','FontSize',14)
xlabel('Pearson correlation coefficient','FontSize',14)
hold on
set(gca,'tickdir','out', 'box','off','FontSize',14)
y = get(gca,'yLim');

% Next we can show the empirical correlation (r). The height of the
% distribution at the value of d_empirical is the probability that that
% difference comes from samples drawn from the same distribution.
plot([r,r],y*.99,'r-','lineWidth',2)

% We can obtain the probability of H0 being true by counting the number of
% randomly obtained correlation coefficient values that are larger than the
% actual observed value and divide this by the total number of simulations
% that were run. Because value can be larger or smaller than the mean we
% coupute the absolute value of the ditribution and r values.
p = sum(abs(r) < abs(r_dist))/N;
text(0,y(2)*.95,sprintf('H0 is true with %2.2f probability.',p),'FontSize',14)


%% 7. Monte Carlo Method.
%
% - While we are on the topic of p-values on correlation values, it is
% convenient to introduce here Monte Carlo methods. Monte Carlo methods are
% a very general class of methods that make use of randomly generated data
% to test various hypotheses. 
%
% Imagine that we have obtained two samples with a correlation of r = 0.4
% for a sample size of 10. We do not know the population or have any other
% information.
%
% But we can assume that the true populations generating the two variables
% were normally distributed. We can then simulate a process in which we
% draw samples of size 10 from two independent Gaussian distributions and
% compute the correlation value observed in each sample.
% 
% Such (Monte Carlo) simulations reveal that it is actually quite likely to
% obtain a correlation of r = 0.4 for a sample size of 10, even when the
% underlying distributions are independent. Thus, we should not believe the
% claim.

% Empirical (measured) correlation. 
r_empirical = .4;

% We simulate two independent and normally distributed variables.
N = 10000; % We set the number of montecarlo simulations
n = 10;    % The sample size
s1_sim = randn(N,n); % siulated distribution of samples
s2_sim = randn(N,n); % simulated distribution of samples

% compute the correlation coefficient
for ii = 1:N
  % Note, here we are using matlab's corr.m to compute the pearon
  % correlation instead of computing it manually.
  r_dist(ii) = corr(s1_sim(ii,:)',s2_sim(ii,:)','type' ,   'Pearson' );
end

% We plot the distribution of correlation coefficients obtained this way.
figure('name','Test Correlation (Monte Carlo Method)','color','w');
% Let's compute the histogram of the distribution and trasform the
% distribution into a probabilities.
%
% We choose the sampling spacing between the bins.
binWidth = .01;
% We make bins that cover the full range fo values in P and are sampled at
% the resolution we decided using binWidth.
bins = min(r_dist):binWidth:max(r_dist);

% We now make a histogram of the distribution of these differences.
% We scale the dstributions properly to transform into probabilites.
[nn, xx] = hist(r_dist,bins);
bar(xx,nn / (N * binWidth),'k');
title('Distribution of the Montecarlo-simulated correlation coefficients','FontSize',14)
ylabel('Probability of occurrence','FontSize',14)
xlabel('Corr(s1,s2)','FontSize',14)
hold on
set(gca,'tickdir','out', 'box','off','FontSize',14)

% Get the limits of the y axis from the plot. We will use these for the
% next plots.
y = get(gca,'yLim');

% We can now show the empirical difference we simulated (d_empirical). The
% height of the distribution at the value of d_empirical is the probability
% that that difference comes from samples drawn from the same distribution.
plot([r_empirical,r_empirical],y*.99,'r-','lineWidth',2)

% We can obtain the probability of H0 being false by counting the number of
% randomly obtained correlation coefficient values that are larger than the
% empirical value and divide this by the total number of simulations that
% were run.  Because the value can be larger or smaller than the mean we
% coupute the absolute value of the ditribution and r values.
p = sum(abs(r_empirical) < abs(r_dist))/N;
text(-0.8,y(2)*.95,sprintf('H0 is false with %2.4f probability.',p),'FontSize',14)

