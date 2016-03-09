%% 1. Correlation. Exploring more complex datasets: two variables, one condition
% Hereafter, we show an example of a slightly more complex data set.
%
% So suppose we have one condition and measure not just one quantity (which
% was the subject of Lecture 1) but measure two distinct quantities. For
% example, suppose we measure both the heights and weights of male adults.
% 
% What can we do with the data? 

% Let's generate two samples rapresenting two variables drawn from gaussian
% distributions. The two variables are correlated with eachother. 

n = 10; % This is the size of our samples
% We set a correlation level between the two samples.
p = .75;  

% We generate the two samples.
u = randn(1,n);
v = randn(1,n);
m.s1  = 5.5; % mean deviation of the first sample
m.s2  = 6; % mean deviation of the second sample
sd.s1 = 1; % standard deviation of the first sample
sd.s2 = 1.2; % standard deviation of the second sample
s1 = sd.s1 * u + m.s1;
s2 = sd.s2 * (p * u + sqrt(1 - p^2) * v) + m.s2;

%% 2. Correlation formula
%
% We can quantify how correlated two variables are by using the metric
% 'correlation.'
%
% Correlation values lie in the range -1 to 1, where -1 indicates a perfect
% negative linear relationship, 0 indicates no relationship, and 1
% indicates a perfect positive linear relationship. (Note that we are
% using Pearson's product-moment correlation, but there are other variants
% of correlation.)
tmp = corrcoef(s1,s2); % MatLab has a formula for correlations
r   = tmp(1,2);

% In words, we z-score each variable (subtract off the mean, divide by the
% standard deviation) and then compute the average product of the
% variables. 
%
% Technical note: in the above formula, std should be computed
% using a version of standard deviation where we normalize by n instead of
% n - 1. Tht is why we write std(v1,1) instea of simply std(v1).
% This is the version of the std that the Book Uses.

% Here we show the results in a plot
figure('name','Two correlated variables','color','w')
% scatter.m is a ueful matlab function, in this case it would be equivalent
% to plot(s1,s2,'bo')
subplot(1,2,1)
scatter(s1,s2); 
ylabel('Height of sons (feet)','fontsize',14)
xlabel('Height pf fathers (feet)','fontsize',14)

% We display the simulted and computed correlations at 80% the max y-value.
y = get(gca,'yLim');
% Please notice here we introduce a new function that allows to write text
% on a plot. try doc text.m 
text(0,y(2)-.2, ...
    sprintf('Simulated %2.2f\nEstimated %2.2f\n',p,r),'fontsize',14)
set(gca,'ylim', [2 8],'xlim', [2 8], 'tickdir','out', ...
        'ytick',[2 4 6 8],'xtick',[2 4 6 8])
axis square

%% 3. Error bars on the correlation coefficient using bootstrap.
%
% - How do we obtain error bars on the correlation observed in a set of
% data?
%
% - How do we test whether the correlation observed in a set of data is
% significantly different from zero? 
%
% Hereafter, we compute errorbars on the correlation coefficient using
% bootstrap. Bootstrap is a non-parametric statistical method that allows
% for computing tests and error without assuming a shape for the underlying
% distribution in the data. We will talk more abot Bootstrap later on.
%
% Hereafter we use randsample.m to bootstrap. But MatLab also has a
% bootstrap function called bootstrp.m, please read abotu it. We will use
% it later one. 
% 
% Th logic to compute the error bars of a correlation coefficient in via
% bootstrap is simple: 1. We sample the data with replacement and compute
% compute the correlation coefficient.
% 
% 2. We repated the operation k times where k is a large number, at least
% 1000 better if 10,000.
%
% 3. We compute mean and standard deviation of the distribution of 10,000
% correlation coefficients we comptued with the data.
%
% The computations will give us an estimate of how well the data sample we
% were given constrains the correlation value we measured. If the data are
% nto very noisy they sampling and computing the correlationw will compute
% very similar values for the correlation. If the data are more noisy every
% time we draw a new set of elements out of the sample the correlation
% value will be different.
% 
% We then compute the correlation coefficient with each sample.
k = 10000; % number of randomizations, normally numbers must be more 1000
r_dist = zeros(k,1);
for ii = 1:k
 % Notice that we resample index and we use the same indexes to address
 % both samples. otherwise we lose the correaltion value. 
 index      = randsample(1:n,n,1);
 tmp        = corrcoef(s1(index),s2(index));
 r_dist(ii) = tmp(2,1);
end

% We copute the two-tailed 95% confidence intervals:
ci = prctile(r_dist,[2.5,97.5]);

% Now we plot the correlation we estimated with the the 95% confidence
% intervals.
subplot(1,2,2)
plot([1 1],ci,'r-','lineWidth',3); % Plot the error on the estimate
hold on
plot(.5,p,'k^','markeredgecolor','g','markerfacecolor','k',...
         'markersize',14,'linewidth',2); % plot the true correlation
plot(1,r,'ko','markeredgecolor','w','markerfacecolor','k',...
         'markersize',14,'linewidth',2); % plot a bar to indicate the etimated correlation
ylabel('Correlation(s1, s2)')
axis([0 2 r-.1 r+.1 ]); % set the limits to the current axis
set(gca, ...
    'tickdir','out','ytick',[r-.1 r r+.1], ...
    'xtick',[0 1 4],'ylim', [r-.1 r+.1], ...
    'box','off','fontsize',14)
axis square
h = legend({'Uncertainty on estimate', ...
        'Real correlation', ...
        'Estimated from sample'}, ...
        'box','off','Location','NorthEastOutside');
