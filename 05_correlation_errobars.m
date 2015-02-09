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

n = 20; % This is the size of our samples
% We set a correlation level between the two samples.
w = .6;  

% We generate the two samples.
s1 = randn(n,1);
s2 = w .* s1 + (1 - w).*randn(n,1);


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
r = sum((s1 - mean(s1))/std(s1,1) .* (s2 - mean(s2)/std(s2,1)))/n;

% In words, we z-score each variable (subtract off the mean, divide by the
% standard deviation) and then compute the average product of the
% variables. 
%
% Technical note: in the above formula, std should be computed
% using a version of standard deviation where we normalize by n instead of
% n - 1. Tht is why we write std(v1,1) instea of simply std(v1).
% This is the version of the std that the Book Uses.

% Here we show the results in a plot
figure('name','Two correlated variables')
% scatter.m is a ueful matlab function, in this case it would be equivalent
% to plot(s1,s2,'bo')
scatter(s1,s2); 
ylabel('Values in variable 2','fontsize',14)
xlabel('Values in variable 1','fontsize',14)

% We display the simulted and computed correlations at 80% the max y-value.
y = get(gca,'yLim');
x = get(gca,'xLim');
% Please notice here we introduce a new function that allows to write text
% on a plot. try doc text.m 
text(x(1)+.15,y(2)*.8, ...
    sprintf('Correlations:\nSimulated %2.2f\nEstimated %2.2f',w,r),'fontsize',14)
set(gca,'tickdir','out','ytick',[-1 0 1],'xtick',[-1 0 1])

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
k = 2000; % number of randomizations
r_dist = zeros(k,1);
for ii = 1:k
 % Notice that we resample index and we use the same indexes to address
 % both samples. otherwise we lose the correaltion value. 
 index = randsample(1:n,n,1);
 r_dist(ii)    = sum((s1(index) - mean(s1(index)))/std(s1(index),1) .* ...
                     (s2(index) - mean(s2(index))/std(s2(index),1)))/n;
end

% We compute the two-tailed 95% confidence intervals:
ci = prctile(r_dist,[2.5,97.5]);

% Now we plot the correlation we estimated with the the 95% confidence
% intervals.
figure('name','Errorbar on the correlation coefficient','color','w');
plot([1 1],ci,'r-','lineWidth',3); % Plot the error on the estimate
hold on
plot(1,r,'ko','markeredgecolor','w','markerfacecolor','k',...
         'markersize',14,'linewidth',2); % plot a bar to indicate the etimated correlation
ylabel('Correlation s1, s2')
axis([0 2 0 1]); % set the limits to the current axis
set(gca,...
    'tickdir','out','ytick',[-1 -.5 0 .5 1], ...
    'xtick',[-1 0 1],'ylim',[-1 1], ...
    'box','off','fontsize',14)
