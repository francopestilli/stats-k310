%% Make a histogram of each variable
% Plot them one against eachother to show correlation
% Compute correlation by hand, explain STd(x,1)
% Show a coupel of weird distributions for which the correaltion does nto
% plot well, like in the book, the half moon
% 

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

n = 100; % This is the size of our samples
% We set a correaltion level between the two samples.
w = .8;  
% We generate the two samples.
s1 = randn(n,1);
s2 = w .* s1 + (1 - w).*randn(n,1);


%% 4. Correlation
%
% We can quantify how correlated two variables are by using the metric
% 'correlation.'
%
% Correlation values lie in the range ?1 to 1, where ?1 indicates a perfect
% negative linear relationship, 0 indicates no relationship, and 1
% indicates a perfect positive linear relationship. (Note that we
% arPearson's product-moment correlation, but there are other variants of
% correlation.) 
r = sum((s1 - mean(s1))/std(s1,1) .* (s2 - mean(s2)/std(s2,1)))/n;

% In words, we z-score each variable (subtract off the mean, divide by the
% standard deviation) and then compute the average product of the
% variables. 
%
% Technical note: in the above formula, std should be computed
% using a version of standard deviation where we normalize by n instead of
% n - 1. Tht is why we write std(v1,1) instea of simply std(v1).

% Here we show the results in a plot
figure
% scatter.m is a ueful matlab function, in this case it would be equivalent
% to plot(s1,s2,'bo')
scatter(s1,s2); 
ylabel('Values in variable 2')
xlabel('Values in variable 1')

% We display the simulted and computed correlations at 80% the max y-value.
y = get(gca,'yLim');
x = get(gca,'xLim');
text(x(1)+.15,y(2)*.8,sprintf('Correlations:\nSimulated %2.2f\nEstimated %2.2f',w,r))


%% 5. Error bars
%
% - How do we obtain error bars on the correlation observed in a set of
% data?
%
% - How do we test whether the correlation observed in a set of data is
% significantly different from zero? 
%

% Hereafter, we compute errorbars on the correlation coefficient using
% bootstrap.
%
% We compute k bootstrap (resampling with replacement) with the same size
% of the original samples. 
%
% We then compute the correlation coefficient with each sample.
k = 1000; % number of randomizations
r_dist = zeros(k,1);
for ii = 1:k
 % Notice that we resample index and we use the same indexes to address
 % both samples. otherwise we lose the correaltion value. 
 index = randsample(1:n,n,1);
 r_dist(ii)    = sum((s1(index) - mean(s1(index)))/std(s1(index),1) .* ...
                     (s2(index) - mean(s2(index))/std(s2(index),1)))/n;
end

% We copute the two-tailed 95% confidence intervals:
ci = prctile(r_dist,[2.5,97.5]);

% Now we plot the correlation we estimated with the the 95% confidence
% intervals.
figure;
bar(r); % plot a bar to indiecate the etimated correlation
hold on
plot([1 1],ci,'r-','lineWidth',3); % Plot the error on the estimate
title('Estimated correlation coefficient with error')
ylabel('Correlation s1, s2')
axis([0 2 0 1]); % set the limits to the current axis

