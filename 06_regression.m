%% 1. Regression. Exploring more complex datasets: two variables, one condition
% Hereafter, we show an example of a slightly more complex data set.
%
% So suppose we have two variables and. For example, suppose we measure
% both the heights and weights of male adults.
% 
% How can we express the relationship between the two variables by
% computing the correlation coefficient. We can use linear regression to
% estimate one variable from the other. Linear regression will go beyond
% correlation because it will allow for predicting values that are not in
% the sample. Prediction is one of the fundamental properties of models. If
% we have a model we are going beyond the data, we can make predictions of
% data we have never collected.
%
% To start with our example we will first simulate a population of linearly
% related data sets, x and y.
%
% Let's generate two samples rapresenting two variables drawn from gaussian
% distributions. The two variables are correlated with eachother. 
n  = 200;   % This is the size of our sample, meaning the number of x values we sample.
sd = 250; % We add noise alittle bit more variable that default
close all

% Hereafter we will introduce the concept of population. A populations is
% the theoretical pool of data out of which our specific sample is drawn.
% In general a population cannot be measured. But we have a sample and we
% hope for the sample to represent the true population. Some samples are
% better some worse at representng the true population.
%
% We will first simulate a population. Using a "true" model of our choice;
% a straight line.
x = repmat(1:n,n,1)';
x = x + sd*randn(size(x));% We add some new noise to X in the population
a = 1.5; % The slope and 
b = 2;  % the intercept of the true popolation

% Polyval is a function that allows computing polinomials. Here after we
% use the function to make a line.
polynomialOrder  = 1; % This is a parameter for the MatLab function polyval.m 
y = polyval([a b],x,polynomialOrder);
y = y + sd*randn(size(y)); % We add some new noise to Y in the population

% Show the simulated true population.
h = figure('Name','Population and sample data','color','w')
plot(x,y,'ro','MarkerFaceColor','r')
xlabel('x','fontsize',14);
ylabel('y','fontsize',14)
set(gca, ...
    'fontsize',14, ...
    'box','off','tickdir','out')

% Now we show the true model representing the population.
% This means the model used to build the population without added noise.
true_m = polyval([a b],x,polynomialOrder);
hold on
plot(x,true_m,'k-', 'LineWidth',4)
axis equal

%% 2. We can now 'sample' data out of this theoretical population
% We will now sample from the population.
sampleSize      = .25*n; % This is the total number of measurements in x
numObservations = 2; % Observations in y at each x
ind1 = randsample(1:size(x,1),sampleSize);
ind2 = randsample(1:size(x,2),numObservations);
x_sample = x(ind1,ind2);
y_sample = y(ind1,ind2);

% Show the sample.
plot(x_sample,y_sample,'bo','MarkerFaceColor','w','MarkerSize',12)

%% 3. Fit a regression model to the sample.
%
% Next we will fit a regression line to the data in the sample.
% The model, a straight-line, will estimate the slope and intercept
% for the population. TO do so we will find the least-square
ab = polyfit(x_sample,y_sample,polynomialOrder);

% Evaluate the estimated model.
y_hat = polyval(ab,x_sample);

% Show the model fit to the sample.
plot(x_sample,y_hat,'b-','LineWidth',4)
axis equal

%% 4. Correlation and coefficient of determination
%
% We can quantify how correlated two variables x andy are by using the metric
% 'correlation.'
%
% Correlation values lie in the range -1 to 1, where -1 indicates a perfect
% negative linear relationship, 0 indicates no relationship, and 1
% indicates a perfect positive linear relationship. (Note that we are
% using Pearson's product-moment correlation, but there are other variants
% of correlation.)
tmp = corrcoef(x_sample,y_sample); % MatLab has a formula for correlations
r   = tmp(1,2);

% We can now compute the coefficient of determination (R2) of the fitted
% model to the sample data. This number shows the percent of variance in
% the data explained by the model.
R2 = 100 * (1 - (sum((y_sample(:) - y_hat(:)).^2) / sum((y_sample(:) - mean(y_sample(:))).^2)));

%% 5. Next we will compute the residuals of the regression model. 
% If the regression model (a line) is a good model of the population the
% residuals should be centered at 0 and normally distributed.
res = y_hat - y_sample;

figure('name','residuals','color','w')
[y,x] = hist(res(:));
plot(x,y,'k-','linewidth',2)
hold on
plot([0 0],[0 max(y)],'k--')
set(gca, ...
    'fontsize',14, ...
    'box','off','tickdir','out')
ylabel('Frequency','fontsize',14)
xlabel('Residuals','fontsize',14)
text(min(min([x,y]))*.99, ...
     max(max([x y]))*.99, ...
    sprintf('Coefficient of determination (R2): %2.3f\nCorrelation coefficient: %2.3f',R2,r), ...
    'fontsize',14)

