%% Bootstrap: What is your confidence, given your data?
%
% - This tutorial will show how to implement a bootstrap method to
%   establish the reliability of a statistics given data samples of different
%   size.
% 
% Pestilli, Franco K310 Spring 2016 Indiana University Bloomington
close all

%% (0) Simulating data.
% We create two small data sets. Imagine measuring the height of two groups of
% students in two classes. One class has 200 students the other 12.

height.mean = 5.6; % feet
height.sd   = 2;

nStudents.class1 = 200;
nStudents.class2 = 12;

% Now imagine if we were to estimate the mean height of the students at IUB
% by measuring only from class 1 or class 2.
measuredHeight.class1 = height.mean + height.sd*randn(nStudents.class1,1); 
measuredHeight.class2 = height.mean + height.sd*randn(nStudents.class2,1); 
measuredHeight.meanClass1 = mean(measuredHeight.class1);
measuredHeight.meanClass2 = mean(measuredHeight.class2);

figure('name','Measured students heights','color','w')
plot([1, 2],[measuredHeight.meanClass1, measuredHeight.meanClass2],'ro')
ylabel('Mean students height (feet)')
xlabel('Class')
set(gca,'box','off','tickdir','out','xlim',[0.5 2.5],'xtick',[1 2])

% The two estimates are somehow different.
% Which one is more precise and why?

% Likely the sample with lower numerosity (class 2) is less precise. Or
% better the estimate of the mean from that sample is less reliable.
% We can demonstrate that this is tru by repeating the above process.
measuredHeight.meanClass1(2) = mean(height.mean + height.sd*randn(nStudents.class1,1));
measuredHeight.meanClass2(2) = mean(height.mean + height.sd*randn(nStudents.class2,1));

% We plot the new means estimated from each sample
hold on
plot([1,2],[measuredHeight.meanClass1(2), measuredHeight.meanClass2(2)],'ko')

% OK the two estimates from class 1 (with more students) are more similar
% than those from Class 2 (with less students). So perhaps class 2 having
% less students have more reliable estimates.

% Let's test this hypothesis competelly. To do so we can repeate the above
% process not just twice but many times. Let's say 100 times.
nRepeats = 100;
for ii = 1:nRepeats
measuredHeight.meanClass1(ii) = mean(height.mean + height.sd*randn(nStudents.class1,1));
measuredHeight.meanClass2(ii) = mean(height.mean + height.sd*randn(nStudents.class2,1));
end

% We plot the new results:
% Which one of the two samples generates more variability around the true 
% mean? Why?
plot(ones(size(measuredHeight.meanClass1)),measuredHeight.meanClass1,'ko')
plot(2*ones(size(measuredHeight.meanClass2)),measuredHeight.meanClass2,'ro')

% Bootstrap allows to ask similar questions for cases where you have 
% no ways to repeate the experiment as we did about in our simulations.
%
% Imagine if you had only access to Class 2 and wanted to know the mean
% height of the students at the university. This is all you have.
nStudents.class2      = 12;
measuredHeight.class2 = height.mean + height.sd*randn(nStudents.class2,1); 
estimatedMean         = mean(measuredHeight.class2);

% Now you can ask the question: "OK I have this mean, how much should 
% I trust it?" The levelof trust really depends (among many things) 
% on how many measurements we made. In this case 12. So how reliable is a
% mean estimated on 12 measurements represent?
%
% We can ask this question empirically. We can use our data to estimate the
% reliability for the estimate of the mean. To do that we are going to
% resample from our sample with replacement.
resampeldData = randsample(measuredHeight.class2,nStudents.class2,true);

% OK how much did our estimate of the mean students' height change by
% simply doing that?
estimatedMean2 = mean(resampeldData);

% What if we were to do this again?
resampeldData  = randsample(measuredHeight.class2,nStudents.class2,true);
estimatedMean3 = mean(resampeldData);

% OK we can go on and notice that everytime we resampel we have a slgihtly
% different estimate of the mean. To really establish how reliable all
% these estimates are we really need to repeat this process many many
% times. By repeating the process many times let's say 10,000 times we will
% get a sense of how variable our estimated mean is given the single data
% set we collected. In other words we would estimate how well those 12
% measurements can constrain the estimate of the mean. How many other means
% could I have gotten if I only had access to these 12 measurements.
%
% Looking at the spread in means given this single data gives us a good
% idea of how reliable our measured mean might be. So let's do it. 

% Below we resample from the original sample 10,000 times and we comptue the
% mean everytime, we collect all the means and plot their distribution. A
% wide distrobution tells us that the 12 measurements do not constrain the
% mean very well, they ar enot very reliable.
numBoots = 10000;
for ib = 1:numBoots
estimatedMeans(ib) = mean(randsample(measuredHeight.class2,nStudents.class2,true));    
end

% We make a histogram to show the spread.
figure('name','Reliability of the estimates of the mean height','color','w')
[y,x] = hist(estimatedMeans,60);
plot(x,y,'r-')
set(gca,'box','off','tickdir','out','xlim',[3 7],'xtick',[3 4 5 6 7])
ylabel('Frequency')
xlabel('Mean student height (feet)')
hold on

% Next let's compare this distribution with the equivalent distribution
% that we would have obtained with a sample of larger size, like Class 1,
% which has 200 students.
for ib = 1:numBoots
estimatedMeans(ib) = mean(randsample(measuredHeight.class1,nStudents.class1,true));    
end
[y,x] = hist(estimatedMeans,x);
plot(x,y,'k-')

% Remember what was the real height:
plot([height.mean height.mean],[0 160],'ko--')

% You can see now that the red distribution under estimates the mean and
% covers a wide range of potential means. This is becuase the original
% samples was small - only 12 measurements. Collecting more measurements
% generates a better estimate yet not perfect black curve.
%
% This is an example of how bootstrapping your data can be used to get a 
% sense abotu the strenght of evidence about your claim. Where you claim 
% is the estimate from a statistics comptued on sample of data collected.

% Matlab comes with a function that computes the Confidence Intervals usign
% bootstrap. The function computes improved methods for estimating
% confidence on an estimate.
ci = bootci(numBoots,{@mean,measuredHeight.class2},'type','bca');
plot(ci, [100 100],'rs--')

ci = bootci(numBoots,{@mean,measuredHeight.class1},'type','bca');
plot(ci, [150 150],'ks--')

% End
