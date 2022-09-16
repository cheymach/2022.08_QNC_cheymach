
%Sorry this part of my homework is late!

clear

%Exercise 1
%Assume that there are 10 quanta available in a nerve terminal, and for a 
% given release event each is released with a probability of 0.2. 
% For one such event, what is the probability that 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, or 10 quanta will be released?

p_quanta = 0.2;
n_quanta = 10;       % number of "trials" per "experiment"
NumExperiments_quanta = 1000;     % number of "experiments"
outcomes_quanta = binornd(n_quanta,p_quanta,NumExperiments_quanta,1);


%%ANSWER depicted in the histogram below

% Make histogram of all possible outcomes. We want bins centered on whole
% numbers so we offset the edges
edges = -0.5:10.5;
counts = histcounts(outcomes_quanta, edges);

% Show a bar plot of the simulated bionimal distribution
clf;
xs = edges(1:end-1)+diff(edges)/2;
bar(xs, counts);
title(sprintf('Histogram of binomial distribution, n=%d, p=%.2f, N=%d simulations', ...
   n,p,Num_experiments_quanta));
xlabel(sprintf('Number of successes in %d tries', n));
ylabel('Count');




%%
%EXERCISE 2

clear

n = 14
xs = [0:14]
likelihood_by_probability_firing = zeros(10,1)
prob_firing = 0.1:0.1:1

for i = 1:10
    true_binomial_distribution = binopdf(xs,n,prob_firing(1,i));
    likelihood_By_probability_firing = true_binomial_distribution(1,9)
    disp(prob_firing(1,i))
    disp(likelihood_By_probability_firing)
end

%%ANSWER: 1.641514960799999e-04 for 1 experiment at .1 probabilty
% ANSWER: looks like highest was at a probability of .6, with  0.2066
%probability of getting 8 out of 14 firing
%The other parts of this question are also answered by running this code


%%
%Exercise 3
clear
%assuming 5 quanta released in both trials out of 14

n = 14
xs = [0:n]
prob_firing = 0:0.1:1
quanta = 5

for i = 1:11
    true_binomial_distribution = binopdf(xs,n,prob_firing(1,i));
    likelihood_by_probability_firing = true_binomial_distribution(1,(quanta + 1))
    %calculate sum of probabilities for 2 trials
    product_probabilities = (likelihood_by_probability_firing)^2
    log_likelihood = 2*log10(likelihood_by_probability_firing)
    disp("Release probability")
    disp(prob_firing(1,i)) 
    disp ("likelihood of 5 quanta firing")
    disp(likelihood_by_probability_firing)
    disp("product probabilties")
    disp(product_probabilities)
    disp("log likelihood")
    disp(log_likelihood)
end


%ANSWER: my max value was for a release probability of .4
%don't totally understand what is meant by increasing resolution in this
%context

%%
%Exercise 4
%calculate the weighted average of number of firing, then calculate p_hat
clear

probability_estimate = (0/14*0 + 1/14*0 + 2/14*3 + 4/14*10 + 5/14*19 + 6/14*26 + 7/14*16 + 8/14*16 + 9/14*5 + 10/14*5)/10
p_hat = probability_estimate/14 
%ANSWER: most likely release probability estimate is .3224


%%
%EXERCISE 5
%new phat is 7/14 (.5)
clear

n = 14
xs = [0:14]
prob_firing = 0.3
true_binomial_distribution = binopdf(xs,n,prob_firing);
disp(true_binomial_distribution(1,8))

%The probability is .0618, so using a p value threshold of .05 you cannot
%reject your null hypothesis that there is no effect of temperature on
%release probability
%%

