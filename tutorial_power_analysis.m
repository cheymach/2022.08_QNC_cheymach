
%obtain null distribution
% computing correlation coefficient on simulated data generated from Poisson
%(LC spiking) and Gaussian (pupil data) distributions

%simulate data from Poisson distribution, with a random mean and stdev and
%an n of 1,000
clear

%create an array to save my corr coefficients
corr_coeff_dist_random = zeros(1000,1)


for i = 1:1000

    LC_random_lambda = randi(5)
    LC_random_dist = poissrnd(LC_random_lambda,1000,1)


%simulate data from Gaussian distribution, with a random mean and stdev and
%an n of 1,000

    pupil_random_mean = randi(8)
    pupil_random_sigma = randi(8)
    pupil_random_dist = normrnd(pupil_random_mean,pupil_random_sigma,1000,1)

%compute correlation coefficients between the two distributions
    corr_coeff_random = corrcoef(pupil_random_dist,LC_random_dist)
%save corrcoef into an array

    corr_coeff_dist_random(i,1) = corr_coeff_random(1,2)
end


%%
%now that i Have generated a null didstribution of correlation
%coefficients, proceed to plot n as a function of effect sizes for 80%
%power
power = .8
%true_mean is the mean of the null distribution, stdev is the std of the
%null distribution
true_mean = mean(corr_coeff_dist_random)
stdev = std(corr_coeff_dist_random)
effect_sizes_neg = -.9:0.1:-.1
effect_sizes_pos = .1:.1:.9
effect_sizes = cat(2,effect_sizes_neg,effect_sizes_pos)
%create array to save the effect size and the associated required n
summary_required_n = zeros(18,2)
summary_required_n(:,1) = effect_sizes


%%
for i = 1:18
    sample_mean = true_mean + effect_sizes(1,i)
%for t test, structure is samplesizepwr('t' [sample_mean
%stdev],true_mean,pwr)
%to test effect sizes, I will keep the true_mean constant and vary the
%sample_mean
    n_required = sampsizepwr('t',[sample_mean stdev],true_mean,power)
    summary_required_n(i,2) = n_required
end


%%
%plot my results
bar(summary_required_n(:,1), summary_required_n(:,2))
xlabel('effect size')
ylabel('n required')

