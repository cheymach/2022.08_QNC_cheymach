%generate a sample from a normal distribution with mu=10 and 
% standard deviation =2, for n=5, 10, 20, 40, 80, 160, 1000 at a 95% confidence level
%generate a sample from this distribution
clear
N_values = [5 10 20 40 80 160 1000]

for i = 1:7
    n = N_values(1,i)
    mu = 10
    sigma = 2
    sample_distribution = normrnd(mu,sigma,20,1)

%%

%Calculate the confidence interval using method 1
    sample_mean = mean(sample_distribution)
    sample_sem = std(sample_distribution)/sqrt(n)
    upper_bound_ci1 = sample_mean + sample_sem*1.96
    lower_bound_ci1 = sample_mean - sample_sem*1.96
    ci1 = [lower_bound_ci1 upper_bound_ci1]
%%

%Calculate the confidence interval using method 2, students T
%uses matlab's tinv function, using code documented here: https://www.mathworks.com/help/stats/tinv.html
    nu = n-1
    lower_bound_probability = .025
    upper_bound_probability = .975
    ci2_bound = tinv([lower_bound_probability upper_bound_probability], nu)
    upper_bound_ci2 = sample_mean + sample_sem*ci2_bound(1,2)
    lower_bound_ci2 = sample_mean - sample_sem*ci2_bound(1,2)
    ci2 = [lower_bound_ci2 upper_bound_ci2]

%%
%Calculate the confidence interval using method 3, bootstrapping
    ci3 = bootci(200,@mean,sample_distribution)


%% calculate the credible interval using method 4, bayesian
%ultimately calculated the same way as method 1
    ci4 = [lower_bound_ci1 upper_bound_ci1]

%create an array comparing confidence intervals for various values of n
%first column is method number, second column is confidence interval lower
%bound, third column
%is confidence interval upper bound
    summary_ci = zeros(4,3)
    methods = [1 2 3 4]
    methods_label = methods'
    summary_ci(:,1) = methods_label
    summary_ci(1,2:3) = ci1
    summary_ci(2,2:3) = ci2
    summary_ci(3,2:3) = ci3
    summary_ci(4,2:3) = ci4

    if i == 1
        concat_summary = summary_ci
    else
        concat_summary = cat(3,concat_summary,summary_ci)

    end
end
%%
%Print nicely
disp("The confidence intervals by method for n=")
disp(N_values(1,1))
disp(concat_summary(:,:,1))

disp("The confidence intervals by method for n=")
disp(N_values(1,2))
disp(concat_summary(:,:,2))


disp("The confidence intervals by method for n=")
disp(N_values(1,3))
disp(concat_summary(:,:,3))


disp("The confidence intervals by method for n=")
disp(N_values(1,4))
disp(concat_summary(:,:,4))

disp("The confidence intervals by method for n=")
disp(N_values(1,5))
disp(concat_summary(:,:,5))

disp("The confidence intervals by method for n=")
disp(N_values(1,6))
disp(concat_summary(:,:,6))

disp("The confidence intervals by method for n=")
disp(N_values(1,7))
disp(concat_summary(:,:,7))