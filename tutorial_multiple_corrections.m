%Exercise
%In this exercise we will run through an example of correcting for multiple comparisons with both the Benjamini-Hochberg procedure and the more conservative Bonferroni correction.

%First, simulate multiple (say, 1000) t-tests comparing two samples with equal means and standard deviations, and save the p-values. Obviously, at p<0.05 we expect that ~5% of the simulations to yield a "statistically significant" result (of rejecting the NULL hypothesis that the samples come from distributions with equal means).
%Second, once you have the simulated p-values, apply both methods to address the multiple comparisons problem.
%Third, set the sample 1 and sample 2 means to be 1 and 2 respectively, and re-run the exercise. What do you notice? What if you make the difference between means even greater?
clear


%%
%simulate t tests
num_tests = 1000
t_test = zeros(num_tests,1);
alpha = .05
sample_size = 75

for i = 1:num_tests
mean_s1 = 1;
stdev_s1 = 1;
mean_s2 = 2;
stdev_s2 = 1;

dist_s1 = normrnd(mean_s1,stdev_s1,sample_size,1);
dist_s2 = normrnd(mean_s2,stdev_s2,sample_size,1);
[h,p] = ttest2(dist_s1,dist_s2)
t_test(i,1) = p

end

%%
%count up number of significant p values without mutliple correction
p_count = 0
for i = 1:num_tests
    if t_test(i,1)<alpha
        p_count = p_count + 1
    end
    i = i+1
end



%%
%apply Bonferoni correction and count significant p values
p_count_bonferroni = 0
for i = 1:num_tests
    if t_test(i,1)<alpha/num_tests
        p_count_bonferroni = p_count_bonferroni + 1
    end
    i = i+1
end



%%
%apply Benjamini–Hochberg procedure
fdr = .05
p_count_bh = 0;
bh_dist = sort(t_test,'descend')
bh_found = 0


%%
%calculate critical value
for i = 1:num_tests
    if bh_found==0
        critical_value = (i/num_tests)*fdr
        if bh_dist(i,1)<critical_value
         bh_p_cutoff = bh_dist(i,1)
            bh_found = 1
        end
    end
end

%calcualte number of significant p values for bh
for i = 1:num_tests
    if t_test(i,1)<bh_p_cutoff
        p_count_bh = p_count_bh+ 1
    end
end

%ANSWER: When I ran this assignment with all values the same except for both means
%were 1, no correction resulted in 60 significant values, bonferroni
%resulted in 0, and BH resulted in 57. I ran the code a few times and there
%was some variability with this count. 

%I re-ran the code with mean_s2 = 2, and with no corrections the p count
%was 1,000, BH was 992, and bonferroni was 964


