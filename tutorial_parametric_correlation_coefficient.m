%Parametric correlation coefficients
%Exercises: 
clear
wing_tail = [10.4 10.8 11.1 10.2 10.3 10.2 10.7 10.5 10.8 11.2 10.6 11.4; 7.4 7.6 7.9 7.2 7.4 7.1 7.4 7.2 7.8 7.7 7.8 8.3]
%% 

%Exercise 1: plot
scatter(wing_tail(1,:),wing_tail(2,:))
%ANSWER: Yes, they appear related
%%
%Exercise 2. Calculate using equations
%used a sample r equation from internet because did nto understand why part of the
%summation was "i-1" in the tutorial
wing_bar = mean(wing_tail(1,:))
tail_bar = mean(wing_tail(2,:))
wing_diff_sq_sum = 0
tail_diff_sq_sum = 0
wing_tail_product_sum = 0
for i = 1:12
    wing_diff = wing_tail(1,i) - wing_bar
    tail_diff = wing_tail(2,i) - tail_bar
    wing_tail_product = wing_diff*tail_diff
    wing_tail_product_sum = wing_tail_product_sum + wing_diff*tail_diff
    wing_diff_sq_sum = wing_diff_sq_sum + wing_diff*wing_diff
    tail_diff_sq_sum = tail_diff_sq_sum + tail_diff*tail_diff
end

r_wing_tail = wing_tail_product_sum/(sqrt(wing_diff_sq_sum)*sqrt(tail_diff_sq_sum))

%my calculated r was .8704
%%
%use matlab's method of calculating
r_matlab = corrcoef(wing_tail(1,:),wing_tail(2,:))
%yes, this was also .8704 -- the same as my calculated value

%%
%Exercise 3: standard error 
n = 12
stderror_r = sqrt((1-r_wing_tail*r_wing_tail)/(n-2))

%Confidence intervals
%Fisher's z transformation
p_ci = .975
z_wing_tail = .5*log((1+r_wing_tail)/(1-r_wing_tail))
stdev_z = sqrt(1/(n-3))
%compute ci in z space
ci_upper = z_wing_tail + norminv(p_ci)*stdev_z
ci_lower = z_wing_tail - norminv(p_ci)*stdev_z
%transform back
r_ci_upper = (exp(1)^(2*ci_upper) - 1)/(exp(1)^(2*ci_upper) + 1)
r_ci_lower = (exp(1)^(2*ci_lower) - 1)/(exp(1)^(2*ci_lower) + 1)

%Exercise 4. ANSWER: Yes, would be considered significant because 95% confidence interval does
%not overlap with 0
%Exercise 5: We cannot reject the null that the r value is truly .75
%because that is within the 95% confidence interval
%%
%exercise 6: 
%don't understand how to get the power AND n out of sampsize pwr because
%each equation depends ont he other value... ...also was unsure which test to use for sampsizepwr
%becuase we don't have a stdev and r is not normally distirbuted
p0 = 0
p1 = .5
beta = .2
pwr = 1 - beta
n_r= sampsizepwr('t',[p1 stderror_r],p0,pwr)




