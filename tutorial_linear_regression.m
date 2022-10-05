%linear regression
age_WingLength = [3 4 5 6 7 8 9 11 12 14 15 16 17; 1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0]
age = age_WingLength(1,:)
WingLength = age_WingLength(2,:)
%plot relationship
scatter(age_WingLength(1,:),age_WingLength(2,:))
%%
% create regression model. This also outputs the p value, SE, and tStat,
% and R squared
mdl_age_Wing = fitlm(age,WingLength)
%%
%plot the regression line with confidence bounds
plot(mdl_age_Wing)
%%
%calculate Pearson's r
r_age_Wing = corrcoef(age,WingLength)

%%add some noise
noise_age_1 = randi (10)
noise_age_2 = randi(10)
noise_age_3 = randi(10)

noise_wing_1 = randi(6)
noise_wing_2 = randi(6)
noise_wing_3 = randi(6)

noise_age_wing = [3 4 5 6 7 8 9 11 12 14 15 16 17 noise_age_1 noise_age_2 noise_age_3; 1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0 noise_wing_1 noise_wing_2 noise_wing_3]
noise_age = noise_age_wing(1,:)
noise_wing = noise_age_wing(2,:)
mdl_noise = fitlm(noise_age,noise_wing)
plot(mdl_noise)

%effect: R Squared decreased
