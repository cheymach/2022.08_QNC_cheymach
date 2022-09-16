
%generate two normal distributions to simulate the the two cell types in
%the uploaded histogram

%Notes regarding the selected distribution
%using Figure 1c from 
% unsure what statistical Sarkar A, Mei A, Paquola ACM, Stern S, Bardy C, Klug JR, et al. Efficient Generation of CA3 Neurons from Human Pluripotent Stem Cells Enables Modeling of Hippocampal Connectivity In Vitro. Cell Stem Cell. 2018;22(5):684-97.e9.tests the paper used two compare these two
%populations, but I am assuming they they assumed normal distributions
%estimating mu and sigma by looking at the approximate mean and SEM of the
%could not find n so assuming an n of 100 for each 

clear

mu_hpnpc = 52
sem_hpnpc = 13
n_hpnpc = 100;
sigma_hpnpc = sem_hpnpc*sqrt(n_hpnpc)

mu_pan = 5
sem_pan= 3
n_pan = 100
sigma_pan = sem_pan*sqrt(n_pan)

% Get samples
samples_hpnpc = normrnd(mu_hpnpc, sigma_hpnpc, n_hpnpc, 1);
samples_pan = normrnd(mu_pan, sigma_pan, n_pan, 1);

%%
% create comparison bar plot with error bars to represent these two groups;
mean_hpnpc_s = mean(samples_hpnpc)
sem_hpnpc_s = std(samples_hpnpc)/sqrt(n_hpnpc)
err_high_hpnpc = mean_hpnpc_s + sem_hpnpc_s
err_low_hpnpc = mean_hpnpc_s - sem_hpnpc_s



mean_pan_s = mean(samples_pan)
sem_pan_s = std(samples_pan)/sqrt(n_pan)
err_high_pan = mean_pan_s + sem_pan_s
err_low_pan = mean_pan_s - sem_pan_s

mean_data = [mean_hpnpc_s mean_pan_s]
err_low = [err_low_hpnpc err_low_pan]
err_high = [err_high_hpnpc err_high_pan]

x = [1 2]

bar(x,mean_data)                

hold on

er = errorbar(x,mean_data,err_low,err_high);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
