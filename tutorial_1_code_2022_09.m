clear

%Exercise 1
%operating via frequentist approach: given that the null hypothesis (person
%does not have HIV) is true, how likely is it that this data (a positive
%result)will happen?
%Answer:not statistically significant because p=.05 is not <.05


%exercise 2
%set my variables for sample number, prevalence range (.1 to 1) and
%specificity
likelihood_false_positive_by_probability = zeros(10,1)
prev_dis = 0:0.1:1
sample_number = 1000
prob_false_positive = [.95 .05]

%create a for loop to cycle through disease prevalence range
for i = 1:11
    %set the probabilty for not having disease given this prevalence and
    %the probability for having this disease
    prob_disease = [(1-prev_dis(1,i)) prev_dis(1,i)]
    %0 means "No," 1 means "yes" has disease
    has_disease = [0 1]
    %create a random distribution of 1000 samples based on prevalence for
    %whether someone has HIV
    disease_distribution = randsrc(sample_number,1,[has_disease; prob_disease]) 

    %create an array that decides whether a given sample should be a false
    %positive
    %1 signifies "yes," a false positive
    false_positive = [0 1]
    false_positive_distribution = randsrc(sample_number,1,[false_positive; prob_false_positive])

    %create an array of the actual results of tests
    %0 means negative test, 1 means positive test
    test_result = zeros(sample_number,1)
    for j = 1:sample_number
     if disease_distribution(j,1) == 1
        test_result(j,1) = 1
        elseif false_positive_distribution(j,1) == 1
            test_result(j,1) = 1
   
        end
    end

%now have an array of disease distribution and test results
%can calculate likelihood of false positive given this disease distribuiton
%by comparing the arrays
false_positive_number = 0
for j = 1:sample_number
    if test_result(j,1) == 1 & test_result(j,1) ~= (disease_distribution(j,1))
        false_positive_number = false_positive_number + 1
        
    end
end
false_positive_percentage = false_positive_number/sample_number
likelihood_false_positive_by_probability(i,1) = false_positive_percentage
end

%% 
%generate a table illustrating results
Prevalence_of_disease = prev_dis'
Proportion_infected_given_pos_result = 1-likelihood_false_positive_by_probability
summary_table = table(Prevalence_of_disease,likelihood_false_positive_by_probability,Proportion_infected_given_pos_result)

        



