%simulation of stem cell clonal expansion
%written by AsahiNakamuta out of Kyoto University
%lightly annotated NP 10/11/2023

tic %starts clock
clear variables
clc %clears all content from command window

% set seed
rng("default")

%initial conditions and parameters

%number of iteration
iter=103000;
%number of clones (K in manuscript, also number of master stem cells)
num_of_clones=10;
% set proliferation rates
epsilon = 0.3;
lambda = 1;

% Specify list of number of competitive stem cells to simulate
n_openniche_list = [10 20 50 100 200 500];

% for each value of n_openniche...
for i_para = 1:length(n_openniche_list)
    % assign n_open_niche
    n_openniche = n_openniche_list(i_para);
    disp("Currently processing " + n_openniche + " competitive stem cells");

%Reinitialize data matrix for every epsilon
    %matrix containing time series data
    x_matrix=zeros(num_of_clones,iter); % row = clone, col = timepoint

    %initional condition
    init_size = floor(n_openniche/num_of_clones); %divide cells in open niche into the 10 clones
    x_matrix(:,1) = init_size*ones(num_of_clones,1); %set initial timepoint to init_size
    x_matrix(1:n_openniche-init_size*num_of_clones,1) = x_matrix(1:n_openniche-init_size*num_of_clones,1) + 1; %divys out remaining stem cells

    %time evolution
    for i=1:iter-1

        x_temp=x_matrix(:,i); %grab that iteration's column
    %lists of the probabilities that each clone will be selected for proliferation/differentiation.
        f_pro_list=(lambda*x_temp+epsilon)/(lambda*n_openniche+num_of_clones*epsilon);
        f_diff_list=x_temp/n_openniche;
    %choose one clone to proliferate
        % create list of cumulative probabilities
        threshold_pro_list=zeros(num_of_clones-1,1);
        threshold_pro_list(1)=f_pro_list(1);
        for k=2:num_of_clones-1
            threshold_pro_list(k)=threshold_pro_list(k-1)+f_pro_list(k);
        end

        % grab random number b/w 0 and 1 and find which interval it fits into
        rand_num1=rand;
        pro_type=sum(threshold_pro_list<=rand_num1)+1;

        % update with additional cell
        x_temp(pro_type)=x_temp(pro_type)+1;

    %choose one clone to differeniate
        % create list of cumulative probabilities
        threshold_diff_list=zeros(num_of_clones-1,1);
        threshold_diff_list(1)=f_diff_list(1);

        for k=2:num_of_clones-1
            threshold_diff_list(k)=threshold_diff_list(k-1)+f_diff_list(k);
        end

        rand_num2=rand;
        diff_type=sum(threshold_diff_list<=rand_num2)+1;

        % update with removing cell
        x_temp(diff_type)=x_temp(diff_type)-1;

        x_matrix(:,i+1)=x_temp;
    end

    %%save
    run_tag = strcat('_ep',num2str(epsilon),'_lm',num2str(lambda),'_K',num2str(num_of_clones),'_N',num2str(n_openniche));
    filename = ['results/hierarchal_neutral_competition/Fig2/matrices/timeseries', run_tag , '.mat'];
    save(filename, 'x_matrix')
    save(strcat('results/hierarchal_neutral_competition/Fig2/matrices/variables',run_tag,".mat"))
end

toc %ends clock