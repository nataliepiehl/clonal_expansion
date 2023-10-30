%data-generation for scaling law in clone size
%goal: test different combinations of epsilon, lambda, and N to see where
%scaling law starts to fall apart
%written by AsahiNakamuta out of Kyoto University
%lightly annotated NP 10/11/2023

tic
clear variables
clc

%initial conditions and parameters

%define output folder
output_dir = "/projects/b1042/GoyalLab/nat/clonal_expansion/results/hierarchal_neutral_competition/probability_curve_parameter_screen/matrices/";
%number of simulation trials for each parameter set
m = 100000;
%number of clones (K in manuscript)
num_of_clones=10;

% Specify list of epsilons (master rate), lambdas (competitive rate), and
% total open niche cells to try
% start with iterations of 0.1 for the rates
epsilon_list = [0:0.1:4 1];
lambda_list = [repelem(1,length(0:0.1:4)) 0];
n_openniche_list = [120 140 200 300 500];

% for each value of n_openniche...
for i_n = 1:length(n_openniche_list)
    % assign n_openniche
    n_openniche = n_openniche_list(i_n);

    % for each value of epsilon...
    for i_para = 1:length(epsilon_list)
        % assign epsilon and lambda
        epsilon = epsilon_list(i_para);
        lambda = lambda_list(i_para);
        disp("Currently processing epsilon " + epsilon + " and lambda " + lambda + " with " + n_openniche + " competitive stem cells");
    
        %number of iteration - curious why they didn't define number of iterations
        %this way in fig 2
        iter=round(15*(epsilon*num_of_clones+lambda*n_openniche));
    
        %matrix containing time series data
        % this is the matrix to update for each simulation trial
        x_matrix=zeros(num_of_clones,iter);
        %timeseries of clone #1
        % this holds the clone sizes of only the first tested clone after the simulation
        % are any of the other clones even used or necessary?
        x_clone1 = zeros(m,iter);
    
        for j = 1:m
            %initional condition - disperse competitive stem cells amongst clones
            init_size = floor(n_openniche/num_of_clones);
            x_matrix(:,1) = init_size*ones(num_of_clones,1);
            x_matrix(1:n_openniche-init_size*num_of_clones,1) = x_matrix(1:n_openniche-init_size*num_of_clones,1) + 1;
        
            % equal probability of a competitive stem cell vs master stem cell gets
            % chosen
            rand_num0 = rand;
            label = (rand_num0 <= 1/(x_matrix(1,1)+1));   %%%% rand_num0 <= 1/(init_size+1) corresponds to random label. 0 corresponds to label on competitive cells. 1 corresponds to label on master stem cell.
            if label
                % If choose master stem cell, then no open niche cells (arbitrarily
                % move them elsewhere to retain same starting open niche number))
                x_matrix(2,1)=x_matrix(2,1) + x_matrix(1,1);
                x_matrix(1,1)=0;
            else
                % If choose competitive stem cell, then start with one open niche
                % cell and arbitrarily move the rest elswehere
                x_matrix(2,1)=x_matrix(2,1) + x_matrix(1,1) - 1;
                x_matrix(1,1)=1;
            end
        
            %time evolution - proceed through the same simulation as in Fig2
            for i=1:iter-1 % why -17? removing the -17 cause it messes with getting the last time point and don't know why it's there anyway
                
                x_temp=x_matrix(:,i);
            %lists of the probabilities that each clone will be selected for proliferation/differentiation.
                if label
                    f_pro_list=(lambda*x_temp+epsilon)/(lambda*n_openniche+num_of_clones*epsilon);
                else
                    f_pro_list=(lambda*x_temp+epsilon)/(lambda*n_openniche+(num_of_clones-1)*epsilon);
                    f_pro_list(1)=lambda*x_temp(1)/(lambda*n_openniche+(num_of_clones-1)*epsilon);
                end
                f_diff_list=x_temp/n_openniche;
            %choose one clone to proliferate
                threshold_pro_list=zeros(num_of_clones-1,1);
                threshold_pro_list(1)=f_pro_list(1);
                
                for k=2:num_of_clones-1
                    threshold_pro_list(k)=threshold_pro_list(k-1)+f_pro_list(k);
                end
                
                rand_num1=rand;
                pro_type=sum(threshold_pro_list<=rand_num1)+1;
                
                x_temp(pro_type)=x_temp(pro_type)+1;
            %choose one clone to differeniate
                threshold_diff_list=zeros(num_of_clones-1,1);
                threshold_diff_list(1)=f_diff_list(1);
                
                for k=2:num_of_clones-1
                    threshold_diff_list(k)=threshold_diff_list(k-1)+f_diff_list(k);
                end
                
                rand_num2=rand;
                diff_type=sum(threshold_diff_list<=rand_num2)+1;
                
                x_temp(diff_type)=x_temp(diff_type)-1;
                
                x_matrix(:,i+1)=x_temp;
                
            end
            x_clone1(j,:) = x_matrix(1,:);
        end
    
        %%save
        run_tag = strcat('_ep',num2str(epsilon),'_lm',num2str(lambda),'_K',num2str(num_of_clones),'_N',num2str(n_openniche));
        filename = strcat(output_dir, 'clone1_timeseries', run_tag , '.mat');
        save(filename, 'x_clone1', '-v7.3')
        save(strcat(output_dir, 'clone1_variables',run_tag,".mat"), '-v7.3')
    
    end

end

toc