%simulation of stem cell clonal expansion
%testing with addition of memory
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
%number of competitive stem cells in the open layer (N in manuscript)
n_openniche = 100;

%define memory parameters
memory_gain_prob=0.5;
% memory_loss_prob=0.5; % add later but right now just fixed at 5
memory_strength_list=[5 10];

% Specify list of epsilons (master rate) and lambdas (competitive rate) to try
epsilon_list = [0 0.3 1];
lambda_list = [1 1 0];

% for each value of epsilon...
for i_mem = 1:length(memory_strength_list)
    % assign memory
    memory_strength = memory_strength_list(i_mem);
    for i_para = 1:length(epsilon_list)
        % assign epsilon and lambda
        epsilon = epsilon_list(i_para);
        lambda = lambda_list(i_para);
        disp("Currently processing epsilon " + epsilon + " and lambda " + lambda + " and memory strength " + memory_strength);
    
    %Reinitialize data matrix for every epsilon
        %matrix containing time series data
        x_matrix=zeros(num_of_clones,iter); % row = clone, col = timepoint
        %define matrix to retain memory info
        % col1 = clone cell belongs to
        % col2 = competitive or master
        % col3 = memory retained (binary)
        % col4 = number of iterations retained
        mem_matrix=zeros(n_openniche+num_of_clones, 4); 
        % assign clone IDs and stem cell label
        mem_matrix(:,1)=repelem(1:num_of_clones, num_of_clones+1);
        mem_matrix(:,2)=repmat([repelem(0,n_openniche/num_of_clones),1], 1, num_of_clones); %0 = competitive, 1 = master
    
        %initional condition
        init_size = floor(n_openniche/num_of_clones); %divide cells in open niche into the 10 clones
        x_matrix(:,1) = init_size*ones(num_of_clones,1); %set initial timepoint to init_size
        x_matrix(1:n_openniche-init_size*num_of_clones,1) = x_matrix(1:n_openniche-init_size*num_of_clones,1) + 1; %divys out remaining stem cells
    
        %time evolution
        for i=1:iter-1
    
            x_temp=x_matrix(:,i); %grab that iteration's column
        %lists of the probabilities that each clone will be selected for proliferation/differentiation.
            % note below summing number of memory retaining cells per clone to
            % multiply by memory strength
            f_pro_list=(lambda*x_temp+epsilon + memory_strength*groupsummary(mem_matrix(:,3), mem_matrix(:,1), @sum))/(lambda*n_openniche+num_of_clones*epsilon+(memory_strength*sum(mem_matrix(:,3))));
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
    
        % 1. change cell identities based on pro & diff events
            % choose random cell from diff_type clone to become pro_type
            diff_cell_randid = randi([1 x_temp(diff_type)+1]);
    
            % get index of that cell in the memory array
            diff_clone_indeces = find(ismember(mem_matrix(:,1), diff_type));
            diff_cell_index = diff_clone_indeces(diff_cell_randid);
    
            % convert to proliferative clone type
            mem_matrix(diff_cell_index,1) = pro_type;
    
            % erase any memory (since not actually same cell)
            mem_matrix(diff_cell_index, 3:4) = 0;
    
        % 2. account for any previously obtained memory
            % Add one iteration to any nonzero memorys
            mem_matrix(mem_matrix(:,4) ~= 0, 4) = mem_matrix(mem_matrix(:,4) ~= 0, 4) + 1; 
    
            % if number of iterations >= 5 -> set to 0
            % otherwise do nothing (constant memory)
            mem_matrix(mem_matrix(:,4) == 5, 3:4) = 0;
    
        % 3. Determine if memory is retained
            rand_num_4=rand; % grab random number b/w 0 and 1
            mem_type=(rand_num_4 <= memory_gain_prob); % 1 = memory retained, 0 = no memory retained
            % if memory retained...
            if mem_type
        % 4. Determine if master or comp stem cell divided 
                % Draw random number and determine if less than epsilon /
                % (epsilon + lambda * n_k) to see if master sc divided
                rand_num_5=rand; % grab random number b/w 0 and 1
                % if master stem cell divided...
                if rand_num_5 <= epsilon/(epsilon+lambda*x_temp(pro_type))
        % 5. Assign memory to a cell
                   % Isolate master stem cell for pro_type clone and update
                   mem_matrix(mem_matrix(:,1) == pro_type & mem_matrix(:,2) == 1, 3) = 1; % yes, memory
                   mem_matrix(mem_matrix(:,1) == pro_type & mem_matrix(:,2) == 1, 4) = 1; % first iteration
                else
                   % choose random competitive stem cell to be the progenitor
                   pro_cell_randid = randi([1 x_temp(pro_type)]);
                   pro_clone_indeces = find(ismember(mem_matrix(:,1), pro_type));
                   pro_cell_index = pro_clone_indeces(pro_cell_randid);
                   
                   % update their memory
                   mem_matrix(pro_cell_index, 3) = 1;
                   mem_matrix(pro_cell_index, 4) = 1;
                end
            end
        end
    
        %%save
        run_tag = strcat('_ep',num2str(epsilon),'_lm',num2str(lambda),'_K',num2str(num_of_clones),'_N',num2str(n_openniche),'_a',num2str(memory_strength));
        filename = ['/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/matrices/timeseries', run_tag , '.mat'];
        save(filename, 'x_matrix')
        save(strcat('/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/matrices/variables',run_tag,".mat"))
    end
end

toc %ends clock