%scaling law in clone size
%visualize RMSE b/w scaled probability curves at different timepoints
%written by AsahiNakamuta out of Kyoto University
%lightly annotated NP 10/11/2023

tic
clear variables
clc

% Find all timeseries matrices
input_dir = "/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/probability_curve_parameter_screen/matrices/";
mats = dir(fullfile(input_dir, "clone1_timeseries*.mat"));

% Initialize output matrix
% row = parameter set
% col = parameter value (l, e, N) + RMSE between probability curves (5-10, 10-15,
% 5-15)
rmse_mat = zeros(length(mats), 6);

% Define color map (Renoir colors from MetBrewer)
color_map = [
    25 18 82;
    % 112 90 163;
    158 154 218;
    255 176 174;
    243 123 108;
    207 37 26;
    % 242 153 0;
    255 186 49;
    % 174 165 23;
    40 90 32
    ];
color_map = color_map / 255;

% For each matrix...
for mat=1:length(mats)
    % Define filename
    disp(mat) % just for progress update
    mat_filename = mats(mat).name;

    % Load matrix
    x_matrix = importdata(input_dir + mat_filename);

    % Load in variables
    load(replace(strcat(input_dir + mat_filename), 'timeseries', 'variables'))

    % Initialize the matrices to store the clone sizes per time point
    temp = zeros(3,n_openniche);
    unscale = zeros(3,n_openniche);
    scale = zeros(3,n_openniche);
    x_axis = zeros(3,n_openniche);
    for t = 1:3
        % 5 because we set tmax to 15 in the simulation
        pickup = round(5*t*(epsilon*num_of_clones+lambda*n_openniche));
        % For each simulation
        for i = 1:m
            if x_clone1(i,pickup)~=0
                temp(t,x_clone1(i,pickup)) = temp(t,x_clone1(i,pickup))+1;
            end
        end
        % For each possible clone size...
        for j=1:n_openniche
            % Divide number in each clone size by number of clones for
            % proportion density
            unscale(t,j) = temp(t,j)/nnz(x_clone1(:,pickup));
            % Multiply by the average clone size
            scale(t,j) = temp(t,j)/nnz(x_clone1(:,pickup))*sum(x_clone1(:,pickup))/nnz(x_clone1(:,pickup));
            % Divide by the average clone size
            x_axis(t,j) = j/(sum(x_clone1(:,pickup))/nnz(x_clone1(:,pickup)));
        end
    end

    % Calculate RMSE between each time point combo (provided each clone
    % size has a nonzero value in both time points)
    prop_5 = scale(1,:);
    prop_10 = scale(2,:);
    prop_15 = scale(3,:);
    rmse_5_10 = rmse(prop_5(scale(1,:) > 0 & scale(2,:) > 0), prop_10(scale(1,:) > 0 & scale(2,:) > 0));
    rmse_5_15 = rmse(prop_5(scale(1,:) > 0 & scale(3,:) > 0), prop_15(scale(1,:) > 0 & scale(3,:) > 0));
    rmse_10_15 = rmse(prop_10(scale(2,:) > 0 & scale(3,:) > 0), prop_15(scale(2,:) > 0 & scale(3,:) > 0));

    % Update RMSE mat with parameters and RMSE values
    rmse_mat(mat,:) = [lambda epsilon n_openniche rmse_5_10 rmse_5_15 rmse_10_15];

    % % Produce scaled figure
    % fig = figure('visible','off');
    % fig.Position = [10 10 500 400]; 
    % semilogy(0,0)
    % hold on
    % for t = 1:3
    %     p = semilogy(x_axis(t,:),scale(t,:), "LineWidth", 2);
    %     p.Color = color_map(t,:);
    % end
    % hold off
    % 
    % % Add formatting
    % sgtitle(strcat("\epsilon = ", num2str(epsilon), ", \lambda = ", num2str(lambda), ", N = ", num2str(n_openniche), ": scaled"))
    % legend("", "time = 5", "time = 10", "time = 15")
    % ylabel('Probability density x avg. clone size');
    % xlabel('Clone size / avg. clone size');
    % 
    % % Save image (clone size across time)
    % fontsize(fig, 20, "points")
    % output_dir = "/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/probability_curve_parameter_screen/plots/";
    % plot_filename = output_dir + erase(mat_filename, ".mat") + "_scaled.png";
    % exportgraphics(fig, plot_filename, 'Resolution', 300)
end

% Find mean RMSE
rmse_mat(:,7) = mean(rmse_mat(:,4:6), 2);

% Calculate epsilon to lambda ratio
rmse_mat(:,8) = rmse_mat(:,2)./rmse_mat(:,1);

% List all tested values of N
n_list = [50 100 120 140 150 200 300 500];

% Remove Inf rows
rmse_mat_noinf = rmse_mat(sum(isfinite(rmse_mat), 2) == 8, :);

% Initialize scaled figure
fig = figure('visible','off');
fig.Position = [10 10 500 400]; 
hold on

for i=1:length(n_list)
    % Specify n
    n = n_list(i);

    % Isolate hierararchal RMSE value (lambda = 0)
    rmse_hi = rmse_mat(rmse_mat(:,1) == 0 & rmse_mat(:,3) == n, 7);

    % Isolate N value in matrix
    rmse_mat_subset = sortrows(rmse_mat_noinf(rmse_mat_noinf(:,3) == n,:),8);

    % Plot RMSE for N
    p = plot(rmse_mat_subset(:,8), rmse_mat_subset(:,7), '-ok', "LineWidth", 2);
    p.Color = color_map(i,:);

    % Add RMSE for hierarchal model as horizontal line
    p = yline(rmse_hi, '--', "LineWidth", 2);
    p.Color = color_map(i,:);
end
hold off

% Add formatting
sgtitle("Avg. RMSE between timepoints 5, 10, and 15")
legend("N = 50", "N = 100", "N = 150")
ylabel('avg. RMSE between scaled probability curves');
xlabel('\epsilon/\lambda');

% Save image (clone size across time)
fontsize(fig, 18, "points")
output_dir = "/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/probability_curve_parameter_screen/plots/";
plot_filename = output_dir + "avg_rmse_bw_scaled_probability_curves_throughE=4.png";
exportgraphics(fig, plot_filename, 'Resolution', 300)

toc