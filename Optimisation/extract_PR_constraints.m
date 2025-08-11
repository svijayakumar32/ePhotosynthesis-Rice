% Read in all 129 optimization results, average the results and export PR
% enzyme Vmax values to PR_constraints to serve as constraints for full optimizations

% Get absolute path of this script
ScriptPath = mfilename('fullpath');

% Locate directory containing this script
ScriptDir = fileparts(ScriptPath);

% Define the folder containing the BestMatrix text files at 129 umol mol-1
results_folder = fullfile(ScriptDir,'Results', 'Enzymes');

% List all text files inside the folder
all_file_list = dir(fullfile(results_folder, '*.txt'));

% Get the list of 10 text files containing enzyme Vmax of optimizations at 129
% umol mol-1 using a reg exp to match names containing _129_ and numbers 1-9 or 10
%optimized_129_list = all_file_list(~cellfun('isempty', regexp({all_file_list.name}, '^outputenz_129_\d{1,2}\.txt$')));
optimized_129_list = all_file_list(~cellfun('isempty', regexp({all_file_list.name}, '^outputenz_129_(10|[1-9])_100\.txt$')));

% Initialize an empty cell array to store the matrices
optimized_129_matrices = cell(numel(optimized_129_list), 1);

% Loop over the files
for i = 1:numel(optimized_129_list)
    % Read the current file using readmatrix
    optimized_129_matrix = readmatrix(fullfile(results_folder, optimized_129_list(i).name));
    
    % Store the matrix in the cell array
    optimized_129_matrices{i} = optimized_129_matrix;
end

% Combine all cells in the array 
combined_129_matrix = horzcat(optimized_129_matrices{:,1});
blank_rows = zeros(1,10);
combined_129_matrix_new = vertcat(combined_129_matrix(2:8,:), ...
    blank_rows,combined_129_matrix(9,:),blank_rows,combined_129_matrix(10:11,:),blank_rows,combined_129_matrix(12:25,:));

% Average across the columns for each enzyme Vmax
avg_optimized_129 = mean(combined_129_matrix,2);

% Extract rows for photorespiratory constraints and save in text file
PR_constraints = avg_optimized_129(12:18);
% Create output filepath
PR_output_filepath = fullfile(ScriptDir,'PR_constraints.txt');
% Save PR constraints to output file
writematrix(PR_constraints,PR_output_filepath);
