%Convert R.matlab's matrices to sparse matrices to save space and speed.
folder_prefix = 'Data'; % Or Data
file_prefix = 'RandomPartition'; %Does not appear to be used.

disp(pwd)

%Does not work in 2015
%dirs = clean_dir(fullfile([folder_prefix, '*'], [file_prefix, '*']), 1);
%Need to instead retrieve folders and concatenate them by hand.
dirs1 = dir(fullfile([folder_prefix, '*']));

if isempty(dirs1) || all(~dirs1.isdir)
    cd ..
    dirs1 = dir(fullfile([folder_prefix, '*']));
end

disp(dirs1)

dirs2 = {};
for i = 1:length(dirs1)
    tempdirs = dir(dirs1(i).name);
    tempdirs = tempdirs(~ismember({tempdirs.name},{'.','..'}));
    tempdirs = {tempdirs.name};
    tempdirs = fullfile(dirs1(i).name, tempdirs);
    dirs2 = [dirs2, tempdirs];
end

%dirs3 = {};
%for i = 1:length(dirs2)
%    tempdirs = dir(dirs2{i});
%    tempdirs = tempdirs(~ismember({tempdirs.name},{'.','..'}));
%    tempdirs = {tempdirs.name};
%    tempdirs = fullfile(dirs2{i}, tempdirs);
%    dirs3 = [dirs3, tempdirs];
%end

skip_rows = 1; %Number of rows to skip when joining matrices.

save_file_suffix = '_S.mat';

disp(['length: ', num2str(length(dirs2))])
%For each .mat in dirs, we need to convert its contents to sparse .mat.
%Then, we need to resave, but with the new name.
for d = 1 : length(dirs2)
    disp(d)
    %Folders do not exist in 2015.
    %dirs_path = fullfile(dirs2(d).folder, dirs2(d).name);
    %Does not work in 2015.
    %files = clean_dir(fullfile(dirs_path, '*.mat'), 0);
    
    files = dir(fullfile(dirs2{d}, '*.mat'));
    files = files(~cell2mat({files.isdir}));    

    disp(['number of files: ', num2str(length(files))])
    for f = 1 : length(files)
        file_name = fullfile(dirs2{d}, files(f).name);
        [file_part_1, file_part_2] = fileparts(file_name);
        file_name_s = [file_part_2, save_file_suffix];

        if(~isempty(strfind(file_name, save_file_suffix)) || ...
           any(ismember({files.name}, {file_name_s})))
            %If the file contains _S.mat or
            %If the file has an _S.mat version
            continue;
        end

        try
            %Load
            loaded = load(file_name);
            lnames = fieldnames(loaded);
            
            %Stitch Together: Pre-allocate. Note we are stacking matrices.
            %Now need three matrices: one for population, one for CFlux,
            %one for FFlux.
            %Rows are the same for each.
            %Flux cols are one less.
            %Number nonzero is obviously different.
            rows = 0;
            cols = 0;
            nonzPop = 0;
            nonzCFx = 0;
            nonzFFx = 0;
            
            %Need to retrieve the matrices that actually correspond to
            %population matrices from the cell array of all.
            %These are the ones without a CFlux or FFlux suffix.
            %https://uk.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string#answer_3240
            %find(contains(lnames, 'CFlux'))
            cfluxIndex = find(not(cellfun('isempty', strfind(lnames, 'CFlux'))));
            ffluxIndex = find(not(cellfun('isempty', strfind(lnames, 'FFlux'))));
            popIndex = (0:(length(cfluxIndex) - 1))*3 + 1;
            
            for n = 1 : length(popIndex)
                sizes = size(loaded.(lnames{popIndex(n)}));
                rows = rows + sizes(1) - skip_rows;
                cols = max(cols, sizes(2));
                nonzPop = nonzPop + nnz(loaded.(lnames{popIndex(n)}));
                nonzCFx = nonzCFx + nnz(loaded.(lnames{cfluxIndex(n)}));
                nonzFFx = nonzFFx + nnz(loaded.(lnames{ffluxIndex(n)}));
            end
            stitchPop = spalloc(rows, cols, nonzPop); %As recommended by MATLAB.
            stitchCFx = spalloc(rows, cols - 1, nonzCFx);
            stitchFFx = spalloc(rows, cols - 1, nonzFFx); 
            
            %Stitch Together: Perform Joining.
            sindex = 1;
            for n = 1 : length(popIndex)
                sizes = size(loaded.(lnames{popIndex(n)}));
                stitchPop((sindex : (sindex + sizes(1) - skip_rows - 1)),...
                       (1 : sizes(2))) = ...
                       loaded.(lnames{popIndex(n)})(1 + skip_rows:end, :);
                stitchCFx((sindex : (sindex + sizes(1) - skip_rows - 1)),...
                       (1 : (sizes(2)- 1))) = ...
                       loaded.(lnames{cfluxIndex(n)})(1 + skip_rows:end, :);
                stitchFFx((sindex : (sindex + sizes(1) - skip_rows - 1)),...
                       (1 : (sizes(2)- 1))) = ...
                       loaded.(lnames{ffluxIndex(n)})(1 + skip_rows:end, :);
                sindex = sindex + sizes(1) - skip_rows;
            end
            
            %Save stitched result.
            save(fullfile(file_part_1, file_name_s), ...
                'stitchPop', 'stitchCFx', 'stitchFFx');
            
            
        %https://www.mathworks.com/matlabcentral/answers/
        %325475-display-error-message-and-execute-catch#answer_255132
        catch e %e is an MException struct
            fprintf(1, 'Error: Identifier:\n%s\n', e.identifier);
            fprintf(1, 'Error: Message:\n%s\n', e.message);
            continue
        end

        clear loaded;
    end
end
