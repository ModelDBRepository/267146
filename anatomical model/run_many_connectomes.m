
% run multiple m anatomical models: 
% for each iteration i of the model output files are saved as "file_i.txt"
% in the connectome files folder

m = 200; % number of anatomical connectomes to be simulated

parfor i=1:m
    main_spinal_cord_long_14(i)
end