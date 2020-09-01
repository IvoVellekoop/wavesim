function wiggle_set = wiggle_perm(wiggle_flags)
    %%% function returns matrix with all possible combinations of enabled
    %%% wiggle directions
    %%% Variables
    %%% wiggle_flags:   vector indicating which direction will be wiggled (size: n_directions x 1)
    %%% wiggle_set:     matrix with every possible combination of enabled
    %%%                 wiggle direction in each column
    n_directions = length(wiggle_flags);    % number of wiggle directions considered
    n_wiggles = sum(wiggle_flags);          % number of enabled wiggle flags
    wiggle_set = zeros(n_directions,2^n_wiggles);

    % Powers of minus one generate rows of alternating ones and minus ones.
    % The alternation period is two on the first row, four on the second,
    % eight on the third, etc. This generates all possible wiggle permutations 
    % where every column is a unique permutation.
    wiggle_set(wiggle_flags == 1,:) = (-1).^(ceil([1:2^n_wiggles]./2.^([1:n_wiggles]' - 1)) + 1);
    
    % change the order of the columns of wiggle_set to remove any shift 
    % symmetries in the sequences
    ind = [1:2:2^n_wiggles-1,2^n_wiggles:-2:2];
    wiggle_set = wiggle_set(:,ind);
end

