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
    wiggle_set(wiggle_flags == 1,:) = (-1).^(ceil((1:2^n_wiggles)./2.^((1:n_wiggles)' - 1)) + 1);
    
    % shuffle the order of the columns of wiggle_set to minimize the
    % second-order wrap-around artefact (see efficient_wiggle_sequences)
    ind = 1:2^n_wiggles;
    if n_wiggles == 2
        ind = [1,2,4,3];
    elseif n_wiggles == 3
        ind = [1,2,3,5,4,6,7,8];
    end
       
    wiggle_set = wiggle_set(:,ind);
end

