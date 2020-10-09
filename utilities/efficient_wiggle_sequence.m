%%% script used to find an efficient set of ACC sequences for an
%%% n-dimensional simulation. A set of sequences should contain all
%%% different combinations of wiggles in all n-dimensions. The number
%%% of 2nd-order wrap-around artefact is minimized

clear
%% Parameters
n_dim = 3;    % number of dimensions the ACC is performed on

%% generate reference sequence (same starting sequence as used in wiggle_perm)
% Powers of minus one generate rows of alternating ones and minus ones.
% The alternation period is two on the first row, four on the second,
% eight on the third, etc. This generates all possible wiggle permutations 
% where every column is a unique permutation.
n_seq = 2^n_dim;
wiggles_ref = (-1).^(ceil((1:n_seq)./2.^((1:n_dim)' - 1)) + 1);

%% generate all posible sequences of column indices ((n_seq-1)! permutations)
all_perms = perms((2:n_seq)); % sequence always starts on 1 

%% iterate over all possible combinations
% initialization
c_best = 0;
var_best = inf;
for p = 1:factorial(n_seq-1) 
    % check one of the possible permutations
    ind = [1,all_perms(p,:)]; % always start sequence with column 1
    wiggles = wiggles_ref(:,ind);    
    
    % generate all 8 sequences in all 3 dimensions
    sequences = zeros(n_seq,n_seq,n_dim);
    for dim = 1:n_dim
        for s = 1:n_seq
            sequences(s,:,dim) = circshift(wiggles(dim,:),s);
        end
    end
    
    % check number of times the sequence cancels 2nd-order wrap-around
    % artefacts
    cancel_set = zeros(n_seq,sum(1:n_dim)); % list of canceled events
    i_dim = 0; % dimension counter
    for dim1 = 1:n_dim
        for dim2 = dim1:n_dim 
            i_dim = i_dim + 1;
            
            % check combination of sequences for number of ++, +-, -+ and -- terms
            for n = 1:n_seq
                seq = [sequences(:,1,dim1),sequences(:,n,dim2)];
                d = sum(abs(sum(seq,2)) == 2); % count number of -- or ++ terms
                if d == n_seq/2
                    % perfect cancelation of second order wrap-around artefact
                    % when number of ++/-- exactly equals +-/-+ terms
                    cancel_set(n,i_dim) = 1;
                elseif d == n_seq / 4 || d == n_seq * 3/4
                    % half cancelation of second order wrap-around artefact
                    cancel_set(n,i_dim) = 1/2;
                end
            end
        end
    end
    
    % check total number of 2nd-order artefacts canceled and check
    % difference in efficiency across different dimensions
    c_tot = sum(cancel_set(:));
    c_var = var(sum(cancel_set,1));
    
    % the wiggle sequence is chosen on basis of two optimization criteria
    % 1) most 2nd-order artefacts canceled
    % 2) balanced performance across different dimensions
    if c_tot > c_best || (c_tot == c_best && c_var <= var_best)
        wiggles_best = wiggles;            % currently best wiggle sequence 
        ind_best = ind;                    % column indices best sequence
        var_best = var(sum(cancel_set,1)); % variation efficiency across different dimensions
        cancel_set_best = cancel_set;      
        c_best = sum(cancel_set(:));
    end
end

