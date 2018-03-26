classdef Source
    % Object to represent one or more sources for a simulation
    %
    % to create a point source at pos x,y:   s = source([x,y], value);
    % to create a planar source at depth y:  s = source([1,y], ones(N,1));
    % etc.
    %
    % Sources can be added together with the + operator
    % for example:
    %
    % two_points = source([x1,y1], 1) + source([x2, y2], 1);
    %
    % Ivo M. Vellekoop 2018
    properties
        positions
        values
    end
    
    methods
        function obj = Source(position, value)
            % Creates a new source at the specified position
            % pos = relative position of source within roi of the simulation
            %       for vector simulations, specify a 4-element vector
            %       where the last element represents the polarization
            %       (1=x, 2=y, 3=z).
            %
            % value = scalar (for a point source), or any array (1-D, 2-D,
            %       3-D, 4-D) describing the amplitude and dimensions
            %       of the source.
            obj.positions = cell(1,1);
            obj.positions{1} = Source.make4(position);
            obj.values = cell(1,1);
            obj.values{1} = value;
        end
        
        function C = plus(A, B)
            C = A;
            C.positions = [A.positions, B.positions];
            C.values = [A.values, B.values];
        end
        
        function source_energy = energy(obj)
            % calculates the total energy of all sources
            % todo: properly account for overlapping sources!
            source_energy = 0;
            for c=1:numel(obj.positions)
                source_energy = source_energy + simulation.energy(obj.values{c});
            end
        end
        
        function source = prepare(obj, sim)
            % Shifts the sources to fit in the region of interest
            % If the source does not fit inside the roi, an error
            % is given.
            % This function also converts all sources to the proper format
            % (single / double precision, cpu / gpu)
            % todo: design issue, this function is coupled too tightly with
            % the simulation class
            %
            source = obj;
            for c=1:numel(source.positions)
                % shift source positions so that they are relative to the start of the roi
                source.positions{c} = source.positions{c} + sim.roi(1,:) - 1;
                if any((source.positions{c} + source.make4(size(source.values{c})) - 1) >  sim.roi(2,:))
                    error('Source does not fit inside the simulation');
                end
                
                % todo: clear source from GPU memory after usage
                % 
                source.values{c} = sim.data_array(source.values{c}); 
            end
        end
        
        function E = add_to(obj, E, roi, A)
            % adds the source to E
            % the source is clipped tothe specified roi, values outside of
            % the roi are discarded
            %
            % A is an amplitude prefactor to multiply the source by
            
            % process all sources
            for c=1:numel(obj.positions)
                % determine size of overlap
                pos = obj.positions{c};
                sz  = Source.make4(size(obj.values{c}));

                %calculate intersection of source and roi
                tlt = max(roi(1,:), pos); %top left corner target
                brt = min(roi(2,:), pos + sz - 1); %bottom right corner target
                if any(tlt > brt)
                    continue; %overlap is empty
                end
                tls = tlt - pos + 1; %top left corner source
                brs = brt - pos + 1; %bottom right corner source
                
                %TODO: the code below is really slow on a gpu because
                % of the indexing. Fortunately we don't call it often in 
                % wavesim, but for PSTD it is _the_ bottleneck
                %
                % it appears that it is already 50% faster if we convert the
                % indices to int64 first! (for a 1x400x1 array set)
                %%
                tls = int64(tls);
                brs = int64(brs);
                tlt = int64(tlt);
                brt = int64(brt);
        
                %add source, multiplied by prefactor A
                E(tlt(1):brt(1), tlt(2):brt(2), tlt(3):brt(3), tlt(4):brt(4)) =...
                        E(tlt(1):brt(1), tlt(2):brt(2), tlt(3):brt(3), tlt(4):brt(4)) + ...
                        A * obj.values{c}(tls(1):brs(1), tls(2):brs(2), tls(3):brs(3), tls(4):brs(4));
             end
        end
    end
    methods (Static)
        function sz = make4(sz)
            % converts the vector sz to a 4-element vector by appending 1's
            % when needed. This function is useful since 'size' removes
            % trailing singleton dimensions of arrays, so a 100x100x1x1
            % array returns a size of [100, 100], whereas a 100x100x1x3
            % array retuns a size of [100, 100, 1, 3]
            % As a workaround for this inconsistency, we always use
            % 4-element size vectors.
            sz = sz(:).';
            if numel(sz) < 4
                sz((end+1):4) = 1;
            end
        end
    end
end

