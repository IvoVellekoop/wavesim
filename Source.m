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
        positions %don't access directly, implementation may change!
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
        
        function source = shift(obj, pos)
            % Shifts the sources so that an original source position of
            % [1,1,1,1] now corresponds to 'pos'
            %
            source = obj;
            for c=1:numel(source.positions)
                source.positions{c} = source.positions{c} + pos - 1;
            end
        end

        function source = crop(obj, roi)
            % Removes all parts of the sources that fall outside of the roi
            source = obj;
            for c=1:numel(source.positions)
                pos = source.positions{c};
                sz  = Source.make4(size(source.values{c}));
                tlt = max(roi(1,:), pos); %top left corner target
                brt = min(roi(2,:), pos + sz - 1); %bottom right corner target
                if any(tlt > brt)
                    continue; %overlap is empty
                end
                tls = tlt - pos + 1; %top left corner source
                brs = brt - pos + 1; %bottom right corner source
                source.values{c} = source.values{c}(tls(1):brs(1), tls(2):brs(2), tls(3):brs(3), tls(4):brs(4));
                source.positions{c} = tlt;
            end
        end

        function source = gpuArray(obj)
            % Places the source data on the gpu
            source = obj;
            for c=1:numel(source.positions)
                source.values{c} = gpuArray(source.values{c});
            end
        end

         function source = cast(obj, type)
            % Converts the source data to the specified type ('single' or 'double')
            source = obj;
            for c=1:numel(source.positions)
                source.values{c} = cast(source.values{c}, type);
            end
        end

        function E = add_to(obj, E, A)
            % pre-multiplies the source by A and adds it to E
            
            % process all sources
            for c=1:numel(obj.positions)
                % determine size of overlap
                pos = obj.positions{c};
                sz  = Source.make4(size(obj.values{c}));
               %TODO: the code below is really slow on a gpu because
                % of the indexing. Fortunately we don't call it often in 
                % wavesim, but for PSTD it is _the_ bottleneck
                %
                % it appears that it is already 50% faster if we convert the
                % indices to int64 first! (for a 1x400x1 array set)
                %%
                tlt = int64(pos);
                brt = int64(pos + sz - 1);
        
                %add source, multiplied by prefactor A
                E(tlt(1):brt(1), tlt(2):brt(2), tlt(3):brt(3), tlt(4):brt(4)) =...
                        E(tlt(1):brt(1), tlt(2):brt(2), tlt(3):brt(3), tlt(4):brt(4)) + ...
                        A * obj.values{c};
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

