function obj = SampleMedium(refractive_index, options, boundaries)
      % SampleMedium - Generation of refractive index map of a medium with absorption boundary layer
      % Saroch Leedumrongwatthanakun 2015
      % SampleMedium object:
           %   boundaries.width = extra space to add around simulation to simulate absorbing boundaries, in pixels
           %   boundaries.xcurve = shape of the damping curve horizontal boundary
           %   boundaries.ycurve = shape of the damping curve vertical boundary
           obj.N = size(refractive_index);
           
            %% Determine constants based on refractive_index map
            obj.PPW = options.lambda/options.pixel_size;
			obj.n_min = min(abs(refractive_index(:))); 
            obj.n_max = max(abs(refractive_index(:))); 
            obj.n_center = sqrt((obj.n_max^2 + obj.n_min^2) / 2); %central refractive index (refractive index that k_r is based on)
          
            
           %% fill in boundary settings
           if ~isfield(boundaries,'free') % if no free gap parameters are given
               boundaries.free=[0 0];
           end
           if ~isfield(boundaries,'width') % if no bonudary parameters are given
               boundaries.width = floor( (2.^nextpow2(size) - size)/4); % width of the absorbing boundary
               boundaries.xcurve = 1-linspace(0, 1, boundaries.width(2)).^2;   % damping curve of horizontal absorbing boundary
               boundaries.ycurve = 1-linspace(0, 1, boundaries.width(1)).^2;   % curve of vertical absorbing layer
           end
                    
           if ~isfield(boundaries,'xcurve') || boundaries.width(2) ~= length(boundaries.xcurve) % if width is given but not xcurve
               boundaries.xcurve = 1-linspace(0, 1, boundaries.width(2)).^2;
           end
           
           if ~isfield(boundaries,'ycurve') || boundaries.width(1) ~= length(boundaries.ycurve) % if width is given but not ycurve
               boundaries.ycurve = 1-linspace(0, 1, boundaries.width(1)).^2;
           end
          
           
            %% Setup grid, taking into account required boundary; adding minimum free space and padding to next power of 2 when needed
            % Determine the problem space boundaries including boundary
            % Determining the problem space size in grid.N
            obj.grid = simgrid(size(refractive_index)+2*boundaries.width+boundaries.free, options.pixel_size); 
            
            %% Padding refractive index to match simulation grid
            obj.grid.padding=obj.grid.padding+boundaries.free;% adding minimum free space
            obj.refractive_index = padarray(refractive_index, round((2*boundaries.width+obj.grid.padding)/2), 'replicate', 'both');
            obj.refractive_index = circshift(obj.refractive_index,round(-(2*boundaries.width+obj.grid.padding)/2));
            
                      
			%% Determine optimum value for epsilon, otherwise use given epsilon
            obj.BCepsmin = (56.43/options.lambda/(max(boundaries.width)*options.pixel_size)-0.2468); %empirical formular for threshold width of g0(r)=1e-3                                      
            
            %% Absorbing boundaries
            obj.damping_x = [ ones(1, obj.grid.N(2)-obj.grid.padding(2)-2*boundaries.width(2)), boundaries.xcurve, zeros(1, obj.grid.padding(2)), fliplr(boundaries.xcurve)];
            obj.damping_y = [ ones(1, obj.grid.N(1)-obj.grid.padding(1)-2*boundaries.width(1)), boundaries.ycurve, zeros(1, obj.grid.padding(1)), fliplr(boundaries.ycurve)];
                        
            %% Apply absorbing boundaries on potential map
            obj.refractive_index = obj.refractive_index(1:size(obj.damping_y,2),1:size(obj.damping_x,2));
            %obj.refractive_index = obj.refractive_index(1:size(obj.damping_y,2),1:size(obj.damping_x,2)) .* (obj.damping_y' * obj.damping_x);
            
           
            %obj.refractive_index = obj.refractive_index(1:size(obj.damping_y,2),1:size(obj.damping_x,2)) + 1.0i .* (obj.damping_y' * obj.damping_x) - 1.0i;
end
  