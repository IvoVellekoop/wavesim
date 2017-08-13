classdef wavesimv < wavesim
    %Simulation of the 2-D or 3-D Maxwell's equations wave equation using a Born series approach
    %Only valid for non-magnetic media.
    % Ivo M. Vellekoop 2017
    properties
        pxr; %p_x / sqrt(k_0^2 + i epsilon)        (used to evaluate dyadic Green's function)
        pyr;
        pzr;
    end
    methods
        function obj=wavesimv(sample, options)
            obj@wavesim(sample, options);
            
            red = 1/sqrt(obj.k^2 + 1.0i * obj.epsilon); 
            obj.pxr = sample.grid.px_range * red;
            obj.pyr = sample.grid.px_range * red;
            obj.pzr = sample.grid.px_range * red;
        end
        
        function state = run_algorithm(obj, state)
            %% Allocate memory for calculations
            state.E  = data_array(obj);    
            state.Ey = data_array(obj);    
            state.Ez = data_array(obj);    
            
            %% simulation iterations
            while state.has_next
                % Apply the potential function and add the source to all 3 field components
                Etx = obj.V.* state.E;
                Etx(state.source_range{1}, state.source_range{2}, state.source_range{3}) = Etx(state.source_range{1}, state.source_range{2}, state.source_range{3}) + state.source(1,:,:,:);
                Etx = fftn(Etx);

                Ety = obj.V.* state.Ey;
                Ety(state.source_range{1}, state.source_range{2}, state.source_range{3}) = Etx(state.source_range{1}, state.source_range{2}, state.source_range{3}) + state.source(2,:,:,:);
                Ety = fftn(Ety);
                
                Etz = obj.V.* state.Ez;
                Etz(state.source_range{1}, state.source_range{2}, state.source_range{3}) = Etx(state.source_range{1}, state.source_range{2}, state.source_range{3}) + state.source(3,:,:,:);
                Etz = fftn(Etz);
                
                %calculate divergence term in dyadic Green function:
                Ediv = Etx .* obj.pxr + Ety .* obj.pyr + Etz .* obj.pzr;
                
                % calculate gradient of the divergence (the p p^T term) and
                % subtract it from E
                Etx = Etx - obj.pxr .* Ediv;
                Ety = Ety - obj.pxy .* Ediv;
                Etz = Etz - obj.pxz .* Ediv;
                
                Ediff = (1.0i/obj.epsilon*obj.V) .* (state.E-ifftn(obj.g0_k .* Etx));
                Ediffy = (1.0i/obj.epsilon*obj.V) .* (state.Ey-ifftn(obj.g0_k .* Ety));
                Ediffz = (1.0i/obj.epsilon*obj.V) .* (state.Ez-ifftn(obj.g0_k .* Etz));
           
                if state.calculate_energy
                   state.last_step_energy = simulation.energy(Ediff(obj.roi{1}, obj.roi{2}, obj.roi{3}));
                end
                
                state.E = state.E - Ediff;
                state.Ey = state.Ey - Ediffy;
                state.Ez = state.Ez - Ediffz;
                state = next(obj, state);
            end
        end
    end
end