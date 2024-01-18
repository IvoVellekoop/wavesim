function [] = SimulationCallback(callbackType,diff_energy, energy, energy_threshold, energy_domain, ax_label)
%Simulation callback function. Callable from Mex
%   callbackType - The name of the callback ('Default','NoCallback' or
%   'CrossImage')
if(strcmp(callbackType, 'NoCallback'))
    return;
end

if(strcmp(callbackType,'EnergyAdded'))
    if(nargin < 4)
        error("Not enough input arguments for EnergyAdded callback.")
    end
    
    figure(1); clf;
    de = diff_energy / diff_energy(1);
    threshold = energy_threshold;
    plot(1:length(de),log10(de),'b',[1,length(de)],log10(threshold)*ones(1,2),'--r');
    title(length(de));  xlabel('# iterations'); ylabel('log_{10}(energy added)');
end

if(strcmp(callbackType,'EnergyAddedDisp'))
    if(nargin < 4)
        error("Not enough input arguments for EnergyAddedDisp callback.")
    end
    
    disp(['Iter: ', int2str(length(diff_energy)), '. Energy added: ', ...
        num2str(diff_energy(end)/diff_energy(1)), '. Threshold: ', num2str(energy_threshold)]);
end

if(strcmp(callbackType,'Default'))
    %default callback function. Shows  total energy evolution and 
    %real value of field along specified dimension (longest
    %dimension by default)
    %
    if(nargin < 4)
        error("Not enough input arguments for Default callback.")
    end

    figure(1); clf;
    de = diff_energy / diff_energy(1);
    threshold = energy_threshold;
    subplot(2,1,1); plot(1:length(de),log10(de),'b',[1,length(de)],log10(threshold)*ones(1,2),'--r');
    title(length(de));  xlabel('# iterations'); ylabel('log_{10}(energy added)');

    subplot(2,1,2); 
    sig = log(abs(energy));
    if(nargin >= 5)
        plot(energy_domain, squeeze(sig),'b'); hold on;
    else
        plot(squeeze(sig),'b'); hold on;
    end
    title('midline cross-section')
    ylabel('log(|E|)');
    if(nargin >= 6)
       xlabel(ax_label); 
    end
    
    drawnow;
end

if(strcmp(callbackType,'CrossImage'))
    if(nargin < 3)
        error("Not enough input arguments for CrossImage callback.")
    end
    
    figure(1); clf;
    imagesc(energy);
    axis image;
    title(['Differential energy ' num2str(diff_energy(end)/diff_energy(1))]);
    drawnow;
end

end

