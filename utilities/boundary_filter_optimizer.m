%%% script used the design the optimal filter window function for wavesim, 
%%% minimizing the reflections at the boundary layer.

%% parameters
N = 128;     % total size window function (in pixels)
B = 32;     % size boundary layer (in pixels)

% weighting function  
% a quadratric weight function is used to minimize the higher frequencies 
% in the boundary layer. As a result, the reflections of the waves with a 
% close to normal incidence to the boundary layer are minimized most
% efficiently
k = fftshift(-N/2:N/2-1); % spatial frequencies
weights = abs(k).^2;      % weighting coefficients

%% Construct absorbing boundaries
% Optimizes coefficients so that a boundary with +1 on both sides
% reflects as little as possible, while having a minimum average
% value over this boundary. When the weight function is exactly quadratic,
% the optimum solution is a linear ramp

% set initial guess
coefs0 = ones(B,1);


% optimize boundary coeffiecients (coefficients can only between 0 and 1)
coefs = fmincon(@(c)f_cost_bump(c, weights), coefs0, [], [], [], [], zeros(size(coefs0)), ones(size(coefs0)));

% calculate cost function
cost0 = f_cost_bump(coefs0, weights); % cost function initial guess
cost = f_cost_bump(coefs, weights); % cost function optimized coefficients
disp(cost0 / cost)

% compare initial guess and optimized coefficients
plot(1:B,coefs0,'--k',1:B,coefs,'-b','LineWidth',2);
legend('initial guess','optimized','Location','SouthEast');
xlabel('boundary pixels'); ylabel('boundary coefficients');
set(gca,'FontSize',14);

% fit linear ramp through optimal coefficients
fit_fun = ['(x-a)./(',num2str(B),'+b)'];
fit((1:B)',coefs,fit_fun,'StartPoint',[0,0])

function cost = f_cost_bump(coefs, weight)
%% Calculate cost function
    N = length(weight);
    cost = weight * abs(fft(f_bump(coefs, N))).^2;
end

function bump = f_bump(coefs, N)
%% Construct bump function
    B = length(coefs);
    f_coefs = coefs(end:-1:1);
    bump = [...
        0
        coefs;...             % from 0 to 1
        ones(N-2*B-1, 1); ... % 1
        f_coefs;...             % from 1 to 0
        ];              % from -1 to 0
end