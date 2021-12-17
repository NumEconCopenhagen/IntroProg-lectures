classdef model
    methods (Static)
        
        % model solution
        function c = c_star(w,e,par)
            c = par.gamma*(1.0-par.tau)*w + par.gamma*e;
        end
        
        function l = l_star(w,e,par)
            l = (1.0-par.gamma) + (1.0-par.gamma)/(1.0-par.tau)*(e./w);
        end
        
        % objective function
        function Q = obj_fun(theta,est_par,w,mom_data,moments_fun,par)
            
            % a. update parameter struct with elements in theta
            for p = 1:numel(theta)
                par.(est_par{p}) = theta(p);
            end
            
            % b. draw random draws
            rng(893245); % note: different from the seed used for the true data
            S = 100; % number of simulation draws for each observed individual

            n = numel(w); % number of individuals
            w = repmat(w,S,1); % stack observed wages S times
            e = par.sigma*randn(n*S,1);
            
            con = model.c_star(w,e,par);
            lab = model.l_star(w,e,par);
            
            % c. calculate moments based on simulated data for the value of theta
            mom_sim = moments_fun(w,con,lab);
            
            % d. calculate objective function as squared difference
            distance = mom_data - mom_sim;
            display(mom_sim)
            Q = distance'*distance;
            
        end
    end
end