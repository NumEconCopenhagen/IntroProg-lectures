classdef appendix
    methods (Static) 
        
    function [] = test_debug()
        
        x = 2;
        keyboard;
        
    end
    function [] = allocation()
        
        N = 3000000;
        
        t1 = tic;
        for i = 1:N
            A(i) = i;
        end
        fprintf('without pre-allocation: %5.3f secs\n',toc(t1));

        t1 = tic;
        B = zeros(N,1);
        for i = 1:N
            B(i) = i;
        end
        fprintf('with pre-allocation:    %5.3f secs\n',toc(t1));
        
    end        
    function [] = loop_order()
        
        N = 5000;
        
        t1 = tic;
        A = zeros(N,N);
        for i = 1:N
            for j = 1:N
                A(i,j) = i+j;
            end
        end
        fprintf('wrong loop order:   %5.3f secs\n',toc(t1));

        t1 = tic;
        B = zeros(N,N);
        for j = 1:N
            for i = 1:N
                B(i,j) = i+j;
            end
        end
        fprintf('correct loop order: %5.3f secs\n',toc(t1));
               
    end        
    function [] = vectorization()
        
        N = 2000;                
        f = @(x,y) log(x+y)+log(y+x.^2);
        
        t1 = tic;
        A = zeros(N,N);
        for j = 1:N
            for i = 1:N
                A(i,j) = f(i,j);
            end
        end
        fprintf('non-vectorized double loop: %5.3f secs\n',toc(t1));

        t1 = tic;
        B_vec = zeros(N,N);
        i_vec = 1:N;
        for j = 1:N
            B_vec(:,j) = f(i_vec,j);           
        end
        fprintf('vectorized, inner loop:     %5.3f secs\n',toc(t1));
        assert(isequal(A,B_vec));
        
        t1 = tic;
        i_vec = 1:N;
        j_vec = 1:N;
        B_vec_full = bsxfun(f,i_vec,j_vec')'; 

        fprintf('vectorized, both loops:     %5.3f secs\n',toc(t1));
        assert(isequal(A,B_vec_full))
        
    end
    
    end
end
