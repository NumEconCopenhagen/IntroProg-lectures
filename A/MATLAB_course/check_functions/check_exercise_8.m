function [success] = check_exercise_8(task,results)

success = 0;
switch task
   
case 'task_1'
        
    % check field
    if isfield(results,'a') == false
        fprintf('there is no a field');
        return;
    end
    if isfield(results,'b') == false
        fprintf('there is b fval field');
        return;
    end
    if isfield(results,'c') == false
        fprintf('there is b fval field');
        return;
    end
    
    seed = 2017;
    rng(seed);
    a = rand(100,1);
    b = randn(100,1);
    vec = [1 8 9];
    c = vec(randi(3,100,1));
    
    % check value
    if isequal(a,results.a) == false
        fprintf('a is not correct');
        return;
    end
    if isequal(b,results.b) == false
        fprintf('b is not correct');
        return;
    end
    if isequal(c,results.c) == false
        if ismatrix(results.c); 
            if size(results.c,2)>1; 
                fprintf('c is not correct (it is %i*%i; should be 100*1.)',size(results.c)); 
                return ;
            end; 
        end; 
        fprintf('c is not correct');
        return;
    end
    
    success = 1; return;
    
otherwise
        
    error('Unknown task, ''%s''.',task)

end
end