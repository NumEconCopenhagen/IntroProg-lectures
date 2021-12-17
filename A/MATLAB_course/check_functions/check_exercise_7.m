function [success] = check_exercise_7(task,results)

success = 0;
switch task
   
case 'task_1'
    
    % check field
    if isfield(results,'x') == false
        fprintf('there is no x field');
        return;
    end
    if isfield(results,'fval') == false
        fprintf('there is no fval field');
        return;
    end
    
    % check value
    if max(abs(results.x)) > 1e-3 || abs(results.fval) > 1e-6
        fprintf('solution is not correct');
        return;
    end
    
    success = 1; return;
    
case 'task_2'
    
    % check field
    if isfield(results,'x') == false
        fprintf('there is no x field');
        return;
    end
    if isfield(results,'fval') == false
        fprintf('there is no fval field');
        return;
    end
    
    % check value
    if max(abs(results.x)) > 1e-4 || abs(results.fval) > 1e-6
        fprintf('solution is not correct');
        return;
    end
    
    success = 1; return;

case 'task_3'
    
    % check field
    if isfield(results,'x') == false
        fprintf('there is no x field');
        return;
    end
    if isfield(results,'fval') == false
        fprintf('there is no fval field');
        return;
    end
    
    % check value
    if max(abs(results.x-[0.1,0.0])) > 1e-3 || abs(results.fval-0.01) > 1e-6
        fprintf('solution is not correct');
        return;
    end
    
    success = 1; return;
    
case 'task_4'
    
    % check field
    if isfield(results,'x') == false
        fprintf('there is no x field');
        return;
    end
    if isfield(results,'fval') == false
        fprintf('there is no fval field');
        return;
    end
    
    % check value
    if max(abs(results.x)) > 1e-4 || abs(results.fval) > 1e-6
        fprintf('solution is not correct');
        return;
    end
    
    success = 1; return;
    
otherwise
        
    error('Unknown task, ''%s''.',task)

end
end