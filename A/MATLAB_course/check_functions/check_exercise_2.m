function [success] = check_exercise_2(task,results)

success = 0;

switch task
   
case  'task_1'
    
    % check for fields
    if (isfield(results,'z')) == false
        fprintf('some fields are missing.\n')
        return;
    end
    
    % check for values
    if (isequal(results.z,-1)) == false
        fprintf('some fields are not correct.\n')
        return;
    end
    
    % return with success
    success = 1; return;

case 'task_2'
 
    % check for fields
    if (isfield(results,'x')) == false
        fprintf('some fields are missing.\n')
        return;
    end
    
    % check for values    
    for i = 1:10
        if i <= 5
            x(i) = 0;
        else
            x(i) = 1;
        end
    end
    if (isequal(results.x,x)) == false
        fprintf('some fields are not correct.\n')
        return;
    end
    
    % return with success
    success = 1; return;
    

case 'task_3'
 
    % check for fields
    if (isfield(results,'x')) == false
        fprintf('some fields are missing.\n')
        return;
    end
    
    % check for values    
    for i = 1:10
        if i <= 5
            x(i) = 0;
        else
            x(i) = 1;
        end
    end
    if (isequal(results.x,x)) == false
        fprintf('some fields are not correct.\n')
        return;
    end
    
    % return with success
    success = 1; return;
    
case 'task_4'

    d = [2,10,10,5,9,2,5,10,8,10];
    e = d > 5;
    f = d(e);

    % check for fields
    if (isfield(results,'e') && isfield(results,'f')) == false
        fprintf('some fields are missing.\n')
        return;
    end
    
    % check for values
    if (isequal(results.e,e) && isequal(results.f,f)) == false
        fprintf('some fields are not correct.\n')
        return;
    end
    
    % return with success
    success = 1; return;
    
otherwise
        
    error('Unknown task, ''%s''.',task)

end

end