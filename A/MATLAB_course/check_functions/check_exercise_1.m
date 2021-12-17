function [success] = check_exercise_1(task,results)

success = 0;

switch task
   
case 'task_1'
    
    % check for fields
    if (isfield(results,'a') && isfield(results,'b') && isfield(results,'c')) == false
        fprintf('some fields are missing.\n')
        return;
    end
    
    % check for values
    if (isequal(results.a,2) && isequal(results.b,'another example') && isequal(results.c,[1,2,3])) == false
        fprintf('some fields are not correct.\n')
        return;
    end
    
    % return with success
    success = 1; return;
    
case  'task_2'
    
    % check existence 
    if (isfield(results,'X') && isfield(results,'x') && isfield(results,'k')) == false
        fprintf('some fields are missing.\n'); 
        return; 
    end; 
       
    % check values 
    X = [1,2;3,4]; 
    x = zeros(2,4);
    k = 3:3:30;
    if (isequal(results.X,X) && isequal(results.x,x) && isequal(results.k,k)) == false 
        fprintf('some fields are not correct.\n');
        return;
    end; 
    
    % return 
    success = 1; return; 

case 'task_3'

    % check existence 
    if (isfield(results,'Y') && isfield(results,'Z')) == false
        fprintf('some fields are missing.\n'); 
        return; 
    end; 
    
    % check values 
    X = [1,2;3,4]; 
    if (isequal(results.Y,X*X) && isequal(results.Z,X.*X)) == false
        fprintf('some fields are not correct.\n');
        return; 
    end; 
    
    % return 
    success = 1; 
    return; 
    
case 'task_4'
    
    % check existence 
    if (isfield(results,'X_row2') && isfield(results,'X_col2')) == false
        fprintf('some fields are missing.\n'); 
        return; 
    end; 
    
    % check values 
    X = [1,2;3,4];     
    X_row2 = X(2,:);
    X_col2 = X(:,2);    
    if (isequal(results.X_row2,X_row2) && isequal(results.X_col2,X_col2)) == false
        fprintf('some fields are not correct.\n');
        return; 
    end; 
    
    % return
    success = 1; 
    return; 

case 'task_5'

    % check existence 
    if (isfield(results,'vectors')) == false
        fprintf('some fields are missing.\n'); 
        return; 
    end; 
    
    % check values 
    vector_short = linspace(0,1,10);
    vector_long = linspace(0,1,100);
    vectors = cell(2,1);
    vectors{1} = vector_short;
    vectors{2} = vector_long;
    if (isequal(results.vectors,vectors))== false
        fprintf('some fields are not correct.\n');
        return; 
    end; 
    
    % return
    success = 1; 
    return; 


otherwise
        
    error('unknown task, ''%s''.',task)

end

end