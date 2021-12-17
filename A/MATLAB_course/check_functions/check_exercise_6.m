function [success] = check_exercise_6(task,results)

success = 0;
switch task
   
case 'task_1'
    
    % check file
    if exist('exercise_6_task_1.mat','file') == false
        fprintf('exercise_6_task_1.mat has not been created')
        return;
    else
        load('exercise_6_task_1.mat')
    end
    
    % check a
    if exist('a','var') == false
        fprintf('variable a is not loaded');
        return;
    end
    if isequal(a,1:10) == false
        fprintf('variable a is not correct');
        return;
    end   
    
    % check b
    if exist('b','var') == false
        fprintf('variable b is not loaded')              
        return;
    end
    if isequal(b,1:5) == false
        fprintf('variable  b is not correct')        
        return;
    end
    
    % return with success    
    success = 1; return;
    
case 'task_2'
    
    % check a
    if isfield(results,'a') == false
        fprintf('a is not a field in results')
        return;
    end    
    if isequal(results.a,1:10) == false
        fprintf('field a is not correct')
        return;
    end
    
    % check b
    if isfield(results,'b') == false
        fprintf('a is not a field in results') 
    	return;
    end
    if isequal(results.b,1:5) == false
        fprintf('field  b is not correct')
        return;
    end    
    
    % return with success    
    success = 1; return;
    
case 'task_3'
    
    % check year
    if isfield(results,'year') == false
        fprintf('year is not a field in results')
        return;
    end
    if isequal(results.year,(1999:2003)') == false
        fprintf('field year is not correct')
        return;
    end
    
    % check GDP
    if isfield(results,'GDP') == false
        fprintf('GDP is not a field in results')
        return;
    end
    if isequal(results.GDP,ones(5,1)) == false
        fprintf('field GDP is not correct')
        return;
    end
    
    % return with success    
    success = 1; return;
    
case 'task_4'   
    
    % check file
    if exist('exercise_6_task_4.xlsx','file') == false
        fprintf('exercise_6_task_4.xlsx does not exist');
        return;
    else
        [num,txt,~] = xlsread('exercise_6_task_4.xlsx','GDP');
    end
    
    % check size
    if isequal(size(txt),[1 2]) == false
        size(txt)
        fprintf('fieldnames are not correct');
        return;
    end
    if isequal(size(num),[5 2]) == false
        fprintf('numeric values are are not correct');
        return;
    end
    
    results.(txt{1}) = num(:,1);
    results.(txt{2}) = num(:,2);
    
    % check year
    if isfield(results,'year') == false
        fprintf('A1 is wrong')
        return;
    end
    if isequal(results.year,(1999:2003)') == false
        fprintf('year is not correct')
        return;
    end
    
    % check GDP
    if isfield(results,'GDP') == false
        fprintf('B1 is wrong')
        return;
    end
    if isequal(results.GDP,ones(5,1)) == false
        fprintf('GDP is not correct')
        return;
    end
    
    % return with success    
    success = 1; return;
    
otherwise
        
    error('Unknown task, ''%s''.',task)

end
end