function [] = check_answer(exercise,task,results)
    
    % a. find task number
    i = strfind(task,'_');
    results.task_number = task(i+1:end);
    fprintf('TASK %s AT HAND:\n',results.task_number);

    % b. call check function
    addpath('check_functions');
    name_str = ['check_' exercise]; 
    path = which(name_str); 
    assert(~isempty(path),sprintf('Unable to find solutions for exercise %s, task %s.',exercise,task));
    success = feval(['check_' exercise],task,results);
    
    if success == 1
        fprintf('SOLUTION IS CORRECT!\n')
    elseif success == 0
        error('SOLUTION IS INCORRECT!')
    end
end