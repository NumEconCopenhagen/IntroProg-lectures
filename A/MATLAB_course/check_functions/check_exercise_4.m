function [success] = check_exercise_4(task,results)

success = 0;

switch task
   
case 'task_1'
    
    % check field
    if isfield(results,'str') == false
        fprintf('there is no str field');
        return;
    end
    
    % check value
    correct_str = 'a = 1, b = 70';
    if isequal(results.str,correct_str) == false
        fprintf('the text is not correct, should be:\n');
        fprintf(correct_str);
        return;
    end
    
    % return with success
    success = 1; return;
    
case 'task_2'
    
    % check for fields
    if isfield(results,'str') == false
        fprintf('there is no str field');
        return;
    end
    
    % check value
    correct_str = 'a = 1\nb = 70';
    if isequal(results.str,sprintf(correct_str)) == false && ...
        isequal(results.str,sprintf('a = 1 \nb = 70')) == false
        fprintf('the text is not correct, should be:\n');
        fprintf(correct_str);
        return;
    end
    
    % return with success
    success = 1; return;
    
case 'task_3'
    
    % check for fields
    if isfield(results,'str') == false
        fprintf('there is no str field');
        return;
    end
    
    % check value
    correct_str = 'a =  1\nb = 70';
    if isequal(results.str,sprintf(correct_str)) == false && ...
        isequal(results.str,sprintf('a =  1 \nb = 70')) == false
        fprintf('the text is not correct, should be:\n');
        fprintf(correct_str);
        return;
    end
    
    % return with success
    success = 1; return;
    
case 'task_4'
    
    % check for fields
    if isfield(results,'str') == false
        fprintf('there is no str field');
        return;
    end
    
    % check value    
    txt_str = '%8.4f\n%8.4f\n%8.4f';
    correct_str = sprintf(txt_str,pi,-(pi+10),0);

    if isequal(results.str,correct_str) == false
        fprintf('the text is not correct, should be:\n');
        fprintf(correct_str);
        return;
    end
    
    % return with success
    success = 1; return;
    
case 'task_5'
    
    % check for fields
    if isfield(results,'str') == false
        fprintf('there is no str field');
        return;
    end
    
    % check value    
    str1 = 'this';
    str2 = 'is';
    str3 = 'simple';
    txt_str = '%s %s %s';
    correct_str = sprintf(txt_str,str1,str2,str3);
    
    if isequal(results.str,correct_str) == false
        fprintf('the text is not correct, should be:\n');
        fprintf(correct_str);
        return;
    end
    
    % return with success
    success = 1; return;
    
otherwise


end

end