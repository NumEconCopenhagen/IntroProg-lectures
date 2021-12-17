function [success] = check_exercise_3(task,results)

success = 0;

switch task
   
case 'task_1'
    
    % check for fields
    if (isfield(results,'f1') && isfield(results,'f2')) == false
        fprintf('fields f1 and/or f2 are missing in the results struct.\n')
        return;
    end

    % check data types 
    if ~strcmp(class(results.f1),'function_handle') 
        fprintf('f1 is not a function handle.\n'); 
        return;
    elseif ~strcmp(class(results.f2),'function_handle') 
        fprintf('f2 is not a function handle.\n'); 
        return;
    end; 

    if nargin(results.f1)~=1 || nargin(results.f2)~=1 
        fprintf('handles must be functions of a single variable.\n'); 
        return; 
    end; 
    
    % check that functions can be evaluated
    try 
        results.f1(0); 
    catch me 
        fprintf('the function f1 cannot be evaluated.\n'); 
        return;
    end; 
    try 
        results.f2(0); 
    catch me 
        fprintf('the function f2 cannot be evaluated.\n'); 
        return;
    end; 
    
    % check for values over a grid 
    xx = linspace(-100,100,100); 
    f1 = @(x) 3 + 2*x - 4*x^2; 
    f2 = @(x) exp(f1(x) - 1); 
    for i=1:numel(xx)
        if f1(xx(i))~=results.f1(xx(i))
            fprintf('f1 function is incorrect.\n'); 
            return; 
        elseif f2(xx(i))~=results.f2(xx(i))
            fprintf('f2 function is incorrect.\n'); 
            return; 
        end; 
    end; 
    
    % return with success
    success = 1; return;
    
case  'task_2'
    
    if isfield(results,'f') == false
        fprintf('the field ''f'' is missing in the results struct.\n')
        return;
    end

    if ~strcmp(class(results.f),'function_handle') 
        fprintf('f is not a function handle.\n'); 
        return;
    end; 
    
    if nargin(results.f)~=2
        fprintf('function handle must be functions of two variables.\n'); 
        return; 
    end; 

    % check that functions can be evaluated
    try 
        results.f(0,0); 
    catch me 
        fprintf('the function f cannot be evaluated.\n'); 
        return;
    end; 

    N = 10; % # of points to test
    xx = 20*(rand(N,1)-.5); % uniform [-10,10]
    yy = 20*(rand(N,1)-.5); 
    a = 1; 
    b = 1; 
    f = @(x,y) (a-x)^2 + b*(y-x^2)^2; 
    for ix=1:N
        for iy=1:N
            err = results.f(xx(ix),yy(iy)) - f(xx(ix),yy(iy)); 
            if abs(err)>1e-08 
                fprintf('incorrect values coming out of f.\n'); 
                return; 
            end; 
        end; 
    end; 

    % Return 
    success = 1; 
    return; 

case 'task_3'

    if exist('fun','file')~=2 
        fprintf('unable to find file fun.m.\n'); 
        fprintf('make sure the file is in the same directory and named ''fun.m''.\n'); 
        return; 
    else 
        fprintf('file fun.m located.\n'); 
    end; 

    try fun(1); 
    catch me 
        fprintf('unable to evaluate fun(1).\n');
        fprintf('function returns this error:\n');
        fprintf('%s\n',me.message); 
        return; 
    end; 

    if nargin('fun')~=1
        fprintf('fun() must take precisely one input.\n'); 
        return; 
    end;

    if nargout('fun')~=1
        fprintf('fun() must return precisely one output.\n'); 
        return; 
    end; 

    fprintf('fun() has correct number of inputs/outputs.\n'); 

    % check values
    N = 100; 
    xx = linspace(-10,10,N); 
    f = @(x) sin(sqrt(exp(x-3))); 
    for i=1:N 
        err = fun(xx(i))-f(xx(i)); 
        if abs(err)>1e-08 
            fprintf('fun() returns incorrect values.\n'); 
            return; 
        end; 
    end; 

    fprintf('fun() returns correct values.\n'); 

    % Return 
    success = 1; 
    return; 
    
case 'task_4'
    
    if exist('rosenbrock','file')~=2 
        fprintf('unable to find file rosenbrock.m.\n'); 
        fprintf('make sure the file is in the same directory and named ''rosenbrock.m''.\n'); 
        return; 
    else 
        fprintf('file rosenbrock.m located.\n'); 
    end; 

    try rosenbrock(1,1); 
    catch me 
        fprintf('unable to evaluate rosenbrock(1,1).\n');
        fprintf('function returns this error:\n');
        fprintf('%s\n',me.message); 
        return; 
    end; 

    if nargin('rosenbrock')~=2
        fprintf('rosenbrock() must take precisely two inputs.\n'); 
        return; 
    end;

    if nargout('rosenbrock')~=1 
        fprintf('rosenbrock() must return precisely one output.\n'); 
        return; 
    end; 

    if ~isfield(results,'ros') 
        fprintf('the struct results has no field named ''ros''.\n'); 
        return; 
    end; 

    if ~strcmp(class(results.ros),'function_handle') 
        fprintf('results.ros is not a function handle.\n'); 
        return;
    end; 

    fprintf('results.ros has correct data type.\n'); 

    if nargin(results.ros)~=1 
        fprintf('function handle, results.ros(), takes %i inputs (should take 1).\n',nargin(results.ros)); 
        return; 
    end; 

    try results.ros(1); 
    catch me 
        fprintf('unable to evaluate results.ros(1).\n'); 
        fprintf('the function handle returns this error:\n'); 
        fprintf('%s\n',me.message); 
        return; 
    end; 

    % testing the function handle, results.ros()
    N = 100; 
    xx = linspace(-10,10,N); 
    a = 1; 
    b = 1; 
    y = 3; 
    f = @(x) (a-x)^2 + b*(y-x^2)^2; 
    for i=1:N
        err = f(xx(i))-results.ros(xx(i)); 
        if abs(err)>1e-08 
            fprintf('results.ros() returns incorrect values (err=%f for i=%i).\n',err,i); 
            return; 
        end; 
    end; 

    fprintf('results.ros() returns correct values.\n'); 

    % testing Rosenbrock function 
    N = 10; % # of points to test
    xx = 20*(rand(N,1)-.5); % uniform [-10,10]
    yy = 20*(rand(N,1)-.5); 
    a = 1; 
    b = 1; 
    f = @(x,y) (a-x)^2 + b*(y-x^2)^2; 
    for ix=1:N 
        for iy=1:N 
            if rosenbrock(xx(ix),yy(iy))~=f(xx(ix),yy(iy)) 
                fprintf('incorrect values coming out of rosenbrock().\n'); 
                return; 
            end; 
        end; 
    end; 

    fprintf('rosenbrock(.,.) returns correct values.\n'); 

    % Return
    success = 1; 
    return; 

case 'task_5'

    if exist('classier','file')~=2
        fprintf('unable to find file classier.m.\n'); 
        fprintf('make sure the file is in the same directory and named ''classier.m''.\n'); 
        return; 
    else 
        fprintf('file classier.m located.\n'); 
    end; 

    try 
        cl = class(classier); 
    catch me 
        fprintf('error in class file:\n'); 
        fprintf('%s\n',me.message); 
        return; 
    end; 

    if ~strcmp(cl,'classier')
        fprintf('classier.m does not appear as aproper class.\n'); 
        return; 
    end; 

    pr = properties(classier); 
    if numel(pr)~=1
        fprintf('classier should have precisely one property.\n'); 
        return; 
    end; 
    if ~strcmp(pr,'A') 
        fprintf('classier does not have the property called ''A''.\n'); 
        return; 
    end; 
    if ~isscalar(classier.A) 
        fprintf('classier.A should be a scalar.\n'); 
        return; 
    end; 
    A = pi; 
    if abs(A-pi)>1e-08
        fprintf('classier.A takes the wrong value.\n'); 
        return; 
    end; 

    me = methods(classier); 
    if numel(me)~=3
        fprintf('classier should have precisely two Static methods.\n'); 
        return; 
    end; 
    if ~strcmp('f1',me{1}) && ~strcmp('f1',me{2}) && ~strcmp('f1',me{3}) 
        fprintf('classier does not have method f1().\n'); 
        return; 
    end; 
    if ~strcmp('f2',me{1}) && ~strcmp('f2',me{2}) && ~strcmp('f2',me{3}) 
        fprintf('classier does not have method f2().\n'); 
        return; 
    end; 

    fprintf('property and methods exist.\n'); 

    try classier.f1(0,0); 
    catch me 
        fprintf('unable to evaluate classier.f1(0,0).\n'); 
        fprintf('error message:\n%s\n',me.message); 
        return; 
    end; 
    try classier.f2(0,0,0); 
    catch me 
        fprintf('unable to evaluate classier.f2(0,0,0).\n'); 
        fprintf('error message:\n%s\n',me.message); 
        return; 
    end; 

    fprintf('no syntax error in f1 or f2.\n'); 

    f1 = @(x,y) max(exp(x),A*y); 
    f2 = @(x,y,z) x + 3*x*f1(y,z); 

    N = 20; 
    xx = linspace(-10,10,N); 
    yy = xx; 
    zz = xx; 
    for ix=1:N 
        for iy=1:N
            err = f1(xx(ix),yy(iy)) - classier.f1(xx(ix),yy(iy)); 
            if abs(err)>1e-08
                fprintf('incorrect value in classier.f1().\n'); 
                return; 
            end; 
            for iz=1:N 
                err = f2(xx(ix),yy(iy),zz(iz)) - classier.f2(xx(ix),yy(iy),zz(iz)); 
                if abs(err)>1e-08 
                    fprintf('incorrect value in classier.f2().\n'); 
                    return; 
                end; 
            end; 
        end; 
    end; 

    fprintf('classifier.f1() is correct.\n'); 
    fprintf('classifier.f2() is correct.\n'); 


    % Return
    success = 1; 
    return; 

otherwise
        
    error('Unknown task, ''%s''.',task)

end

end