function [success] = check_exercise_(task,results)

success = 0;
switch task
   
case 'task_1'
    
    x = linspace(0,2*pi,100);
    y = sin(x);
    figure();
    plot(x,y);
    
    success = -1; return;
    
case 'task_2'
    
    x = linspace(0,2*pi,100);
    y1 = sin(x);
    y2 = cos(x);
    fig = figure();
    plot(x,y1)
    hold('on')
    plot(x,y2)
    hold('off')

    success = -1; return;
    
case 'task_3'    
    
    x = linspace(-1,1,100)';
    y = [sqrt(1-x.^2);-sqrt(1-x.^2)];
    x = [x;x];
    fig = figure();
    scatter(x,y)

    success = -1; return;
    
case 'task_4'   

    x = [1 2 1 1 2 2 2 3 3 1 2 2 4 2 3 4 1 1 2 3 2 1 1 1 2];
    fig = figure();
    histogram(x)

    success = -1; return;
    
case 'task_5'    

    x = linspace(0,1,100);
    y = linspace(0,1,100);
    [x_mat, y_mat] = meshgrid(x,y);
    z = (x_mat-0.25).^2 + (y_mat - 0.5).^2;
    fig = figure();
    mesh(x,y,z)
    fig = figure();
    surf(x,y,z)
    
    success = -1; return;
    
case 'task_6'    

    x1 = linspace(0,2*pi,100);
    y1 = sin(x1);

    x2 = linspace(-1,1,100)';
    y2 = [sqrt(1-x2.^2);-sqrt(1-x2.^2)];
    x2 = [x2+2;x2+2];

    fig = figure();
    plot(x1,y1,...
         '--','linewidth',2,'color','red',...
         'DisplayName','sin');
    hold('on')
    scatter(x2,y2,x2.^2,'blue','filled',...
            'DisplayName','circle');
    hold('off');

    xlim([-0.5, 7]);
    ylim([-1.5, 1.5]);

    title('TITLE','fontsize',16);
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',8);
    box('on');
    grid on;

    legend('Location','northeast');

    success = -1; return;
    
case 'task_7'    

    if exist('exercise_5_task_7.png','file') == false
        fprintf('png files does not exist')
        return;
    end
    
    if exist('exercise_5_task_7.pdf','file') == false
        fprintf('pdf files does not exist')
        return;
    end
    
    success = 1; return;
    
otherwise
        
    error('Unknown task, ''%s''.',task)

end
end