function [outputArg1,outputArg2] = motionPlot(x,y,z)
%%
fig = gcf;
if fig.Number == 1
    morph = fig;
else
    num = fig.Number+1;
    morph = figure(1);
end

ax = axes;

for k = 1:(length(t)-1)

    if ishandle(num) == false
      break;
    end 
    surf(ax,x{k},y{k},z{k});
    zlim(ax,[-5 5])
    xlim(ax,[-5 5])
    ylim(ax,[-5 5])
    drawnow
    pause(0.2)
end


end

