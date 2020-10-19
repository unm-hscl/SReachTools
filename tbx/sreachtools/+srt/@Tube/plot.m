function ph = plot(obj, ax)

arguments
    obj
    ax (1, 1) = gca
end

next = ax.NextPlot;

for k = 1:length(obj.tube_)
    ax.NextPlot = 'add';

    if k == 1
        ph = plot(obj.tube_(k), 'alpha', 0.5, 'color', 'y');
    else
        plot(obj.tube_(k), 'alpha', 0.08, 'LineStyle', ':', 'color', 'y');
    end

end

ax.NextPlot = next;

end
