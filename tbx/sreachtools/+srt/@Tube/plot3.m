function ph = plot3(obj, ax)

if isempty(ax)
    ax = gca
end

next = ax.NextPlot;

colors = get(ax, 'ColorOrder');
color_index  = get(ax, 'ColorOrderIndex');

if color_index > size(colors, 1)
  color_index = 1;
end

next_color = colors(color_index, :);

for k = 1:length(obj.tube_)
    ax.NextPlot = 'add';

    if k == 1
        ph = plot(obj.tube_(k), 'Alpha', 0.5, 'Color', next_color);
    else
        plot(obj.tube_(k), 'Alpha', 0.1, 'LineStyle', ':', 'Color', next_color);
    end

end

ax.NextPlot = next;

next_color_index = color_index + 1;

if next_color_index > size(colors, 1)
  next_color_index = 1;
end

set(ax, 'ColorOrderIndex', next_color_index);

end
