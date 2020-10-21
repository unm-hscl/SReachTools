function plot(obj, ax)

if isempty(ax)
    ax = gca;
end

plot(obj.ConstraintTube);
plot(obj.TargetTube);

end
