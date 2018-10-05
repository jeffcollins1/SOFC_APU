function plot_case(data,options,parameter_name,x_cat,y_cat)
Z = data.(parameter_name);
if isempty(y_cat)
    plot(options.(x_cat),Z);
else
    surf(options.(x_cat),options.(y_cat),Z);
end
end