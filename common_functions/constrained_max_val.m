function [x,y,z] = constrained_max_val(X,Y,Z,P,p)
%finds coordinates x,y coresponding to maximum Z value, subject to
%constraint P = p. Uses vectors X & Y, and matrices Z, P. p is the input condition
%limitation: doesn't do well with local maxima of nearly equal height
[x_p,y_p] = constrain_points(X,Y,P,p);

%sort for most efficient points
n = length(x_p);
z_p = zeros(n,1);
for i = 1:1:n
    z_p(i) = griddata(X,Y,Z',x_p(i),y_p(i));
end
[~,I] = sort(z_p,1,'descend');
x_s = x_p(I);
y_s = y_p(I);
if x_s(1) == x_s(2)%search only in y-direction at x = x_s(1);
    x_p = ones(100,1)*x_s(1);
    y_p = linspace(0,1)'*(y_s(2) - y_s(1)) + y_s(1);
elseif y_s(1) == y_s(2)%search only in x-direction at y = y_s(1);
    x_p = linspace(0,1)'*(x_s(2) - x_s(1)) + x_s(1);
    y_p = ones(100,1)*y_s(1);
else %search in x dimension between x_s(1) and x_s(2), subject to y_i that satisfies P = p
    X2 = linspace(0,1)'*(x_s(2) - x_s(1)) + x_s(1);
    if x_s(1) == X(end)
        j = nnz(X<=x_s(2));
    else
        j = nnz(X<=x_s(1));
    end
    P2 = linspace(0,1)'*(P(j+1,:)-P(j,:)) + ones(100,1)*P(j,:);
    [x_p,y_p] = constrain_points(X2,Y,P2,p);
end
n = length(x_p);
z_p = zeros(n,1);
for i = 1:1:n
    z_p(i) = griddata(X,Y,Z',x_p(i),y_p(i));
end
[z_s,I] = sort(z_p,1,'descend');
x = x_p(I(1));
y = y_p(I(1));
z = z_s(1);
end

function [x_p,y_p] = constrain_points(X,Y,P,p)
%find Y coordinates that satisfy P =p at X(i)
n = length(X);
n2 = length(Y);
x_p = [];
y_p = [];
for i = 1:1:n
    x_i = X(i);
    if any(P(i,:)<=p) && any(P(i,:)>=p)
        rem = P(i,:)-p;
        zero_rem = nonzeros((1:n2-1).*(rem(1:end-1)<=0 & rem(2:end)>=0));%increasing value past p
        if ~isempty(zero_rem)
            for jj = 1:1:length(zero_rem)
                j = zero_rem(jj);
                x_p(end+1,1) = x_i;
                r = -rem(j)/(rem(j+1)-rem(j));%fraction in Y direction for P = p at coordinate X(i)
                y_p(end+1,1) = Y(j) + r*(Y(j+1)-Y(j));%y-coordinate that satisfies P = p @ X(i) 
            end
        end
        zero_rem = nonzeros((1:n2-1)'.*(rem(1:end-1)>=0 & rem(2:end)<=0));%decreasing value past p
        if ~isempty(zero_rem)
            for jj = 1:1:length(zero_rem)
                j = zero_rem(jj);
                x_p(end+1,1) = x_i;
                r = -rem(j)/(rem(j+1)-rem(j));%fraction in Y direction for P = p at coordinate X(i)
                y_p(end+1,1) = Y(j) + r*(Y(j+1)-Y(j));%y-coordinate that satisfies P = p @ X(i) 
            end
        end
    end
end
end



				
					
				

								

				
					
					
					
				

				
					
				

				

				
					