function t = centri_param(points,mu)
% points(1:dim,1:N)
% t = [0,sqrt( sum( (points(:,2:end) - points(:,1:end-1)).^2, 1) ).^mu].';
% t = cumsum(t);
% t = t/t(end);
t = [ 0; ...
      cumsum( ...
        sqrt( sum( abs(points(:,2:end) - points(:,1:end-1)).^2, 1) ).^mu ).'];
% t = t/t(end);
end

% function t = centripetal_param(points,N,mu)
% % points(1:dim,1:N)
% t = zeros(N,1);
% for i = 2:N
%     d = sqrt( sum( (points(:,i) - points(:,i-1)).^2 ) );
%     t(i) = t(i-1) + d^mu;
% end
% 
% end