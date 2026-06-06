function val = qh_interp(t,x1,x2,y1,y2,dy1,dy2,ddy1,ddy2)
x   = [x1,x2];
yin = [  y1,  y2;...
        dy1, dy2;...
       ddy1,ddy2];
[val, ~, ~] = arrayfun(@(t)pqhermite(x,yin,t),t);
end