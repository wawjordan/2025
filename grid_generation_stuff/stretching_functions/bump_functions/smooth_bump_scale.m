function val = smooth_bump_scale(x,xloc,w,wm,wp,scale_factor)
a = xloc - max(wm,w/2);
b = xloc - min(wm,w/2);
c = xloc + min(wp,w/2);
d = xloc + max(wp,w/2);
h = 1 - scale_factor;
val = 1 - h*smooth_transition(x,a,b,c,d);
end