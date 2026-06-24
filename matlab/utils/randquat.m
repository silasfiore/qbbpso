function q = randquat()
%RANDQUAT Summary of this function goes here
%   Detailed explanation goes here

while 1
  x = 2*(rand()-0.5);
  y = 2*(rand()-0.5);

  z = x*x + y*y;



  if ~(z > 1.0)
      break;
  end
end

while 1
  u = 2*(rand()-0.5);
  v = 2*(rand()-0.5);

  w = u*u + v*v;



  if ~(w > 1.0)
      break;
  end
end
s = sqrt((1-z) / w);

q = [x;y;s*u;s*v];

end

