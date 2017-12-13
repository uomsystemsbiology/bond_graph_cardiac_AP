function d = gcd_array(A)
% Check if all elements are integers
if any(floor(A)-A)
    error('The input matrix has non-integer elements.');
end

d = min(A);
while true
  r = mod(A,d);
  if ~any(r)
    break
  end
  r(r == 0) = inf;
  d = min(r);
end

end