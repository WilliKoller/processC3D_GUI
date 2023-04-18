
function [e1,e2,e3] = segmentorientation_r(V1,V3)

for i=1:length(V1)
   e1=V1/sqrt(dot(V1,V1));
   e3=V3/sqrt(dot(V3,V3));
end

for i=1:length(V1)
   e2 = cross(e1,e3);
   e2 = e2/sqrt(dot(e2,e2));
end

for i=1:length(V1)
   e3 = cross(e2,e1);
end
end

