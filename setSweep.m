function [incs,inits,ends]=setSweep(incs,inits,ends, dimsize)
NDims=3;
for i = 1:NDims
    if((incs(i) + 2) <=1)
        incs(i)=incs(i)+2;
        break;
    else
        incs(i) = -1;
    end
end

% Setting inits and ends.
for i = 1:NDims
   if (incs(i) == 1)
       inits(i) = 1;%0;
       ends(i) = dimsize(i);
   else
       inits(i) = dimsize(i);%-1;
       ends(i) = 1;%-1;
   end
end
  