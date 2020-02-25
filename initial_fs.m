function [incs_,inits_,ends_]=initial_fs(n)
%Initializes
inits_=zeros(n,1);
ends_=zeros(n,1);
incs_=zeros(n,1);
for i = 1:n
    incs_(i) = 1;
    inits_(i) = 0;
    ends_(i) = 1;
end