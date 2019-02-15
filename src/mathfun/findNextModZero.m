function d_closest = findNextModZero(n, d, mode)
% FINDNEXTMODZERO finds the next modulo zero denominator
% 
%              Given a nominator denominator pair (n,d) the function finds
%              the next denominator d_closest such that mod(n,d_closest)=0.
%              The function returns -1 if no denominator with mod = 0 was
%              found.
%               
%
% SYNOPSIS    [d_closest] = findNextModZero(n, d, mode)
%
% INPUT       n   : nominator
%             d   : denominator
%             mode: -1 : next smaller
%                    0 : absolute next
%                    1 : next bigger
% 
% OUTPUT     d_closest: closest denominator with mod(n,d)=0
%          
%                           
% DEPENDENCES   findNextModZero uses { 
%                                       }
%
%               findNextModZero is used by {                                     
%                                       }
%
%
% Matthias Machacek 05/27/04

if mod(n,d) == 0
    d_closest = d;
    return
end

if mode == 1 | mode == 0
    d_bigger =d;
    while mod(n,d_bigger) ~= 0 & d_bigger <= n
        d_bigger = d_bigger + 1;
    end
    diff_bigger = d_bigger - d;
end

if mode == -1 | mode == 0
    d_smaller =d;
    while mod(n,d_smaller) ~= 0 & d_smaller > 0 
        d_smaller = d_smaller - 1;
    end
    diff_smaller = d - d_smaller;
end


if mode == 0
    if diff_bigger > diff_smaller & d_smaller > 0 
        d_closest = d_smaller;
    elseif d_bigger <= n
        d_closest = d_bigger;
    else
        d_closest = -1;
    end
elseif mode == 1 & d_bigger <= n
    d_closest = d_bigger;
elseif d_smaller > 0
    d_closest = d_smaller;  
end


