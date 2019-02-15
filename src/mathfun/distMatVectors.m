function dm3=distMatVectors(A,B)
%calculates a distance matrix filled with the differences between each element of A and each element of B
%
%
%JD 10/02

if isempty(A) | isempty (B)
    dm3=[];
else

    LA=length(A);
    LB=length(B);
    

    [MA,MB]=meshgrid(1:LA,1:LB);
    
    AMA=A(MA);
    BMB=B(MB);
    if size(AMA)~=size(BMB)
        AMA=AMA';
    end
    dm3=BMB-AMA;
end
  