function [v] = SAD(v1, v2) 

    Num = dot(v1,v2);
    Den1 = norm(v1);
    Den2 = norm(v2);

    Den = Den1*Den2 + eps;
    v = acos( Num / Den );
    
end