function prod = prod_quat(p,q)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

p_w = p(1);
p_vec = p(2:4);

q_w = q(1);
q_vec = q(2:4);

prod = [
    p_w*q_w - p_vec'*q_vec
    p_w*q_vec + q_w*p_vec + skew(p_vec)*q_vec
    ];

end

