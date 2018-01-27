function sim = label_gs(A,B,L,c)
% graph kernel computation for attributes graphs
% input:
%   A: first graph
%   B: second graph
%   L: label matrix
%   c: decay factor
% output:
%   sim: similarity between A and B

    n1 = size(A,1);
    n2 = size(B,1);
    q = {ones(n1,1)/n1,ones(n2,1)/n2};
    p = {ones(n1,1)/n1,ones(n2,1)/n2};
    p1 = p{1};
    p2 = p{2};
    q1 = q{1};
    q2 = q{2};
    
    X = kron(A,B);
    qx = kron(q1,q2);
    px = kron(p1,p2);
    sim = qx' * inv(eye(n1 * n2) - c * L * X) * L * px;
end