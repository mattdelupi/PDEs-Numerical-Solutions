function W = PesiDer(xs, xc, p)
    N = length(xs);
    xs = xs(:).';
    csi = xs - xc;
    M = zeros(N);
    I = eye(N);

    for j = 1:N
        M(j, :) = csi.^(j-1) ./ factorial(j-1);
    end

    W = M\I;
    W = W.';

    if nargin == 3
        W = W(p+1, :);
    end
end