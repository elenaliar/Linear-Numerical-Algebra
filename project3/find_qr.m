function [Q, R] = find_qr(A, m, n)
    Q = eye(m);
    R = A;
    for i=1:min(m - 1, n)
        [Q, R] = qr_step(Q, R, i, m);
    end
end

function [Q, R] = qr_step(Q, R, i, m)
    v = get_v_vector(R, i);
    H = eye(m);
    H(i:end, i:end) = house_holder(v);
    R = H * R;
    Q = Q * H;
end

function v = get_v_vector(R, i)
    x = R(i:end, i);
    e1 = zeros(size(x));
    e1(1) = 1;
    v = x - norm(x) * e1;
end

function H = house_holder(v)
    I = eye(length(v));
    H = I - 2 * (v * v') / (v' * v);
end