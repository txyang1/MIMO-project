function [rate,T,F] = MSEmac2bc(H,W,G)

K = length(H);
WHG = zeros(K,K);
whg = zeros(K,1);
trg = zeros(K,1);
Gamma = zeros(K,K);

for i=1:K
    for k=1:K
        WHG(i,k) = norm(W{i}'*H{i}*G{k}','fro')^2;
    end
end
for k=1:K
    whg(k) = sum(WHG(1:end ~=k,k));
    trg(k) = trace(G{k}'*G{k});
end

for k=1:K
    for j=1:K
        Gamma(k,j) =-norm(W{k}'*H{k}*G{j}','fro')^2;
    end
end
for k=1:K
    Gamma(k,k) = whg(k) + trg(k);
end

tau = zeros(K,1);
for k=1:K
    tau(k) = trace(W{k}*W{k}');
end

a = sqrt(inv(Gamma)*tau);

T = cell(K,1);
F = cell(K,1);
for k=1:K
    T{k} = a(k)*G{k}';
    F{k} = (1/a(k))*W{k}'; 
end

rate = MSErate(H,T);

end

function rate = MSErate(H,T)

K = length(H);
[M,~] = size(H{1});

r = zeros(K,1);
for k=1:K
    C = eye(M);
    for i=1:K
        if i==k
            continue
        end
        C = C + H{k}*T{i}*T{i}'*H{k}';
    end
    r(k) = real(log2(det(eye(M) + C\(H{k}*T{k}*T{k}'*H{k}'))));
end

rate = sum(r);

end