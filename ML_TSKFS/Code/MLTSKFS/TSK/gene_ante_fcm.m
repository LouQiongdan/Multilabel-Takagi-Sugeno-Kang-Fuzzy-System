function [v,b] = gene_ante_fcm(data,options)
k = options.k;
h = options.h;
[n_examples, d] = size(data);
[v,U,~] = fcm(data,k,[2,NaN,1.0e-6,0]);
for i=1:k
    v1 = repmat(v(i,:),n_examples,1);
    u = U(i,:);
    uu = repmat(u',1,d);
    b(i,:) = sum((data-v1).^2.*uu,1)./sum(uu)./1;
end
b = b*h+eps;

end


