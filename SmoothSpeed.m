function v = SmoothSpeed(x,t,lag)
[s1, s2] = size(x);
if s1<s2
    x = x';
    s2 = s1;
end
for k = 1:s2
    x(:,k) = conv(x(:,k), ones(lag,1)/lag, 'same');
end
if length(t)>1
samplingrate = ceil(1/diff(t(1:2)));
v = bsxfun(@rdivide, -[zeros(lag,s2);x]+[x;zeros(lag,s2)], ...
([t;t(end) + [1:lag]'/samplingrate] - [t(1) + ([1:lag]'-lag-1)/samplingrate;t]));
else
v = -[repmat(x(1,:),lag,1);x]+[x;repmat(x(end,:),lag,1)];
end
v = v([(lag+1):end]-floor(lag/2),:);% middle point speed.