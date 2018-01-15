function prd = StartEnding(x)
% x should be a vector, or in column
[s1, s2] = size(x);
if s2>s1
    x=x';
    s2 = s1;
end
prd = cell(s2, 1);
x(end,:)=0;
tmp = diff([zeros(1, s2); x],1,1);
for k = 1:s2
    prd{k} = [find(tmp(:,k)>0), find(tmp(:,k)<0)];
end
if s2 ==1
    prd = cell2mat(prd);
end
