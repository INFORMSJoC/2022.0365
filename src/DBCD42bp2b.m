function ytest=DBCD42bp2b(aa,aax,bb,r,th,x,x0,lanf,c,qq,p,b)
n1=length(x(1,1,1,:,1));
n2=length(x(1,1,:,1,1));
n20=length(x0(1,1,:,1,1));
n3=length(x(1,:,1,1,1));
n4=length(aa(:,1,1,1,1));
c=c(1:n4,:,1:n1);
x=x(1:n4,:,:,1,:);
x0=x0(1:n4,:,:,1,:);
p=p(1:n4,1:n1);
r=r(1:n4,:);
tic
yy2=aa-aax+bb-repmat(r,1,n3,n2,n1)-permute(repmat(p,1,1,n3,n2),[1 3 4 2]).*permute(repmat(th,1,n4,n1,n2),[2 1 4 3])-lanf*(permute(repmat(c,1,1,1,n2),[1 2 4 3]));
for k=1:n20
    K=sum(exp(-(repmat(x0(:,:,k,:,:),1,1,n2,n1)-repmat(x,1,1,1,n1)).^2./qq),5);
    yy3(:,:,:,:,k)=yy2.*K;
end 
yy3=permute(yy3,[1 2 5 4 3]);
ytest=sum(yy3,5)+permute(repmat(b,1,1,1,1),[1 2 4 3]);
toc