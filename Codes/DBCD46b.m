function [ytest1,ytest]=DBCD46b(aa,aax,bb,x,x0,lanf,c,qq,b)
n1=length(aa(1,1,1,:,1));
n2=length(x(1,1,:,1,1));
n20=length(x0(1,1,:,1,1));
n3=length(x(1,:,1,1,1));
n4=length(x(:,1,1,1,1));
c=c(1:n4,:,1:n1);
x=x(1:n4,:,:,1:n1,:);
x0=x0(1:n4,:,:,1:n1,:);
tic
yy2=aa-aax+bb-lanf*(permute(repmat(c,1,1,1,n2),[1 2 4 3]));
for k=1:n20
    K=sum(exp(-(repmat(x0(:,:,k,:,:),1,1,n2)-x).^2./qq),5);
    yy3(:,:,:,:,k)=yy2.*K;
end 
yy3=permute(yy3,[1 2 5 4 3]);
ytest=sum(yy3,5);
ytest1=sum(yy3,5)+permute(repmat(b,1,1,1,1),[1 2 4 3]);
toc