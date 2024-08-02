function [ytrain1,b]=DBCD32bp2bb(aa,aax,bb,r,th,x,y,lanf,c,qq,p)
%Training set, matrix computation������ģ�ͳ�ʼֵ���㣬�ಽԤ�⣬����ƫ��,b����ʱ�������ڵ���
n1=length(aa(1,1,1,:));
n2=length(x(1,1,:,1,1));
n3=length(x(1,:,1,1,1));
n4=length(aa(:,1,1,1));
c=c(1:n4,:,1:n1);%����ǰ������ά���޸�
x=x(1:n4,:,:,1,:);%����ǰ������ά���޸�
p=p(1:n4,1:n1);%����ǰ������ά���޸�
r=r(1:n4,:);%����ǰ������ά���޸�
tic
yy2=aa-aax+bb-repmat(r,1,n3,n2,n1)-permute(repmat(p,1,1,n3,n2),[1 3 4 2]).*permute(repmat(th,1,n4,n1,n2),[2 1 4 3])-lanf*(permute(repmat(c,1,1,1,n2),[1 2 4 3]));
for k=1:n2
    K=sum(exp(-(repmat(x(:,:,k,:,:),1,1,n2,n1)-repmat(x,1,1,1,n1)).^2./qq),5);
    yy3(:,:,:,:,k)=yy2.*K;
end 
ytrain=sum(yy3,5);
for i=1:n4
    for j=1:n3
        for tt=1:n1
            t=n1-tt+1;
            for k=1:n2
                if(t==n1)
                    y1=0;
                else
                    y1=y(i,j,k,t+1:n1)-ytrain1(i,j,k,t+1:n1);
                end
                b0(k)=y(i,j,k,t)+sum(y1)-ytrain(i,j,k,t);
            end
            b(i,j,t)=sum(b0)/n2;
            ytrain1(i,j,:,t)=ytrain(i,j,:,t)+permute(repmat(b(i,j,t),1,1,1,n2),[1 2 4 3]);
        end
    end
end
toc