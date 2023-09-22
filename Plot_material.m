function [y]=Plot_material(Length,Iter)
    f=fopen('Elastic.dat','r');
    Elastic=[];
    while ~feof(f)
        line=fgetl(f);
        a=str2double(line);
        Elastic=[Elastic;a];
    end
    Elastic=Elastic/1e9;
    fclose(f);
    
    f=fopen('Object.txt','r');
    Object=[];
    while ~feof(f)
        line=fgetl(f);
        a=str2double(line);
        Object=[Object;a];
    end
    Object=Object/1e9;
    fclose(f);
 
    figure(1);
    hold on;
    dis=Length/size(Elastic,1);
    for i=1:size(Elastic,1)
        a1=[dis*(i-1),dis*(i)];
        a2=[Elastic(i),Elastic(i)];
        plot(a1,a2,'r','linewidth',1.5);
    end
    for i=1:size(Elastic,1)-1
        b1=[dis*(i),dis*(i)];
        b2=[Elastic(i),Elastic(i+1)];
        h=plot(b1,b2,'r','linewidth',1.5);
    end
        
    for i=1:size(Object,1)
        a1=[dis*(i-1),dis*(i)];
        a2=[Object(i),Object(i)];
        plot(a1,a2,'k','linewidth',1);
    end
    for i=1:size(Object,1)-1
        b1=[dis*(i),dis*(i)];
        b2=[Object(i),Object(i+1)];
        f=plot(b1,b2,'k','linewidth',1);
    end
    set(gca,'fontname','Times New Roman','fontsize',16,'fontweigh','bold');
    legend([f,h],{'目标分布','反演结果'},'fontname','宋体','fontsize',15,'location','northwest');
    xlabel('\fontname{宋体}距开挖端距离\fontname{Times New Roman}(m)','fontsize',16.5);
    ylabel('\fontname{宋体}弹性模量\fontname{Times New Roman}(GPa)','fontsize',16.5);
    if Iter==0
        name=['plot/Origin','.png'];
    else
        name=['plot/Iter',num2str(Iter),'.png'];
    end
    axis square;
    set(gca,'box','on');
    print(1,'-dpng','-r300',name);
    close(1);
    y=0;
end