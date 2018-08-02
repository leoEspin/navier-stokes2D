function create_movie(plot_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%script for creating an animation of oxygen concentration or averaged 
%oxygen concentration, depending on the string argument passed to 
%the function (either 'average' or any other string).
%requires the files input, average.dat, and the script asymptoticC.m
%plus all the solute files 'sltXX.dat' for making the movie frames.
%NOTE: if an error indicating matrices with non-matching sizes, that
%happens because of that weird fortran file storage error. since 
%the solution is time periodic, the files slt1.dat and slt26.dat are 
%exactly equal. erase slt1, and rename slt26 as slt1. run again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n,re,del,~,~,sig,lamb,aratio,kap,~]=parameters;
dy=1/m;
y=1:(3*m+4);
y(1:2*m+2)=y(1:2*m+2)-1;
y(1:2*m+2)=dy*y(1:2*m+2)-dy/2 -3;
y(2*m+3:3*m+4)=y(2*m+3:3*m+4)-2*m-3;
y(2*m+3:3*m+4)=dy*y(2*m+3:3*m+4)-dy/2 -1;
dx=1/n;
x=(1:n+2)-1;
x=dx*x-dx/2;
h=load('frame_times.dat');
times=h(:,2);
times=times-2*pi*floor(0.5*times/pi);
times=times/pi;
switch plot_type
    case 'average'
        h=figure('Visible','off');
    otherwise
        scrsz = get(0,'ScreenSize');
        h=figure('Position',[1 1 5*scrsz(3)/8 5*scrsz(4)/8],'Visible','off');
end
hax=axes;
set(hax,'Fontsize',14)
switch plot_type
    case 'average'
        aviname=['average R=' num2str(re) ' D=' num2str(del) ' S=' ...
            num2str(sig) ' L=' num2str(lamb) ' d=' num2str(aratio) ' k='...
             num2str(kap) '.avi'];
        aviobj=avifile(aviname,'fps',4,'compression','none'); %creates AVI file
        av=fopen('average.dat');
        v=fread(av,[n+2,1],'float64');
        fclose(av); 
    otherwise
        aviname=['oxygen R=' num2str(re) ' D=' num2str(del) ' S=' ...
            num2str(sig) ' L=' num2str(lamb) ' d=' num2str(aratio) ' k='...
             num2str(kap) '.avi'];
        aviobj=avifile(aviname,'fps',8,'compression','none'); %creates AVI file
end
i=1;
cond=1;
while cond
    name=['slt' num2str(i) '.dat'];
    cond = exist(name,'file');
    if (cond == 0)
        disp(['file slt' num2str(i) '.dat not found'])
        break
    end
    file=fopen(name);
    s=fread(file,[3*m+4,n+2],'float64');
    fclose(file);
    switch plot_type
        case 'average'
            average=0.5*dy*(2*sum(s(2*m+5:3*m+2,:))+s(2*m+4,:)+s(3*m+3,:));
            average=average+dy*(3*s(3*m+3,:)+s(3*m+4,:))/8;
            average=average+dy*(3*s(2*m+4,:)+s(2*m+3,:))/8;                        
            plot(hax,x,average,'r');
            hold on
            plot(x,v,'r-.')
            %asymptoticC
            %axis(hax,[0.9 1 0.95 1]);
            axis(hax,[0 1 0 1]);
            xlabel(hax,'\xi','fontsize',20)
            ylabel(hax,'\eta -average','fontsize',20)      
            title(hax,['time=' num2str(times(i))],'fontsize',20);
            hold off
        otherwise
            pcolor(hax,x,y,s);shading flat  %pseudocolor plot
            caxis([0 1])
            xlabel(hax,'\xi','fontsize',20)
            ylabel(hax,'\eta','fontsize',20)     
            title(hax,['time=' num2str(times(i))],'fontsize',20);
    end
    aviobj=addframe(aviobj,h); %adds frames to the AVI file
    %     name=['frame' num2str(i) '.png'];
    %     saveas(hax,name,'png')
    i=i+1;
end
aviobj=close(aviobj); %closes the AVI file
close(h); %closes the handle to invisible figure

