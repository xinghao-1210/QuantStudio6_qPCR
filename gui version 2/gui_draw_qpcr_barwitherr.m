function [hBar,hErr,mean_used,sd_used] = gui_draw_qpcr_barwitherr(package,graph_param,display_Ct,normed_ecad,lnwdth)
forfigure=true;
if forfigure
    fontweight='normal';
else
    fontweight='bold';
end

askfontsize=false;
fontsize=24; %default if don't ask
asklogscale=true;
logscale=false; %if don't ask, will still be false
%normed_ecad=false; %if you don't want the extra label

sel_tar=graph_param.tars2plot;
%sel_tar=select_from_list_GUI('Target Selection',package.TargetNames,'Please select targets to display (or type "all")');
temp_mean=package.mean(:,sel_tar);
if graph_param.geo
    temp_sd=package.sd(:,sel_tar,:);
else
    temp_sd=package.sd(:,sel_tar);
end

sel_samp=graph_param.grps2plot;
%sel_samp=select_from_list_GUI('Sample Selection',package.SampleNames,'Please select samples to display (or type "all")');
mean_used=temp_mean(sel_samp,:);
mean_used=mean_used';
if graph_param.geo
    sd_used_untransposed=temp_sd(sel_samp,:,:);
    sd_used_low_trans=sd_used_untransposed(:,:,1)';
    sd_used_high_trans=sd_used_untransposed(:,:,2)';
    sd_used(:,:,1)=sd_used_low_trans;
    sd_used(:,:,2)=sd_used_high_trans;
else
    sd_used=temp_sd(sel_samp,:);
    sd_used=sd_used';
end


if askfontsize
   try 
       fontsize=select_from_list_GUI('Font Size',package.TargetNames,'Please type in font size:');
   catch
   end
end


[hBar,hErr]=barwitherr(sd_used,mean_used);colormap('gray')
if size(mean_used,1)>1
    set(gca,'XTickLabel',package.TargetNames(sel_tar),'FontName','Arial','Fontsize',fontsize,'FontWeight',fontweight)
    if graph_param.legend
        legend(package.SampleNames(sel_samp));legend('boxoff')
    end
else
    set(gca,'XTickLabel',package.SampleNames(sel_samp),'FontName','Arial','Fontsize',fontsize,'FontWeight',fontweight)
    title(package.TargetNames(sel_tar),'Fontsize',fontsize+4,'FontWeight','bold')
end
set(gcf,'color',[1,1,1])

yname=['Relative Expression'];
if normed_ecad
    yname=[yname,' norm. to E-cad']; 
end
if display_Ct
    yname=[yname,' (log2 scale)'];
end
ylabel(yname,'FontWeight','bold')

%tidying up the graph
%{
if asklogscale
    temp_bttn = questdlg('Do you want to display y-axis as linear or log?','Y-Axis','Linear','Log','Linear');
    if strcmp(temp_bttn,'Linear')
        logscale=false;
    else
        logscale=true;
    end
end
%}
logscale=graph_param.log;
if logscale
    set(gca,'yscale','log')
else
    yl=ylim;
    ylim([0,yl(2)]);
end


if forfigure
    set(hBar,'BarWidth',0.7)
    lnwdth=2;
    set(gcf,'Position',[300 300 850 600])
end


set(hBar,'LineWidth',lnwdth)
set(hErr,'LineWidth',lnwdth)
set(gca,'LineWidth',lnwdth)

set(gca,'box','off')

end

function whichones = select_from_list_GUI(gui_title,list,questionprompt)

running_list=['1. ',list{1}];
for a=2:max(size(list))
    running_list=[running_list,sprintf('\n'),num2str(a),'. ',list{a}];
end
temp=inputdlg([questionprompt,':  ',sprintf('\n'),running_list],gui_title);
tempstr=temp{1};
if strcmpi(tempstr,'all')
    whichones=[1:max(size(list))];
else
    whichones=str2num(tempstr);
end

end

%{

[whichone,ok]=listdlg('ListString',list,'PromptString',questionprompt,'Name',gui_title,'SelectionMode','multiple');

%}