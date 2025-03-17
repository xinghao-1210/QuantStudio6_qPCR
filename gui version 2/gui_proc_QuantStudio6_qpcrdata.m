function [package] = gui_proc_QuantStudio6_qpcrdata(raw,ctrls_Tar,rem_tarnames,ctrl_Samp,grpnames,grpsamps,graph_param)


%%
%$Hard-code
display_Ct=false;
accepted_controls={'CphA','GAPDH','18S','GUSB','HPRT1'};
NaN_is_num=true;
NaNnum=40;
linewidth=3.5;
%display_Ct=true;

ctrl_Tar=ctrls_Tar{1};

Ct_samp_by_tar=raw.Ct_samp_by_tar;
SampName_list=raw.SampleNames;
TarName_list=raw.TarNames;


%%
%{
if nargin<1 || isempty(filename)
    UI_later=true;
    
    temp=dir;
    ind_xls_files=[];
    for n=1:length(temp)
        tempname=temp(n).name;
        len_tn=length(tempname);
        if len_tn>4
            if strcmpi(tempname((len_tn-3):(len_tn)),'.xls') ||  strcmpi(tempname((len_tn-4):(len_tn)),'.xlsx')
                ind_xls_files=[ind_xls_files,n];
            end
        end
    end
    
    if length(ind_xls_files)==1
        dirn=pwd;
        dirn=[pwd,'\'];
        fn=temp(ind_xls_files).name;
    else
        [fn,dirn] = uigetfile('*.xls;*.xlsx','All Excel files');
    end
    filename=[dirn,fn];
else
    UI_later=false;
end
filename2save=[filename(1:(length(filename)-4))];

if exist([filename2save '.mat'])==2
    load([filename2save '.mat']);
    savelater=false;
    
    Ct_samp_by_tar=package.raw.Ct_samp_by_tar;
    SampName_list=package.raw.SampleNames;
    TarName_list=package.raw.TarNames;
    if length(fieldnames(package))>6
        SampName_list(package.prevsettings.ctrlSample)
        TarName_list(package.prevsettings.ctrlTarget)
        temp_bttn = questdlg('Do you want to use previously used controls?','Prev. Settings','Yes','No','Yes');
        if strcmp(temp_bttn,'Yes')
            ctrl_Tar_ind = package.prevsettings.ctrlTarget;
            ctrl_Samp_ind = package.prevsettings.ctrlSample;
        else
            savelater=true;
        end
    else
        savelater=true;
    end  
    if savelater
        raw=package.raw;
    end
else
    [dat,head,raw]=xlsread(filename,'Results','A41:AK650');
    savelater=true;
    rn=check_range_row(raw);
    
    WellPos=raw(rn,2);
    Omit=raw(rn,3);
    SampName=raw(rn,4);
    TarName=raw(rn,5);
    Ct_raw=ones(max(size(rn)),1);
    Ct_raw=-1*Ct_raw;
    for k=1:max(size(rn))
        if isa(raw{rn(k),15},'double')
            Ct_raw(k)=raw{rn(k),15};
        else
            if NaN_is_num
                Ct_raw(k)=NaNnum;
            else
                Ct_raw(k)=NaN;
            end
        end
    end
    
    
    [SampName_s1,ti]=sort(SampName);
    TarName_s1=TarName(ti);
    WellPos_s1=WellPos(ti);
    Ct_s1=Ct_raw(ti);
    Omit_s1=Omit(ti);
    
    [TarName_s2,ti2]=sort(TarName_s1);
    SampName_s2=SampName_s1(ti2);
    WellPos_s2=WellPos_s1(ti2);
    Ct_s2=Ct_s1(ti2);
    Omit_s2=Omit_s1(ti2);
    
    SampName_list=unique(SampName_s2);
    TarName_list=unique(TarName_s2);
    
    Ct_samp_by_tar=[];r=[];c=[];
    track_depth=ones(max(size(SampName_list)),max(size(TarName_list)));
    for k=1:max(size(Ct_s2))
        r=findchareq(SampName_s2{k},SampName_list);
        c=findchareq(TarName_s2{k},TarName_list);
        Ct_samp_by_tar(r,c,track_depth(r,c))=Ct_s2(k);
        track_depth(r,c)=track_depth(r,c)+1;
    end
    for a=1:size(Ct_samp_by_tar,1)
        for b=1:size(Ct_samp_by_tar,2)
            for c=1:size(Ct_samp_by_tar,3)
                if Ct_samp_by_tar(a,b,c)==0
                    Ct_samp_by_tar(a,b,c)=NaN;
                end
            end
        end
    end
    raw.Ct_samp_by_tar=Ct_samp_by_tar;
    raw.SampleNames=SampName_list;
    raw.TarNames=TarName_list;
end
%}

ctrl_Tar_ind=findchareq(ctrl_Tar,TarName_list);


[size_mat(1),size_mat(2),size_mat(3)]=size(Ct_samp_by_tar);
dCt_samp_by_tar=NaN(size_mat(1)-1,size_mat(2)-1,size_mat(3));
firsttime=true; newTarName_inds=[]; %makes new Target Name list minus the control
newSampName_inds=[];
for a=1:size_mat(1)
    
    newSampName_inds=[newSampName_inds,a];
    
    aa=a;
    
    for b=1:size_mat(2)
        if length(ctrl_Tar_ind)==1
            if b<ctrl_Tar_ind
                dCt_samp_by_tar(aa,b,:)=Ct_samp_by_tar(a,ctrl_Tar_ind,:)-Ct_samp_by_tar(a,b,:);
                if firsttime
                    newTarName_inds=[newTarName_inds,b];
                end
            elseif b>ctrl_Tar_ind
                dCt_samp_by_tar(aa,b-1,:)=Ct_samp_by_tar(a,ctrl_Tar_ind,:)-Ct_samp_by_tar(a,b,:);
                if firsttime
                    newTarName_inds=[newTarName_inds,b];
                end
            end
        else
            dCt_samp_by_tar(aa,b,:)=mean(Ct_samp_by_tar(a,ctrl_Tar_ind,:),2)-Ct_samp_by_tar(a,b,:);
            if firsttime
                newTarName_inds=[newTarName_inds,b];
            end
        end
    end
    if firsttime
        firsttime=false;
    end

end

newTarName_list=TarName_list(newTarName_inds);
newSampName_list=SampName_list(newSampName_inds);
if ~exist('ctrl_Samp_ind','var')
    ctrl_Samp_ind=ctrl_Samp;
    %{
    for n=1:length(ctrl_Samp)
        ctrl_Samp_ind(n)=find(strcmp(newSampName_list,ctrl_Samp{n}));
    end
    %}
    %ctrl_Samp_ind=findchareq(ctrl_Samp,newSampName_list);
end
for n=1:length(grpsamps)
    curgrp=grpsamps{n};
    for m=1:length(curgrp)
        grp_Samp_ind{n}(m)=find(strcmp(newSampName_list,curgrp{m}));
    end
end

dFold_samp_by_tar=2.^dCt_samp_by_tar;


%%
%{
%% norm to E-cad option - only if ran once (and saved) normally (i.e. without this)
normed_ecad=false;
if ~savelater
    ecad_newTar_ind1=findchareq('e-cad',newTarName_list);
    ecad_newTar_ind2=findchareq('cdh1',newTarName_list);
    if ecad_newTar_ind1<=ecad_newTar_ind2
        ecad_newTar_ind=ecad_newTar_ind1;
    else
        ecad_newTar_ind=ecad_newTar_ind2;
    end
    if ecad_newTar_ind<=length(newTarName_list)
        temp_bttn = questdlg('Do you want to normalize to E-cadherin?','Norm to E-cad','Yes','No','No');
    else
        temp_bttn = 'No';
    end
    if strcmp(temp_bttn,'Yes')
        [newsize_mat(1),newsize_mat(2),newsize_mat(3)]=size(dFold_samp_by_tar);
        new_dFold_samp_by_tar=zeros(newsize_mat(1),newsize_mat(2),newsize_mat(3));
        for a=1:newsize_mat(1)
            for b=1:newsize_mat(2)
                if b==ecad_newTar_ind %keeps old E-cad so you can see comparison if you want
                    new_dFold_samp_by_tar(a,b,:)=dFold_samp_by_tar(a,b,:);
                else
                    new_dFold_samp_by_tar(a,b,:)=dFold_samp_by_tar(a,b,:)./dFold_samp_by_tar(a,ecad_newTar_ind,:);
                end
            end
        end
        normed_ecad=true;
        old_dFold_samp_by_tar=dFold_samp_by_tar;
        dFold_samp_by_tar = new_dFold_samp_by_tar;
        savelater=false;
    end
end
%%
%}
%norm to any/multiple secondary targets
normed_ecad=false;
%{
ctrl_Tar_2nd_ind=select_ctrl_GUI('Sample Control',newTarName_list,'Choose secondary normalization:');
if any(ctrl_Tar_2nd_ind>0)
    [newsize_mat(1),newsize_mat(2),newsize_mat(3)]=size(dFold_samp_by_tar);
    new_dFold_samp_by_tar=zeros(newsize_mat(1),newsize_mat(2),newsize_mat(3));
    for a=1:newsize_mat(1)
        for b=1:newsize_mat(2)
            if any(b==ctrl_Tar_2nd_ind) %keeps old secondary controls so you can see comparison if you want
                new_dFold_samp_by_tar(a,b,:)=dFold_samp_by_tar(a,b,:);
            else
                new_dFold_samp_by_tar(a,b,:)=dFold_samp_by_tar(a,b,:)./geomean(dFold_samp_by_tar(a,ctrl_Tar_2nd_ind,:),2);
            end
        end
    end
    normed_ecad=false;
    old_dFold_samp_by_tar=dFold_samp_by_tar;
    dFold_samp_by_tar = new_dFold_samp_by_tar;
    %savelater=false;
end
%}
%%
if graph_param.geo
    tempmean=geomeanNaN3(dFold_samp_by_tar);
else
    tempmean=meanNaN3(dFold_samp_by_tar);
end
    
tempsd=dFold_samp_by_tar;%stdNaN3(dFold_samp_by_tar);
mean_dFold_samp_by_tar=[];
sd_dFold_samp_by_tar=[];
sd_dFold_samp_by_tar_m=[];sd_dFold_samp_by_tar_p=[];
for a=1:length(grp_Samp_ind)
    if graph_param.geo
        mean_dFold_samp_by_tar=[mean_dFold_samp_by_tar;geomean(tempmean(grp_Samp_ind{a},:),1)];
        [tempsdp,tempsdm]=geostd_bounds(dFold_samp_by_tar(grp_Samp_ind{a},:));
        sd_dFold_samp_by_tar_m=[sd_dFold_samp_by_tar_m;tempsdm];
        sd_dFold_samp_by_tar_p=[sd_dFold_samp_by_tar_p;tempsdp];
        
        if a==length(grp_Samp_ind)
            sd_dFold_samp_by_tar(:,:,1)=sd_dFold_samp_by_tar_m;
            sd_dFold_samp_by_tar(:,:,2)=sd_dFold_samp_by_tar_p;
        end
    else
        mean_dFold_samp_by_tar=[mean_dFold_samp_by_tar;nanmean(tempmean(grp_Samp_ind{a},:),1)];
        sd_dFold_samp_by_tar=[sd_dFold_samp_by_tar;nanstd(tempsd(grp_Samp_ind{a},:),0,1)];
    end
end

%{
if length(ctrl_Samp_ind)==1
    mean_ddFold_samp_by_tar=zeros(size_mat(1),size_mat(2)-1);
    sd_ddFold_samp_by_tar=zeros(size_mat(1),size_mat(2)-1);
else
    mean_ddFold_samp_by_tar=zeros(size_mat(1),size_mat(2)-1);
    sd_ddFold_samp_by_tar=zeros(size_mat(1),size_mat(2)-1);
end
%}
for m=1:size(mean_dFold_samp_by_tar,1)
    if graph_param.geo
        mean_ddFold_samp_by_tar(m,:)=mean_dFold_samp_by_tar(m,:)./geomean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
        sd_ddFold_samp_by_tar(m,:,1)=sd_dFold_samp_by_tar(m,:,1)./geomean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
        sd_ddFold_samp_by_tar(m,:,2)=sd_dFold_samp_by_tar(m,:,2)./geomean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
    else
        mean_ddFold_samp_by_tar(m,:)=mean_dFold_samp_by_tar(m,:)./mean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
        sd_ddFold_samp_by_tar(m,:)=sd_dFold_samp_by_tar(m,:)./mean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
    end
end



package.SampleNames=grpnames;
package.TargetNames=rem_tarnames;%newTarName_list;
package.mean=mean_ddFold_samp_by_tar;
package.sd=sd_ddFold_samp_by_tar;
%{
if savelater
    package.raw=raw;
end
%}
%
samp_listlen=size(mean_ddFold_samp_by_tar,1);%samp_listlen=max(size(newSampName_list));
tar_listlen=max(size(newTarName_list));
%{
temp_cp=cell(samp_listlen*2+3,tar_listlen+1);
temp_cp(2:samp_listlen+1,1)=newSampName_list;
temp_cp(4+samp_listlen:(3+2*samp_listlen),1)=newSampName_list;
temp_cp(1,2:tar_listlen+1)=newTarName_list';
temp_cp(3+samp_listlen,2:tar_listlen+1)=newTarName_list';
meancell=cell(samp_listlen,tar_listlen);
sdcell=meancell;
for m=1:samp_listlen
    for n=1:tar_listlen
        meancell{m,n}=mean_ddFold_samp_by_tar(m,n);
        sdcell{m,n}=sd_ddFold_samp_by_tar(m,n);
    end
end
temp_cp(2:samp_listlen+1,2:tar_listlen+1)=meancell;
temp_cp(4+samp_listlen:(3+2*samp_listlen),2:tar_listlen+1)=sdcell;
temp_cp{1,1}='mean';
temp_cp{3+samp_listlen,1}='std-dev';
package.copy2excel=temp_cp;

package.prevsettings.ctrlTarget=ctrl_Tar_ind;
package.prevsettings.ctrlSample=ctrl_Samp_ind;
%}

%{
if savelater
    save(filename2save,'package')
end
%}
%%
if display_Ct
    mean_dCt_samp_by_tar=meanNaN3(dCt_samp_by_tar);
    sd_ddCt_samp_by_tar=stdNaN3(dCt_samp_by_tar);
    for m=1:size(dCt_samp_by_tar,1)
        mean_ddCt_samp_by_tar(m,:)=mean_dCt_samp_by_tar(m,:)-mean_dCt_samp_by_tar(ctrl_Samp_ind,:);
    end

    package.mean=mean_ddCt_samp_by_tar;
    package.sd=sd_ddCt_samp_by_tar;
end
if ~isfield(graph_param,'no_graph')
    figure;
    [hbar,herr]=gui_draw_qpcr_barwitherr(package,graph_param,display_Ct,normed_ecad,linewidth);
    package.handles.hbar=hbar;
    package.handles.herr=herr;
end
end


function rn = check_range_row(c)
k=1;
while true %for k=2:size(c,1)
   if isa(c{k,1},'double')
       if c{k,1}==1
           b=k;
           break
       end
   end
   k=k+1;
end
while true %for k=2:size(c,1)
   if isnan(c{k,1})
       break
   end
   k=k+1;
end
rn=[b:k-1];
end

function ind = findchareq(one,uq_list)
found=false;
for ind=1:max(size(uq_list))
    if strcmpi(one,uq_list{ind})
        found=true;
        break
    end
end
if ~found
    ind=ind+1;
end
end

function out = meanNaN3(in)
%assume dim>2
out=zeros(size(in,1),size(in,2));
for a=1:size(in,1)
    for b=1:size(in,2)
        tempdat=[];
        for c=1:size(in,3)
            if ~isnan(in(a,b,c))
                tempdat=[tempdat,in(a,b,c)];
            end                
        end
        out(a,b)=mean(tempdat);
    end
end
end

function out = geomeanNaN3(in)
%assume dim>2
out=zeros(size(in,1),size(in,2));
for a=1:size(in,1)
    for b=1:size(in,2)
        tempdat=[];
        for c=1:size(in,3)
            if ~isnan(in(a,b,c))
                tempdat=[tempdat,in(a,b,c)];
            end                
        end
        out(a,b)=geomean(tempdat);
    end
end
end

function out = stdNaN3(in)
%assume dim>2
out=zeros(size(in,1),size(in,2));
for a=1:size(in,1)
    for b=1:size(in,2)
        tempdat=[];
        for c=1:size(in,3)
            if ~isnan(in(a,b,c))
                tempdat=[tempdat,in(a,b,c)];
            end                
        end
        out(a,b)=std(tempdat);
    end
end
end

function [outp,outm] = geostd_bounds(infold)
minfold=mean(infold,1);
inlog=log2(infold);
minlog=mean(inlog,1);
sinlog=std(inlog,0,1);

outplog=minlog+sinlog;
outmlog=minlog-sinlog;

outp=(2.^outplog)-minfold;
outm=minfold-(2.^outmlog);
end

function whichone = select_ctrl_GUI(gui_title,list,questionprompt)
if nargin<1 || isempty(gui_title)
    gui_title = '';
end
if nargin<2 || isempty(questionprompt)
    questionprompt = 'What Color Would You Like To Use?';
end

[whichone,ok]=listdlg('ListString',list,'PromptString',questionprompt,'Name',gui_title,'SelectionMode','multiple');

if ~ok
    whichone=-1;
end

end
