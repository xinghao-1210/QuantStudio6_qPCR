function [outraw] = gui_read_QuantStudio6_qpcrdata(filename,filename2save)


%%
%$Hard-code
display_Ct=false;
accepted_controls={'CphA','GAPDH','18S','GUSB','HPRT1'};
NaN_is_num=true;
NaNnum=40;
linewidth=3.5;
%display_Ct=true;

%%

if exist([filename2save '.mat'])==2
    load([filename2save '.mat']);
    savelater=false;
    raw=package.raw;
    Ct_samp_by_tar=package.raw.Ct_samp_by_tar;
    SampName_list=package.raw.SampleNames;
    TarName_list=package.raw.TarNames;
    %{
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
    %}
else
    [dat,head,raw]=xlsread(filename,'Results');
    savelater=true;
    rn=check_range_row(raw);
    
    WellPos=raw(rn,2);
    Omit=raw(rn,3);
    SampName=raw(rn,4);
    TarName=raw(rn,5);
    Ct_raw=ones(max(size(rn)),1);
    Ct_raw=-1*Ct_raw;
    ct_colidx=find(strcmpi('CT',raw(rn(1)-1,:)));
    if isempty(ct_colidx)
        temp=find(strncmpi('CT',raw(rn(1)-1,:),1));
        ct_colidx=temp(1);
    end
    for k=1:max(size(rn))
        if isa(raw{rn(k),ct_colidx},'double')
            Ct_raw(k)=raw{rn(k),ct_colidx};
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
    outraw.Ct_samp_by_tar=Ct_samp_by_tar;
    outraw.SampleNames=SampName_list;
    outraw.TarNames=TarName_list;

end
%{
if UI_later
    for n=1:length(accepted_controls)
        tempctrl_Tar_ind(n)=findchareq(accepted_controls{n},TarName_list);        
    end
    validtempctrl_Tar_ind=tempctrl_Tar_ind<=length(TarName_list);
    if ~exist('ctrl_Tar_ind','var')
        if sum(validtempctrl_Tar_ind)~=1
            ctrl_Tar_ind=select_ctrl_GUI('Target Control',TarName_list,'What target do you want to use as a control?');
        else
            tempind=find(validtempctrl_Tar_ind==1);
            ctrl_Tar_ind=tempctrl_Tar_ind(tempind);
        end
    end
else
    ctrl_Tar_ind=findchareq(ctrl_Tar,TarName_list);
end
noRT_Samp_ind=findchareq('no RT',SampName_list);
if noRT_Samp_ind>length(SampName_list)
   noRT_Samp_ind=findchareq('noRT',SampName_list); 
end

[size_mat(1),size_mat(2),size_mat(3)]=size(Ct_samp_by_tar);
dCt_samp_by_tar=NaN(size_mat(1)-1,size_mat(2)-1,size_mat(3));
firsttime=true; newTarName_inds=[]; %makes new Target Name list minus the control
newSampName_inds=[];
for a=1:size_mat(1)
    if a==noRT_Samp_ind
    else
        newSampName_inds=[newSampName_inds,a];
        if a>noRT_Samp_ind
            aa=a-1;
        else
            aa=a;
        end
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
end

newTarName_list=TarName_list(newTarName_inds);
newSampName_list=SampName_list(newSampName_inds);
if ~exist('ctrl_Samp_ind','var')
    if UI_later
        ctrl_Samp_ind=select_ctrl_GUI('Sample Control',newSampName_list,'What sample do you want to use as a control?');
    else
        ctrl_Samp_ind=findchareq(ctrl_Samp,newSampName_list);
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
if ~savelater
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
        savelater=false;
    end
end
%%


mean_dFold_samp_by_tar=meanNaN3(dFold_samp_by_tar);
sd_dFold_samp_by_tar=stdNaN3(dFold_samp_by_tar);
if length(ctrl_Samp_ind)==1
    mean_ddFold_samp_by_tar=zeros(size_mat(1)-1,size_mat(2)-1);
    sd_ddFold_samp_by_tar=zeros(size_mat(1)-1,size_mat(2)-1);
else
    mean_ddFold_samp_by_tar=zeros(size_mat(1)-1,size_mat(2));
    sd_ddFold_samp_by_tar=zeros(size_mat(1)-1,size_mat(2));
end
for m=1:size(dFold_samp_by_tar,1)
    mean_ddFold_samp_by_tar(m,:)=mean_dFold_samp_by_tar(m,:)./mean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
    sd_ddFold_samp_by_tar(m,:)=sd_dFold_samp_by_tar(m,:)./mean(mean_dFold_samp_by_tar(ctrl_Samp_ind,:),1);
end

package.SampleNames=newSampName_list;
package.TargetNames=newTarName_list;
package.mean=mean_ddFold_samp_by_tar;
package.sd=sd_ddFold_samp_by_tar;
if savelater
    package.raw=raw;
end
%
samp_listlen=max(size(newSampName_list));
tar_listlen=max(size(newTarName_list));
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

[hbar,herr]=draw_qpcr_barwitherr(package,display_Ct,normed_ecad,linewidth);
package.handles.hbar=hbar;
package.handles.herr=herr;
%}
end


function rn = check_range_row(c)
k=1;
while true %for k=2:size(c,1)
   if isa(c{k,3},'char')
       k=k+1;
       b=k;
       break
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
