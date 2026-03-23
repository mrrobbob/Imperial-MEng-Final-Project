function [Spiked_data,Spike_index]=AutoSpike_Matrix(data,WvnAxis,WvnCropIdx)

Spiked_data=zeros(size(data(:,WvnCropIdx:end)));
Spike_idx=zeros(length(data(:,1)),1);
h=waitbar(0,'Processing spectra...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

for ii=1:length(data(:,1));
     if getappdata(h,'canceling')
        break
    end
    waitbar(ii/length(data(:,1)),h,sprintf('%2.2f %%',100*ii/length(data(:,1))));
    
    [Spiked_data(ii,:),Spike_idx(ii,1)]=AutoSpike_v2(data(ii,:),WvnAxis,WvnCropIdx);
    
    Spike_index=find(Spike_idx>0);
end
delete(h);

function [test_flat,spike_flag]=AutoSpike_v2(test,WvnAxis,WvnCropIdx)
% clear test d1 s ii jj w2 threshold count zz reduced order index reorder repeat base1 base2 test_flat
% %%
% test=Data1(2378,:)';
% test=D1{14295,1};
d0=zeros(length(test),1); % initialize peak vector
d2=d0;d3=d0;
s=length(test);
w2=1; % define window for peak detection (2*w2+1 == 7 pixels)
% % % % threshold=2.5e3;
% % % % threshold=.1*max(test);
% % % % threshold=.1*(max(test)-min(test));
% % threshold=max([250;.1*(max(test)-min(test))]);
threshold=max([250;.1*(max(test(WvnCropIdx:end))-min(test(WvnCropIdx:end)))]);

% if sum(test)/length(WvnAxis)>75
    % % % % %     Threshold>100
    for ii =w2+1:(s-w2-1) % moving window peak locator: if pixel intensity = maximum intensity within +/- 3 pixels, record as a peak
        if max(test(ii-w2:ii+w2))==test(ii) && test(ii)>=mean([test(ii-1),test(ii+1)])+threshold %peak at ii
            d0(ii)=ii;
        end
    end
    d0=d0(d0~=0); % remove zeros from peak list
    %%%% Check for nearby cosmic echo at lower intensity
    for ii =w2+1:(s-w2-1) % moving window peak locator: if pixel intensity = maximum intensity within +/- 3 pixels, record as a peak
        if max(test(ii-w2:ii+w2))==test(ii) && test(ii)>0.75*threshold+min(test)%peak at ii
            d2(ii)=ii;
        end
    end
    d2=d2(d2~=0); % remove zeros from peak list
    % d1=d0;
    % for ii=1:length(d0)
    %     if find(abs(d0(ii)-d2)<10 & abs(d0(ii)-d2)>0)
    %         d1=[d1; d2(find(abs(d0(ii)-d2)<10 & abs(d0(ii)-d2)>0))];
    %     end
    % end
    
    %%% Check for evenly distributed peak aboe threshold
%     for ii =2+1:(s-2-1) % moving window peak locator: if pixel intensity = maximum intensity within +/- 3 pixels, record as a peak
%         if max(test(ii-2:ii+2))==test(ii) && test(ii)>=mean([test(ii-2),test(ii+2)])+0.75*threshold %peak at ii
%             d3(ii)=ii;
%         end
%     end
%     d3=d3(d3~=0); % remove zeros from peak list
%     d1=[d0;d3(:)];
d1=[d0];
if isempty(d1) && threshold>500
    d1=d2;
end
    for ii=1:length(d1)
        if find(abs(d1(ii)-d2)<7 & abs(d1(ii)-d2)>0)
            d1=[d1; d2(find(abs(d1(ii)-d2)<7 & abs(d1(ii)-d2)>0))];
        end
    end
    d1=unique(d1);
    %% Appropriate Wavenumber check
    d1=d1(WvnAxis(d1)>WvnAxis(WvnCropIdx));
    %%
    clear d0 d2 d3
    if ~isempty(d1)
        spike_flag=1;
    else
        spike_flag=0;
    end
    % %%
    % test_orig=test;
    test_flat=test;
    if length(d1)<10
        for jj=1:length(d1)
            if (~isempty(diff(d1)) && max(diff(d1))>5) && d1(jj)<=s-15 || isempty(diff(d1)) && d1(jj)>6 && d1(jj)<s-5
                %         base1=find((abs(diff(test(d1(jj)-4:d1(jj)+4)))>threshold)==1,1,'first');
                %         base2=find((abs(diff(test(d1(jj)-4:d1(jj)+4)))>threshold)==1,1,'last');
                base1=find((abs(diff(test(d1(jj)-4:d1(jj)+4)))>0.45*threshold)==1,1,'first');
                base2=find((abs(diff(test(d1(jj)-4:d1(jj)+4)))>0.45*threshold)==1,1,'last');
                if ~isempty(base1) && ~isempty(base2)
                    test_flat(d1(jj)-5+base1:d1(jj)-4+base2)=mean([test(d1(jj)-6+base1);test(d1(jj)-3+base2)])+std(detrend(test(d1(jj)-3+base2:d1(jj)+5)),[],2).*randn(length(d1(jj)-5+base1:d1(jj)-4+base2),1);
                end
            elseif (~isempty(diff(d1)) && max(diff(d1))>5) && d1(jj)<=s-15 || isempty(diff(d1)) && d1(jj)<=6
                %         base1=find((abs(diff(test(d1(jj)-5:d1(jj)+5)))>threshold)==1,1,'first');
                base2=find((abs(diff(test(1:d1(jj)+5)))>threshold)==1,1,'last');
                test_flat(1:d1(jj)+base2)=mean([test(d1(jj)+6);test(d1(jj)+6+base2)])+std(detrend(test(d1(jj)-4+base2:d1(jj)+4)),[],2).*randn(length(1:d1(jj)+base2),1);
            elseif (~isempty(diff(d1)) && max(diff(d1))>5) && d1(jj)<=s-15 || isempty(diff(d1)) && d1(jj)>=s-5
                base1=find((abs(diff(test(d1(jj)-5:s)))>threshold)==1,1,'first');
                %         base2=find((abs(diff(test(1:d1(jj)+5)))>threshold)==1,1,'last');
                test_flat(d1(jj)-5+base1:s)=mean([test(d1(jj)-8+base1);test(d1(jj)-6)])+std(detrend(test(d1(jj)-4:d1(jj)-4+base1)),[],2).*randn(length(d1(jj)-5+base1:s),1);
                %     elseif min(d1)<15
                %         %         base1=find((abs(diff(test(d1(jj)-15:d1(jj)+15)))>threshold)==1,1,'first');
                %         base2=find((abs(diff(test(1:d1(end)+15)))>threshold)==1,1,'last');
                %         test_flat(1:d1(jj)-12+base2)=mean([test(1);test(d1(jj)-5+base2)])+100*randn(length(1:d1(jj)-12+base2),1);
            elseif d1(jj)>=s-15 %(~isempty(diff(d1)) && max(diff(d1))<5) ||
                %         base1=find((abs(diff(test(d1(jj)-5:d1(jj)+5)))>threshold)==1,1,'first');
                %         %         base2=find((abs(diff(test(1:d1(end)+15)))>threshold)==1,1,'last');
                %         test_flat(d1(jj)-5+base1:s)=mean([test(d1(jj)-8+base1);test(d1(jj)-6)])+100*randn(length(d1(jj)-5+base1:s),1);
                base1=find((abs(diff(test(d1(jj)-5:d1(jj)+5)))>0.45*threshold)==1,1,'first');
                %         base2=find((abs(diff(test(1:d1(end)+15)))>threshold)==1,1,'last');
                test_flat(d1(jj)-5+base1:s)=mean([test(d1(jj)-8+base1);test(d1(jj)-6)])+std(detrend(test(d1(jj)-12:d1(jj)-4)),[],2)*randn(length(d1(jj)-5+base1:s),1);
                break;
            else
                base1=find((abs(diff(test(d1(jj)-15:d1(jj)+15)))>0.45*threshold)==1,1,'first');
                %         base2=find((abs(diff(test(d1(end)-15:d1(end)+14)))>threshold)==1,1,'last')
                base2=find((abs(diff(test(d1(jj)-15:d1(jj)+14)))>0.45*threshold)==1,1,'last');
                if ~isempty(base1) && ~isempty(base2)
                    test_flat(d1(jj)-18+base1:d1(jj)-14+base2)=abs(mean([test(d1(jj)-20+base1);test(d1(jj)-5+base2)])+std([detrend(test(d1(1)-23+base1:d1(1)-18+base1)),detrend(test(d1(end)-14+base2:d1(end)-9+base2))]).*randn(length(d1(jj)-18+base1:d1(jj)-14+base2),1));
                end
            end
            test=test_flat;
            clear base1 base2
        end
    end
    test_flat=test_flat(WvnCropIdx:end);
    % % % % figure;plot(X5{1,1},test,'b-',X5{1,1}(d1),test(d1),'ro',X5{1,1},test_flat,'k');
% else
%     clear d0 d2 d3
%     test_flat=test;
%     spike_flag=0;
% end