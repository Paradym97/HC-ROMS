function[beta,r2,SST,ATOM,w]=HCGCM_CFSR_coupling_mrgwr(varname,d_lon,d_lat,d_time,sstname)
%d_lon: e-folding scale (in deg) of the weight in zonal direction 
%d_lat: e-folding scale (in deg) of the weight in meridional direction 
%d_time: e-folding scale (in deg) of the weight in time domain 


if nargin<5
  sstname='SST';
end

ts=5;% subsmapling factor in time domain

% proprecess the data
% dataPath='/home/jingzhao/ym/project/HCGCM/CFSR/NP_daily/';
% filename_SST=[dataPath 'CFSR_SST_NP_2006_daily.nc'];
% filename_ATOM=[dataPath 'CFSR_' varname '_NP_2006_daily.nc'];

dataPath='/home/jingzhao/ym/project/HCGCM/CFSR/';
filename_SST=[dataPath 'CFSR_SST_NP_2006_daily_bigger.nc'];
filename_ATOM=[dataPath 'CFSR_' varname '_NP_2006_daily_bigger.nc'];

t0=tic;
js=1;
tic_toc(js)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SST=ncread(filename_SST,sstname);

q=size(SST,3);
time=1:q;

fprintf(['size of SST= ' num2str(size(SST)) '\n'])
js=js+1;
tic_toc(js)=toc(t0);
[m,n]=size(SST(:,:,1));
SST(SST<-9000)=nan;
SST=double(SST);

parfor i=1:m
    for j=1:n
        SST(i,j,:)=smooth2a(squeeze(SST(i,j,:)),ts,0);% filter out high-frequency variability
%        SST(i,j,:)=squeeze(SST(i,j,:))-nanmean(SST(i,j,:));% remove the time-mean value
    end
end
SST=SST(:,:,1:ts:end);% subsampling in time

%% read ATOM 
ATOM=ncread(filename_ATOM,varname);

ATOM(ATOM<-9000)=nan;
ATOM=double(ATOM);

ATOM(SST<0)=nan;
SST(SST<0)=nan;

parfor i=1:m
    for j=1:n
        ATOM(i,j,:)=smooth2a(squeeze(ATOM(i,j,:)),ts,0);% filter out high-frequency variability
%        ATOM(i,j,:)=squeeze(ATOM(i,j,:))-nanmean(ATOM(i,j,:));% remove the time-mean value
    end
end
ATOM=ATOM(:,:,1:ts:end);% subsampling in time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
js=js+1;
tic_toc(js)=toc(t0);
fprintf(['read data using time ' num2str(tic_toc(js)-tic_toc(js-1)) ' in ' num2str(tic_toc(js)) '\n'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try whether we need this mask
masked = ncread('/home/jingzhao/ym/project/HCGCM/CFSR/cfsr_mask_bigger.nc','mask');
masked = repmat(masked,[1,1,size(SST,3)]);
fprintf(['size of masked = '])
size(masked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recompute SST size
%time=1-d_time*2:q+d_time*2;
fprintf(['size of time= ' num2str(size(time)) '\n'])
time=time(1:ts:end)
q=length(time);
[m,n]=size(SST(:,:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:q
    SST(:,:,t)=SST(:,:,t)-smooth2a(SST(:,:,t),4,4);%isolate mesoscale anomaly
    ATOM(:,:,t)=ATOM(:,:,t)-smooth2a(ATOM(:,:,t),4,4);%isolate mesoscale anomaly
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SST(masked==2)=nan;
ATOM(masked==2)=nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
js=js+1;
tic_toc(js)=toc(t0);
fprintf(['smooth2a using time ' num2str(tic_toc(js)-tic_toc(js-1)) ' in ' num2str(tic_toc(js)) '\n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the weight
index_lon=round(2*d_lon/0.3125);
index_lat=round(2*d_lat/0.3125);% should be modifed if the grid size is not equal to 0.25 deg
index_time=round(2*d_time/(time(2)-time(1)));

w=nan(2*index_lon+1,2*index_lat+1,2*index_time+1);
wlon=linspace(-2,2,2*index_lon+1);
wlat=linspace(-2,2,2*index_lat+1);
wtime=linspace(-2,2,2*index_time+1);

% exponential weighting function
for t=1:2*index_time+1
    w(:,:,t)=exp(-abs(wtime(t))^2)*(exp(-abs(wlon).^2)'*exp(-abs(wlat).^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
js=js+1;
tic_toc(js)=toc(t0);
fprintf(['w using time ' num2str(tic_toc(js)-tic_toc(js-1)) ' in ' num2str(tic_toc(js)) '\n'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defin output
%beta=nan(m-2*index_lon,n-2*index_lat,q-2*index_time,3);
%r2=nan(m-2*index_lon,n-2*index_lat,q-2*index_time);
beta=nan(m,n,q,3);
r2=nan(m,n,q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute beta and r2
fprintf(['time[0]= ' num2str(time(1+index_time)) 'time[end]= ' num2str(time(q-index_time)) '\n'])
for t=1:q
    t
tic
    parfor i=1:m
    % for i=1:m
        for j=1:n
            aaa=max(1,i-index_lon);
            bbb=min(m,i+index_lon);
            ccc=max(1,j-index_lat);
            ddd=min(n,j+index_lat);         
            eee=max(1,t-index_time);
            fff=min(q,t+index_time);
            i2=i-index_lon;
            j2=j-index_lat;
            t2=t-index_time;

            aa = aaa - (i-index_lon)+1;
            bb = bbb - (i+index_lon)+2*index_lon+1;
            cc = ccc - (j-index_lat)+1;
            dd = ddd - (j+index_lat)+2*index_lat+1;
            ee = eee - (t-index_time)+1;
            ff = fff - (t+index_time)+2*index_time+1;

            xxx = nan(2*index_lon+1,2*index_lat+1,2*index_time+1);
            yyy = nan(2*index_lon+1,2*index_lat+1,2*index_time+1);

            xxx(aa:bb, cc:dd, ee:ff) = SST(aaa:bbb,ccc:ddd,eee:fff);
            yyy(aa:bb, cc:dd, ee:ff) = ATOM(aaa:bbb,ccc:ddd,eee:fff);

            % fprintf([num2str(i) ' ' num2str(j) ' ' num2str(t) ' || ' ...
            % num2str(aaa) ' ' num2str(bbb) ' '  num2str(ccc) ' '  num2str(ddd) ' '  num2str(eee) ' '  num2str(fff) ' || '...
            % num2str(aa)  ' ' num2str(bb) ' '  num2str(cc) ' '  num2str(dd) ' '  num2str(ee) ' '  num2str(ff) '\n'])

            % size(xxx)
            % size(yyy)
            % size(w)

%            if nanvar(squeeze(SST(i,j,:)))>0.01  by yman
            if ~all(isnan(squeeze(SST(i,j,t)))) && ~all(isnan(squeeze(ATOM(i,j,t)))) 
                [beta(i,j,t,:),r2(i,j,t)]=HCGCM_st_weighted_regression(xxx, yyy,w);
            end
        end
    end
toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
js=js+1;
tic_toc(js)=toc(t0);
fprintf(['regression using time ' num2str(tic_toc(js)-tic_toc(js-1)) ' in ' num2str(tic_toc(js)) '\n']) 
end