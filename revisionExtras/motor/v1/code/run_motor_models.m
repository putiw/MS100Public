function rows = run_motor_models(statsDir, metricsDir, motorDir, outputFile)
% Manuscript-matched ridge models for the 13-tract Motor union.
% Inputs are local immutable aggregate tables; output is a new scratch CSV.
assert(~isfile(outputFile), 'Refusing to overwrite model results.');
Tclin = readtable(fullfile(statsDir,'clinicalScore.xlsx'));
assert(height(Tclin)==132 && numel(unique(string(Tclin.SubjectID)))==132);
Tclin = add_msfc_sdmt(Tclin);
if iscell(Tclin.Gender) || isstring(Tclin.Gender)
    empty = cellfun(@isempty,cellstr(Tclin.Gender));
    Tclin.GenderNum = double(strcmp(Tclin.Gender,'M'));
    Tclin.GenderNum(empty)=NaN;
else
    Tclin.GenderNum = double(Tclin.Gender==1);
end
metrics={'T1','MTR','FA','MD','NDI','ODI','FWF'};
outcomes={'EDSS','MSPro','T25FW','9HPT-D','9HPT-ND','SDMT','MSFC-SDMT'};
source={'EDSS','MSPro','T25FW','x9HPTD','x9HPTND','SDMTcorrect','MSFC_SDMT'};
groups={'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};
motor={'MotorL_Tail','MotorR_Tail','MotorL_NAWM','MotorR_NAWM'};
standard={};
for gi=1:4
    g=groups{gi};
    standard=[standard,{[g 'L_Tail'],[g 'R_Tail'],[g 'L_NAWM'],[g 'R_NAWM']}]; %#ok<AGROW>
end
template=struct('Outcome','','SourceField','','Metric','','MetricDisplay','','Model','', ...
    'Adjusted',0,'Demographics','','CasePolicy','','N',NaN,'PositiveN',NaN, ...
    'NegativeN',NaN,'PerformanceName','','Performance',NaN,'CI_lower',NaN, ...
    'CI_upper',NaN,'Lambda',NaN,'p_raw',NaN,'AIC',NaN, ...
    'SubjectIDSetSHA256','','Status','','Blocker','');
rows=repmat(template,196,1); cursor=0;
Tclass=readtable(fullfile(metricsDir,'ClassicalRegionAllMetrics.csv'));
for mi=1:7
    metric=metrics{mi};
    if mi<=4, standardFile=fullfile(metricsDir,['GroupTract' metric '_All.xlsx']);
    else, standardFile=fullfile(metricsDir,['GroupTract' metric '_All.csv']); end
    Tstd=readtable(standardFile);
    Tmot=readtable(fullfile(motorDir,['Motor' metric '_All.csv']));
    T=innerjoin(innerjoin(Tclin,Tstd,'Keys','SubjectID'),Tmot,'Keys','SubjectID');
    assert(height(T)==131 && numel(unique(string(T.SubjectID)))==131, ...
        'Expected 131 unique clinical/imaging overlap subjects for %s.',metric);
    classicalFields=cellfun(@(x)[x metric], ...
        {'periventricular','juxtacortical','infratentorial','deepwhitematter'}, ...
        'UniformOutput',false);
    for oi=1:7
        binary=oi<=2; yraw=double(T.(source{oi}));
        for adjusted=0:1
            if binary, demoNames={'Age','GenderNum','DurationOfDisease'};
            else, demoNames={'Age','GenderNum'}; end
            if adjusted==0, demoNames={}; end
            Xmot=table2array(T(:,motor)); Xstd=table2array(T(:,standard));
            if isempty(demoNames), Demo=zeros(height(T),0);
            else, Demo=table2array(T(:,demoNames)); end
            valid=all(isfinite(Xmot),2)&all(isfinite(Xstd),2)&isfinite(yraw)&all(isfinite(Demo),2);
            if mi<=4, valid=valid&all(Xmot~=0,2)&all(Xstd~=0,2); end
            if binary
                [tf,loc]=ismember(string(T.SubjectID),string(Tclass.SubjectID));
                good=false(height(T),1);
                found=find(tf);
                C=table2array(Tclass(loc(found),classicalFields));
                good(found)=all(isfinite(C),2);
                if mi<=4, good(found)=good(found)&all(C~=0,2); end
                valid=valid&good&~startsWith(string(T.SubjectID),'sub-C');
            elseif oi<=5
                valid=valid&yraw>0;
            end
            ids=string(T.SubjectID(valid)); n=numel(ids);
            if binary
                if oi==1, y=double(yraw(valid)>3);
                else, y=double(yraw(valid)>1); end
            elseif oi<=5, y=log(yraw(valid));
            else, y=yraw(valid); end
            D=Demo(valid,:); X1=Xmot(valid,:); X2=[Xstd(valid,:),X1];
            models={'Motor-only','Standard+Motor'};
            for modelIndex=1:2
                cursor=cursor+1; r=template;
                r.Outcome=outcomes{oi};r.SourceField=source{oi};r.Metric=metric;
                if strcmp(metric,'FWF'),r.MetricDisplay='ISOVF';else,r.MetricDisplay=metric;end
                r.Model=models{modelIndex};r.Adjusted=adjusted;
                if isempty(demoNames),r.Demographics='none';else,r.Demographics=strjoin(demoNames,'+');end
                if binary,r.CasePolicy='MS joint standard+motor+classical complete cases';
                else,r.CasePolicy='standard+motor complete cases; controls retained';end
                r.N=n;r.SubjectIDSetSHA256=hash_ids(ids);
                if binary
                    r.PositiveN=sum(y);r.NegativeN=n-sum(y);r.PerformanceName='apparent AUC';
                else,r.PerformanceName='apparent R2';end
                r.Status='PENDING'; rows(cursor)=r;
            end
            fprintf('%s %s adjusted=%d N=%d\n',metric,outcomes{oi},adjusted,n);
            if n<20 || (binary && numel(unique(y))~=2)
                for ri=(cursor-1):cursor
                    rows(ri).Status='UNFIT';rows(ri).Blocker='Fewer than 20 valid participants or only one outcome class';
                end
                continue;
            end
            rng(42,'twister');
            xx={X1,X2}; fitted=cell(1,2); errs=cell(1,2);
            for modelIndex=1:2
                try
                    fitted{modelIndex}=fit_one(xx{modelIndex},y,D,binary);
                catch ME
                    errs{modelIndex}=ME.message;
                end
            end
            for modelIndex=1:2
                ri=cursor-2+modelIndex;
                if ~isempty(errs{modelIndex})
                    rows(ri).Status='FIT_ERROR';rows(ri).Blocker=errs{modelIndex};continue;
                end
                try
                    f=fitted{modelIndex};
                    if binary
                        rng(42,'twister');
                        aucs=bootstrp(2000,@(yy,ss) perf_auc(yy,ss),y(:),f.score(:));
                        f.performance=mean(aucs);f.ci=quantile(aucs,[0.025 0.975]);
                    end
                    rows(ri).Performance=f.performance;
                    rows(ri).CI_lower=f.ci(1);rows(ri).CI_upper=f.ci(2);
                    rows(ri).Lambda=f.lambda;rows(ri).p_raw=f.p;rows(ri).AIC=f.AIC;
                    rows(ri).Status='COMPLETED';
                catch ME
                    rows(ri).Status='EVALUATION_ERROR';rows(ri).Blocker=ME.message;
                end
            end
        end
    end
end
assert(cursor==196);
Trows=struct2table(rows);writetable(Trows,outputFile);
fprintf('Wrote %d Motor model rows (%d completed).\n',height(Trows),nnz(strcmp(Trows.Status,'COMPLETED')));
end

function f=fit_one(X,y,Demo,binary)
Xz=zscore(X);
if any(~isfinite(Xz),'all'),error('Nonfinite standardized predictor');end
if binary
    grid=logspace(-6,6,50);
    cv=fitclinear(Xz,y,'Learner','logistic','Regularization','ridge','Lambda',grid,'KFold',10);
    [~,idx]=min(kfoldLoss(cv));lambda=grid(idx);
    ridge=fitclinear(Xz,y,'Learner','logistic','Regularization','ridge','Lambda',lambda);
    LP=Xz*ridge.Beta+ridge.Bias;
    warning('off','all');
    if isempty(Demo),mdl=fitglm(LP,y,'Distribution','binomial');score=LP;
    else,mdl=fitglm([LP,Demo],y,'Distribution','binomial');score=mdl.Fitted.LinearPredictor;end
    warning('on','all');
    coef=mdl.Coefficients.Estimate(2);se=mdl.Coefficients.SE(2);
    p=2*normcdf(-abs(coef/se));
    f=struct('lambda',lambda,'p',p,'AIC',mdl.ModelCriterion.AIC, ...
        'performance',NaN,'ci',[NaN,NaN],'score',score);
else
    grid=logspace(-6,6,100);
    cv=fitrlinear(Xz,y,'Learner','leastsquares','Regularization','ridge','Lambda',grid,'KFold',10);
    [~,idx]=min(kfoldLoss(cv));lambda=grid(idx);
    ridge=fitrlinear(Xz,y,'Learner','leastsquares','Regularization','ridge','Lambda',lambda);
    LP=Xz*ridge.Beta+ridge.Bias;
    if isempty(Demo),mdl=fitlm(LP,y);R2=mdl.Rsquared.Ordinary;
    else,mdl=fitlm([LP,Demo],y);R2=fitlm(mdl.Fitted,y).Rsquared.Ordinary;end
    coef=mdl.Coefficients.Estimate(2);se=mdl.Coefficients.SE(2);
    t=coef/se;p=2*tcdf(-abs(t),mdl.DFE);
    if p==0,p=2*normcdf(-abs(t));end
    f=struct('lambda',lambda,'p',p,'AIC',mdl.ModelCriterion.AIC, ...
        'performance',R2,'ci',[NaN,NaN],'score',[]);
end
end

function value=perf_auc(y,score)
[~,~,~,value]=perfcurve(y,score,1);
end

function value=hash_ids(ids)
ids=sort(string(ids(:)));bytes=unicode2native(char(strjoin(ids,newline)),'UTF-8');
engine=java.security.MessageDigest.getInstance('SHA-256');engine.update(bytes);
raw=typecast(engine.digest(),'uint8');value=lower(reshape(dec2hex(raw,2).',1,[]));
end

function T=add_msfc_sdmt(T)
v=[T.T25FW,T.x9HPTD,T.x9HPTND,T.SDMTcorrect];
valid=all(isfinite(v),2)&all(v(:,1:3)>0,2)&v(:,4)>=0;
ref=valid&~startsWith(string(T.SubjectID),'sub-C');
assert(sum(valid)==127&&sum(ref)==84,'MSFC reference cohort mismatch');
arm=nan(height(T),1);arm(valid)=(1./T.x9HPTD(valid)+1./T.x9HPTND(valid))/2;
zarm=(arm-mean(arm(ref)))/std(arm(ref));
zleg=-(T.T25FW-mean(T.T25FW(ref)))/std(T.T25FW(ref));
zcog=(T.SDMTcorrect-mean(T.SDMTcorrect(ref)))/std(T.SDMTcorrect(ref));
T.MSFC_SDMT=mean([zarm,zleg,zcog],2,'omitmissing');T.MSFC_SDMT(~valid)=NaN;
end
