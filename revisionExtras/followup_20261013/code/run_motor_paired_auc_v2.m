function Trows = run_motor_paired_auc_v2(statsDir,metricsDir,motorDir,motorReference,edssReference,msproReference,outputFile)
%RUN_MOTOR_PAIRED_AUC_V2 Compare Motor scores with the canonical baseline fits.
% Motor-only and Standard+Motor are fitted in their v1 order with rng(42).
% Standard and Classical are a separate baseline pair, refitted with their
% original Standard-then-Classical ridge order after a fresh rng(42).
% Scores and subject IDs remain in memory; the CSV contains aggregates only.
assert(~isfile(outputFile),'Refusing to overwrite paired result CSV.');
for path={motorReference,edssReference,msproReference}
    assert(isfile(path{1}),'Missing reference: %s',path{1});
end
TmotorRef=readtable(motorReference,'TextType','string');
TedssRef=readtable(edssReference,'TextType','string');
TmsproRef=readtable(msproReference,'TextType','string');
assert(height(TedssRef)==28 && all(TedssRef.Pass==1));
TmsproRef=TmsproRef(strcmp(string(TmsproRef.Outcome),'MSPro'),:);
assert(height(TmsproRef)==28);
Tclin=readtable(fullfile(statsDir,'clinicalScore.xlsx'));
assert(height(Tclin)==132 && numel(unique(string(Tclin.SubjectID)))==132);
if iscell(Tclin.Gender) || isstring(Tclin.Gender)
    blank=cellfun(@isempty,cellstr(Tclin.Gender));
    Tclin.GenderNum=double(strcmp(Tclin.Gender,'M'));
    Tclin.GenderNum(blank)=NaN;
else
    Tclin.GenderNum=double(Tclin.Gender==1);
end
Tclass=readtable(fullfile(metricsDir,'ClassicalRegionAllMetrics.csv'));
assert(height(Tclass)==87 && numel(unique(string(Tclass.SubjectID)))==87);
metrics={'T1','MTR','FA','MD','NDI','ODI','FWF'};
groups={'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};
regions={'periventricular','juxtacortical','infratentorial','deepwhitematter'};
motorFields={'MotorL_Tail','MotorR_Tail','MotorL_NAWM','MotorR_NAWM'};
standardFields={};
for gi=1:numel(groups)
    g=groups{gi};
    standardFields=[standardFields,{[g 'L_Tail'],[g 'R_Tail'], ...
        [g 'L_NAWM'],[g 'R_NAWM']}]; %#ok<AGROW>
end
template=struct('Outcome','','Metric','','MetricDisplay','','Adjusted',0, ...
    'Comparison','','ModelA','','ModelB','','N',NaN,'PositiveN',NaN, ...
    'NegativeN',NaN,'AUC_A_Point',NaN,'AUC_B_Point',NaN, ...
    'AUC_A_BootstrapMean',NaN,'AUC_B_BootstrapMean',NaN, ...
    'DeltaAUC_Point',NaN,'DeltaAUC_BootstrapMean',NaN, ...
    'DeltaAUC_CI_Lower',NaN,'DeltaAUC_CI_Upper',NaN, ...
    'AIC_A',NaN,'AIC_B',NaN,'Lambda_A',NaN,'Lambda_B',NaN, ...
    'SubjectIDSetSHA256','','BootstrapIndexSHA256','', ...
    'BaselineReferenceMatch',false,'Status','');
rows=repmat(template,84,1); cursor=0;
for mi=1:numel(metrics)
    metric=metrics{mi};
    if mi<=4, suffix='.xlsx'; else, suffix='.csv'; end
    Tstd=readtable(fullfile(metricsDir,['GroupTract' metric '_All' suffix]));
    Tmot=readtable(fullfile(motorDir,['Motor' metric '_All.csv']));
    T=innerjoin(innerjoin(Tclin,Tstd,'Keys','SubjectID'),Tmot,'Keys','SubjectID');
    assert(height(T)==131 && numel(unique(string(T.SubjectID)))==131);
    Xstd=table2array(T(:,standardFields));
    Xmot=table2array(T(:,motorFields));
    classicalFields=cellfun(@(r)[r metric],regions,'UniformOutput',false);
    [present,loc]=ismember(string(T.SubjectID),string(Tclass.SubjectID));
    Xclass=nan(height(T),4);
    Xclass(present,:)=table2array(Tclass(loc(present),classicalFields));
    for oi=1:2
        if oi==1
            outcome='EDSS'; raw=double(T.EDSS); threshold=3;
        else
            outcome='MSPro'; raw=double(T.MSPro); threshold=1;
        end
        for adjusted=0:1
            if adjusted
                Demo=table2array(T(:,{'Age','GenderNum','DurationOfDisease'}));
            else
                Demo=zeros(height(T),0);
            end
            valid=all(isfinite(Xmot),2) & all(isfinite(Xstd),2) & ...
                all(isfinite(Xclass),2) & isfinite(raw) & ...
                all(isfinite(Demo),2) & ~startsWith(string(T.SubjectID),'sub-C');
            if mi<=4
                valid=valid & all(Xmot~=0,2) & all(Xstd~=0,2) & all(Xclass~=0,2);
            end
            ids=string(T.SubjectID(valid)); n=numel(ids);
            y=double(raw(valid)>threshold); D=Demo(valid,:);
            assert(n>=20 && numel(unique(y))==2);
            idHash=hash_ids(ids);
            % Preserve the completed Motor package's fitting sequence.
            rng(42,'twister');
            motorOnly=fit_binary(Xmot(valid,:),y,D);
            standardMotor=fit_binary([Xstd(valid,:),Xmot(valid,:)],y,D);
            % This reset is essential: the original baseline pair started
            % independently with Standard, then fitted Classical.
            rng(42,'twister');
            [standardLP,standardLambda]=fit_ridge_lp(Xstd(valid,:),y);
            [classicalLP,classicalLambda]=fit_ridge_lp(Xclass(valid,:),y);
            standard=eval_binary(standardLP,y,D,standardLambda);
            classical=eval_binary(classicalLP,y,D,classicalLambda);
            fits={motorOnly,standardMotor,standard,classical};
            names={'Motor-only','Standard+Motor','Standard','Classical'};
            for modelIndex=1:2
                ref=TmotorRef(strcmp(TmotorRef.Metric,metric) & ...
                    strcmp(TmotorRef.Outcome,outcome) & ...
                    TmotorRef.Adjusted==adjusted & ...
                    strcmp(TmotorRef.Model,names{modelIndex}),:);
                assert(height(ref)==1 && strcmp(ref.Status,'COMPLETED'));
                assert(ref.N==n && ref.PositiveN==sum(y) && ...
                    strcmp(char(ref.SubjectIDSetSHA256),idHash));
                assert(abs(ref.AIC-fits{modelIndex}.AIC)<1e-6);
                assert(abs(ref.Lambda-fits{modelIndex}.lambda)<1e-10);
            end
            if oi==1
                baselineRef=TedssRef(strcmp(TedssRef.Metric,metric) & ...
                    TedssRef.Adjusted==adjusted,:);
            else
                baselineRef=TmsproRef(strcmp(TmsproRef.Metric,metric) & ...
                    TmsproRef.Adjusted==adjusted,:);
            end
            for modelIndex=3:4
                if modelIndex==3, refName='Tract-based'; else, refName='Classical-region'; end
                ref=baselineRef(strcmp(baselineRef.Model,refName),:);
                assert(height(ref)==1 && ref.N==n && ref.PositiveN==sum(y));
                if oi==1
                    assert(strcmp(char(ref.SubjectIDSetSHA256),idHash));
                    referenceAIC=ref.CanonicalAIC;
                    referenceLambda=ref.CanonicalLambda;
                    assert(abs(ref.ApparentAUCPoint-auc_rank(y,fits{modelIndex}.score))<1e-12);
                else
                    referenceAIC=ref.AIC; referenceLambda=ref.Lambda;
                end
                assert(abs(referenceAIC-fits{modelIndex}.AIC)<1e-6, ...
                    '%s %s: baseline AIC differs from canonical reference',metric,outcome);
                lambdaTolerance=1e-12*max(1,abs(referenceLambda));
                assert(abs(referenceLambda-fits{modelIndex}.lambda)<lambdaTolerance, ...
                    '%s %s: baseline lambda differs from canonical reference',metric,outcome);
            end
            pairs=[1 3;2 3;1 4]; bootstrapHash='';
            for pi=1:size(pairs,1)
                a=pairs(pi,1); b=pairs(pi,2);
                paired=paired_stratified_auc(y,fits{a}.score,fits{b}.score,2000,.05,42);
                assert(paired.all_two_class && paired.positive_per_resample==sum(y));
                if isempty(bootstrapHash)
                    bootstrapHash=paired.index_sha256;
                else
                    assert(strcmp(bootstrapHash,paired.index_sha256));
                end
                cursor=cursor+1; r=template;
                r.Outcome=outcome; r.Metric=metric;
                if strcmp(metric,'FWF'), r.MetricDisplay='ISOVF';
                else, r.MetricDisplay=metric; end
                r.Adjusted=adjusted;
                r.Comparison=[names{a} ' minus ' names{b}];
                r.ModelA=names{a}; r.ModelB=names{b};
                r.N=n; r.PositiveN=sum(y); r.NegativeN=n-sum(y);
                r.AUC_A_Point=paired.tract.point_auc;
                r.AUC_B_Point=paired.classical.point_auc;
                r.AUC_A_BootstrapMean=paired.tract.mean_bootstrap_auc;
                r.AUC_B_BootstrapMean=paired.classical.mean_bootstrap_auc;
                r.DeltaAUC_Point=paired.delta_point;
                r.DeltaAUC_BootstrapMean=paired.delta_mean;
                r.DeltaAUC_CI_Lower=paired.delta_ci(1);
                r.DeltaAUC_CI_Upper=paired.delta_ci(2);
                r.AIC_A=fits{a}.AIC; r.AIC_B=fits{b}.AIC;
                r.Lambda_A=fits{a}.lambda; r.Lambda_B=fits{b}.lambda;
                r.SubjectIDSetSHA256=idHash;
                r.BootstrapIndexSHA256=bootstrapHash;
                r.BaselineReferenceMatch=true;
                r.Status='COMPLETED'; rows(cursor)=r;
            end
            fprintf('%s %s adjusted=%d N=%d: canonical baseline matched, three paired deltas complete\n', ...
                metric,outcome,adjusted,n);
        end
    end
end
assert(cursor==84);
Trows=struct2table(rows);
assert(~isfile(outputFile),'Output appeared during run; refusing overwrite.');
writetable(Trows,outputFile);
end

function fit=fit_binary(X,y,D)
[LP,lambda]=fit_ridge_lp(X,y);
fit=eval_binary(LP,y,D,lambda);
end

function [LP,lambda]=fit_ridge_lp(X,y)
Xz=zscore(X);
assert(all(isfinite(Xz),'all'));
grid=logspace(-6,6,50);
cv=fitclinear(Xz,y,'Learner','logistic','Regularization','ridge', ...
    'Lambda',grid,'KFold',10);
[~,index]=min(kfoldLoss(cv)); lambda=grid(index);
ridge=fitclinear(Xz,y,'Learner','logistic','Regularization','ridge','Lambda',lambda);
LP=Xz*ridge.Beta+ridge.Bias;
end

function fit=eval_binary(LP,y,D,lambda)
warningState=warning; cleanup=onCleanup(@()warning(warningState));
warning('off','all');
if isempty(D)
    mdl=fitglm(LP,y,'Distribution','binomial'); score=LP;
else
    mdl=fitglm([LP,D],y,'Distribution','binomial');
    score=mdl.Fitted.LinearPredictor;
end
fit=struct('score',score,'AIC',mdl.ModelCriterion.AIC,'lambda',lambda);
end

function paired=paired_stratified_auc(y,tract_score,classical_score,n_boot,alpha,seed)
y=y(:); tract_score=tract_score(:); classical_score=classical_score(:);
pos=find(y==1); neg=find(y==0);
assert(~isempty(pos)&&~isempty(neg));
rng(seed,'twister');
indices=zeros(numel(y),n_boot,'uint32');
for bi=1:n_boot
    draw_pos=pos(randi(numel(pos),numel(pos),1));
    draw_neg=neg(randi(numel(neg),numel(neg),1));
    indices(:,bi)=uint32([draw_pos;draw_neg]);
end
tract_auc=zeros(n_boot,1); classical_auc=zeros(n_boot,1);
all_two_class=true;
for bi=1:n_boot
    idx=double(indices(:,bi)); yb=y(idx);
    all_two_class=all_two_class && numel(unique(yb))==2;
    tract_auc(bi)=auc_rank(yb,tract_score(idx));
    classical_auc(bi)=auc_rank(yb,classical_score(idx));
end
delta=tract_auc-classical_auc;
paired=struct();
paired.tract=struct('point_auc',auc_rank(y,tract_score), ...
    'mean_bootstrap_auc',mean(tract_auc));
paired.classical=struct('point_auc',auc_rank(y,classical_score), ...
    'mean_bootstrap_auc',mean(classical_auc));
paired.delta_point=paired.tract.point_auc-paired.classical.point_auc;
paired.delta_mean=mean(delta);
paired.delta_ci=quantile(delta,[alpha/2,1-alpha/2]);
paired.index_sha256=sha256_bytes(typecast(indices(:),'uint8'));
paired.all_two_class=all_two_class;
paired.positive_per_resample=numel(pos);
end

function value=auc_rank(y,score)
y=y(:); score=score(:); np=sum(y==1); nn=sum(y==0);
assert(np>0 && nn>0);
ranks=tiedrank(score);
value=(sum(ranks(y==1))-np*(np+1)/2)/(np*nn);
end

function value=hash_ids(ids)
ids=sort(string(ids(:)));
bytes=unicode2native(char(strjoin(ids,newline)),'UTF-8');
value=sha256_bytes(bytes);
end

function value=sha256_bytes(bytes)
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes(:));
raw=typecast(engine.digest(),'uint8');
value=lower(reshape(dec2hex(raw,2).',1,[]));
end
