
function [ BestParameter, BestResult ] = ML_TSKFS_adaptive_validate( data, target, oldOptmParameter,TSKoptions)
optmParameter         = oldOptmParameter;
alpha_searchrange     = oldOptmParameter.alpha_searchrange;
beta_searchrange      = oldOptmParameter.beta_searchrange;
gamma_searchrange     = oldOptmParameter.gamma_searchrange;

k_searchrange=TSKoptions.k_searchrange;
h_searchrange=TSKoptions.h_searchrange;
total = length(alpha_searchrange)*length(beta_searchrange)*length(gamma_searchrange)*length(k_searchrange)*length(h_searchrange);
index = 1;
parameter_cell=zeros(total,35);
ii=1;
for p=1:length(k_searchrange)
    for q=1:length(h_searchrange)
        TSKoptions.k=k_searchrange(p);
        TSKoptions.h=h_searchrange(q);
        [v,b] = gene_ante_fcm(data,TSKoptions);
        [G_data] = calc_x_g(data,v,b);
        
        train_data=G_data;
        num_train = size(train_data,1);
        randorder = randperm(num_train);
        
        BestResult = zeros(15,1);
        num_cv = 5;
   
        
        for i=1:length(alpha_searchrange) 
            for j=1:length(beta_searchrange) 
                for k = 1:length(gamma_searchrange) 
                    fprintf('\n-   %d-th/%d: search parameter,TSK_k= %f,TSK_h= %f, alpha = %f, beta = %f, and gamma = %f',index, total,k_searchrange(p),h_searchrange(q), alpha_searchrange(i), beta_searchrange(j), gamma_searchrange(k));
                    index = index + 1;
                    optmParameter.alpha   = alpha_searchrange(i);
                    optmParameter.beta    = beta_searchrange(j); 
                    optmParameter.gamma   = gamma_searchrange(k);
                    
                    optmParameter.maxIter           = 100;
                    optmParameter.minimumLossMargin = 0.01;
                    optmParameter.outputtempresult  = 0;
                    optmParameter.drawConvergence   = 0;
                    
                    Result = zeros(15,1);
                    cv_index=1;
                    TempResult=zeros(num_cv,15);
                    
                    for cv = 1:num_cv
                        [cv_train_data,cv_train_target,cv_test_data,cv_test_target ] = generateCVSet( train_data,target',randorder,cv,num_cv);
                        [model_LLSF]  = ML_TSKFS( cv_train_data, cv_train_target,optmParameter);
                        Outputs     = (cv_test_data*model_LLSF)';
                        Pre_Labels  = round(Outputs);
                        Pre_Labels  = (Pre_Labels >= 1); 
                        Pre_Labels  = double(Pre_Labels);
                        TempResult(cv_index,:)=EvaluationAll(Pre_Labels,Outputs,cv_test_target')';
                        cv_index=cv_index+1;
                    end
                    Result=(mean(TempResult))';
                    STD=std(TempResult);
                    if optmParameter.bQuiet == 0
                        parameter_cell(ii,1:15)=Result(1:15,1)';
                        parameter_cell(ii,16:30)=STD(1,1:15);
                        parameter_cell(ii,31)=alpha_searchrange(i);
                        parameter_cell(ii,32)=beta_searchrange(j);
                        parameter_cell(ii,33)=gamma_searchrange(k);
                        parameter_cell(ii,34)=k_searchrange(p);
                        parameter_cell(ii,35)=h_searchrange(q);
                        ii=ii+1;
                        save('Parameter_cell','parameter_cell')
                    end
                    r = IsBetterThanBefore(BestResult,Result);
                    if r == 1
                        BestResult = Result;
                        PrintResults(Result);
                        BestParameter.Optm_Parameter = optmParameter;
                        BestParameter.TSK_options=TSKoptions;
                    end
                end
            end
        end
    end
end

end


function r = IsBetterThanBefore(Result,CurrentResult)
% 1 HammingLoss
% 2 ExampleBasedAPCCuracy
% 3 ExampleBasedPrecision
% 4 ExampleBasedRecall
% 5 ExampleBasedFmeasure
% 6 SubsetAPCCuracy
% 7 LabelBasedAPCCuracy
% 8 LabelBasedPrecision
% 9 LabelBasedRecall
% 10 LabelBasedFmeasure
% 11 MicroF1Measure
% 12 Average_Precision
% 13 OneError
% 14 RankingLoss
% 15 Coverage
%
%  the combination of Accuracy, F1, Macro F1 and Micro F1. Of course, any evaluation metrics or the combination of them can be used.

a = CurrentResult(2,1) + CurrentResult(5,1)  + CurrentResult(10,1) + CurrentResult(11,1);
b = Result(2,1) + Result(5,1) + Result(10,1) + Result(11,1);

if a > b
    r =1;
else
    r = 0;
end
end
