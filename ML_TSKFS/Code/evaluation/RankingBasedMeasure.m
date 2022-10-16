
function [Average_Precision,OneError,RankingLoss,Coverage]=RankingBasedMeasure(Outputs,test_target)


    RankingLoss=Ranking_loss(Outputs,test_target);
    OneError=One_error(Outputs,test_target);
    Coverage=coverage(Outputs,test_target);
    %Coverage=Coverage/num_class;
    Average_Precision=Average_precision(Outputs,test_target);
end