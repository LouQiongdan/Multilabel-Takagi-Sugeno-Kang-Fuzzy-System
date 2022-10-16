
function Accuracy=AccuracyEvaluation(predict_target,test_target)

    num_test=size(test_target,2);
    correctones=(predict_target==test_target);
    Accuracy=sum(correctones)/num_test;

end