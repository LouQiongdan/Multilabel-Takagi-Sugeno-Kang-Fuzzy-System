
function Accuracy=MultiLabelAccuracyEvaluation(predict_target,test_target)

    [num_label,num_test]=size(test_target);
    count=0;
    for i=1:num_test
        numerator=test_target(:,i)'*predict_target(:,i);
        denominator=sum(or(test_target(:,i),predict_target(:,i)));
        if denominator~=0
            count=count + numerator/denominator;
        end
    end
    Accuracy=count/num_test;
end