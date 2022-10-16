
function SubsetAccuracy=SubsetAccuracyEvaluation(test_target,predict_target)

    [~,num_test]=size(test_target);
    count=0;
    for i=1:num_test
        if isequal(test_target(:,i),predict_target(:,i))==1
            count=count+1;
        end
    end
    SubsetAccuracy=count/num_test;
end