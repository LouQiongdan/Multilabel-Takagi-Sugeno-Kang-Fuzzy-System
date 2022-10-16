
function recall=Recall(Pre_Labels,test_target)

    [num_class,num_instance]=size(Pre_Labels);
    total=0;
    for i=1:num_instance
        numerator=Pre_Labels(:,i)'*test_target(:,i);
        denominator=test_target(:,i)'*test_target(:,i);
        if denominator ~=0
            total=total+numerator/denominator;
        end
    end
    recall=total/num_instance;
    