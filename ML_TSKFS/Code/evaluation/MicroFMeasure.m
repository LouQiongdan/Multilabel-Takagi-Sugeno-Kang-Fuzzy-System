
function MicroF1Measure=MicroFMeasure(test_targets,predict_targets)

    test_targets=double(test_targets==1);
    predict_targets=double(predict_targets==1);
    [L,num_test]=size(test_targets);
    groundtruth=reshape(test_targets,1,L*num_test);
    predict=reshape(predict_targets,1,L*num_test);
    intersection=groundtruth*predict';
    precision = intersection/(sum(predict)+eps);
    recall = intersection/(sum(groundtruth)+eps);
    MicroF1Measure=2*precision*recall/(precision+recall+eps);
    
end