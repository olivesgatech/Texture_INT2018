function descriptor = CLBP_DV(img, numBins)

radSample = [2, 16]; % the best performing out of all CLBP params
mapping=getmapping(radSample(2),'riu2');

descriptor = clbpJoint(img,radSample(1),radSample(2),mapping, 'nh'); %LBP histogram in (R, N) neighborhood

