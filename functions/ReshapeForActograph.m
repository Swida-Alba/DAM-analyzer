function [dataOut, dataOut_R] = ReshapeForActograph(dataIn,options)
    arguments
        dataIn
        options.BlockNumPerDay = 48; % block number per day, 24 h / binWidth
    end
    dayNum = ceil(length(dataIn) / options.BlockNumPerDay);
    dataMat = nan(1,dayNum*options.BlockNumPerDay);
    dataMat(1:length(dataIn)) = dataIn;
    dataMat = reshape(dataMat,options.BlockNumPerDay,dayNum)';
    dataOut_L = [nan(1,options.BlockNumPerDay); dataMat];
    dataOut_R = [dataMat; nan(1,options.BlockNumPerDay)];
    dataOut = [dataOut_L dataOut_R];
    % dataOut = dataOut(1:end-1,:);
end