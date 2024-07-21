function colorPatch = AddPatchMap(colorPatch, options)
    arguments
        colorPatch;
        options.BlockNumPerDay = 48;
    end

    for i = 1:length(colorPatch)
        PatchMap = [];
        protocol_len = size(colorPatch{i}.ProtocolMat,2);
        for j = 1:protocol_len
            tmap = repmat(colorPatch{i}.ProtocolMat(2,j),1,colorPatch{i}.ProtocolMat(1,j)*2);
            PatchMap = cat(2,PatchMap,tmap);
        end
        [PatchMap, HalfMap] = ReshapeForActograph(PatchMap, "BlockNumPerDay", options.BlockNumPerDay);
        colorPatch{i}.PatchMap = PatchMap;
        colorPatch{i}.HalfMap = HalfMap;
    end
end