function [unwrPhase, mask, texture, codeWord, ph] = UnwrapPhaseGrayCode_citedi(fileEx, FringePitch, NStep, NBits, MedFilt,format) 
    % Fringe periods in pixels
    % NOTE: fringe periods must be in increasing order and the equivalent
    % fringe period for the first two must be smaller than the longest one
    %FringePitch = 18;
    %NStep = 18;
    %NBits = 7;

    % file folder directory
    %folder = 'ref_lightcrafter';

    %FringeWidth = 640;
    %FringeHeight = 480;
    
    if nargin<5,
        MedFilt = [3 3];
    end

    %Iavg = zeros([ FringeHeight, FringeWidth]);
    % Load three phase-shifted patterns 
    Iavg = 0;
    for Step = 1 : NStep
        fileName = sprintf(['%s%02d' format], fileEx, (Step-1));
%         fileName = sprintf('%s%02d.bmp', fileEx, (Step-1));
        Id = imread(fileName);

        Iavg =Iavg + double(Id)/NStep;
        I(:, :, Step) = double(Id);
    end
    % compute wrapped phase
    [ph, mask, texture] = ComputeNStepPhase(I, 0.1, NStep);



    %Initial codeword to be zeros 
    bcode = zeros(size(Id));
    for nBin = 1:NBits
        % Read binary patterns out
        fileName = sprintf(['%s%02d' format], fileEx, NStep + nBin -1);
%         fileName = sprintf('%s%02d.bmp', fileEx, NStep + nBin - 1);
        I0 = imread(fileName);

        % Binarize the image
        Icode = zeros(size(I0));
        Icode(I0 > Iavg) = 1;

        % Adjust codeword 
        bcode = bcode * 2 + Icode;% 
    end
    
%     figure,imagesc(bcode);

    codeLUT = GenInverseGrayCodeLUT(NBits);

    [ FringeHeight, FringeWidth] = size(bcode);
    % generate codeword by looking up the graycode table
    codeWord = bcode;
    for i = 1: FringeHeight
        for j = 1:FringeWidth 
            codeWord(i, j) = codeLUT(bcode(i, j)+1);
        end
    end

    % map wrapped phase starting from -pi instead of 0
    % NOTE: one pixel phase shift between coded pattern and the sinuosidal
    % pattern
    shift = -pi + pi/FringePitch; 
    ph = atan2(sin(ph  + shift), cos(ph + shift));

    % unwrap the phase value
    uph = ph + codeWord * 2*pi;
    
%     figure,imagesc(codeWord);

    % Shift phase back to the original values
    uph = uph - shift; 

    suph = medfilt2(uph, MedFilt); %modf
    unwrPhase = uph - round((uph - suph)/(2*pi))*2*pi;
end
