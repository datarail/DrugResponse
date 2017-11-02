function CyCIFdata = Merge_CyCIFcycles(SingleCelldata)
% CyCIFdata = Merge_CyCIFcycles(SingleCelldata)

Commonfields = {'Barcode' 'Well' 'Date' 'Background'};
datafields = setdiff(fieldnames(SingleCelldata), Commonfields);
datafields = datafields(strfindcell(datafields, 'X_c')~=1);
datafields = datafields(strfindcell(datafields, 'Y_c')~=1);
datafields = datafields(strfindcell(datafields, 'Field_c')~=1);

Ncycle = max(cellstr2mat(regexpcelltokens(datafields, 'cycle([0-9]*)$')));

clear CyCIFdata

fprintf('\n---------------------\n  CyCIF merging:\n\n');
for iW = 1:length(SingleCelldata);
    fprintf('Well %i/%i: ', iW, length(SingleCelldata));

    pos = cell(1,Ncycle);
    for iC=1:Ncycle
        pos{iC} = [SingleCelldata(iW).(['X_cycle' num2str(iC)]) ...
            SingleCelldata(iW).(['Y_cycle' num2str(iC)]) ...
            SingleCelldata(iW).(['Field_cycle' num2str(iC)]) ...
            NaN(length(SingleCelldata(iW).(['Y_cycle' num2str(iC)])), 1)];
    end

    maxpos = cellfun2(@(x) max(x(:,[1 2]))', pos);
    maxpos = max([maxpos{:}]')+5;

    %%
    Nfields = unique(pos{1}(:,3));
    for iF = Nfields'
        fprintf('F%i(1', iF);
        Fidx = iF/(10^ceil(log10(length(Nfields)+1)));
        %%
        ref_pos = pos{1}(pos{1}(:,3)==iF,[1 2]);

        ref_mask = zeros(maxpos);
        ref_mask(sub2ind(maxpos, ref_pos(:,1), ref_pos(:,2))) = 1;
        H = fspecial('disk',2);
        ref_mask = imfilter(ref_mask,H,'symmetric');
        H = fspecial('gaussian',10);
        ref_mask = imfilter(ref_mask,H,'symmetric');

        pos{1}(pos{1}(:,3)==iF,4) = (1:sum(pos{1}(:,3)==iF)) + Fidx;

        for iC=2:Ncycle
            fprintf('|%i', iC);
            %%
            test_pos = pos{iC}(pos{iC}(:,3)==iF,[1 2]);

            test_mask = zeros(maxpos);
            test_mask(sub2ind(maxpos, test_pos(:,1), test_pos(:,2))) = 1;
            H = fspecial('disk',2);
            test_mask = imfilter(test_mask,H,'symmetric');
            H = fspecial('gaussian',10);
            test_mask = imfilter(test_mask,H,'symmetric');

            [optimizer, metric] = imregconfig('monomodal');
            optimizer.GradientMagnitudeTolerance = 2e-6;
            optimizer.MinimumStepLength = 1e-4;
            optimizer.MaximumStepLength = .5;
            optimizer.MaximumIterations = 300;
            tform = imregtform(test_mask, ref_mask, 'translation', optimizer, metric);

            %
            test_pos_corr = test_pos;
            test_pos_corr(:,[2 1]) = max(round(tform.transformPointsForward(test_pos(:,[2 1]))),1);
            test_mask_corr = zeros(maxpos);
            test_mask_corr(sub2ind(maxpos, min(test_pos_corr(:,1),maxpos(1)), ...
                min(test_pos_corr(:,2),maxpos(2)))) = 1;
            H = fspecial('disk',4);
            test_mask_corr = imfilter(test_mask_corr,H,'symmetric');

%             figure(4)
%             imshowpair(test_mask_corr, ref_mask, 'Scaling','joint');
%             pause

            d2 = pdist2(ref_pos, test_pos_corr);
            % hard cutoff
            d2(d2>20) = Inf;
            %%
            ref2test = argmin(d2,[],2);
            test2ref = argmin(d2,[],1);
            test2ref_idx = NaN*test2ref';
            for i=1:length(test2ref)
                if d2(test2ref(i),i)<Inf && ref2test(test2ref(i))==i
                    test2ref_idx(i) = test2ref(i) + Fidx;
                end
            end

            pos{iC}(pos{iC}(:,3)==iF,4) = test2ref_idx;

        end
        fprintf('), ');
    end
    %%

    pos_intersect = pos{1}(:,4);
    for iC=2:Ncycle
        pos_intersect = intersect(pos_intersect, pos{iC}(:,4));
    end
    idxs = NaN(length(pos_intersect), Ncycle);
    for iC=1:Ncycle
        temp = memberidx(pos_intersect, pos{iC}(:,4));
        idxs(:,iC) = temp(temp~=0);
    end

    temp = struct();
    for i=1:length(Commonfields)
        temp.(Commonfields{i}) = SingleCelldata(iW).(Commonfields{i});
    end

    temp.data = table(pos{1}(idxs(:,1),3),'variableName',{'Field'});
    for i=1:length(datafields)
        iC = cellstr2mat(regexpcelltokens(datafields(i), 'cycle([0-9]*)$'));
        temp.data = [temp.data table(SingleCelldata(iW).(datafields{i})(idxs(:,iC)), ...
            'variablename', datafields(i))];
    end

    CyCIFdata(iW) = temp;
    fprintf(' done\n');
end
