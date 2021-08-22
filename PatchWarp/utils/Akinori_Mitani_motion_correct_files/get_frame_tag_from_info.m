function frame_tag = get_frame_tag_from_info(info)
    frame_tag = NaN(size(info));
    for i=1:numel(info)
        if(isfield(info,'ImageDescription'))
            tmp = str2double(regexp( ...
                    info(i).ImageDescription,...
                    'Frame Tag = (\d*)','tokens','once'));
            if(~isempty(tmp))
                frame_tag(i) = tmp;
            end
        end
    end
end