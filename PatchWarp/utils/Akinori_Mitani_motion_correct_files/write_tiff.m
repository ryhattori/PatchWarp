function write_tiff(fn,stack,info)
    assert(isnumeric(stack))
    assert(isreal(stack))
    assert(ndims(stack)<=3);
    
    switch class(stack)
        case 'double'
            tagstruct.BitsPerSample=64;
            tagstruct.SampleFormat=Tiff.SampleFormat.IEEEFP;
        case 'single'
            tagstruct.BitsPerSample=32;
            tagstruct.SampleFormat=Tiff.SampleFormat.IEEEFP;
        case 'int8'
            tagstruct.BitsPerSample=8;
            tagstruct.SampleFormat=Tiff.SampleFormat.Int;
        case 'uint8'
            tagstruct.BitsPerSample=8;
            tagstruct.SampleFormat=Tiff.SampleFormat.UInt;
        case 'int16'
            tagstruct.BitsPerSample=16;
            tagstruct.SampleFormat=Tiff.SampleFormat.Int;
        case 'uint16'
            tagstruct.BitsPerSample=16;
            tagstruct.SampleFormat=Tiff.SampleFormat.UInt;
        case 'int32'
            tagstruct.BitsPerSample=32;
            tagstruct.SampleFormat=Tiff.SampleFormat.Int;
        case 'uint32'
            tagstruct.BitsPerSample=32;
            tagstruct.SampleFormat=Tiff.SampleFormat.UInt;
        otherwise
            error('type not supported.');
    end
            
    
    temp_fn = [tempname_if_on_network2(fn)];
    if(~isempty(temp_fn))
        file_to_write = temp_fn;
    else
        file_to_write = fn;
    end
    
    if(nargin<3)
        info = [];
    end
    
    if(iscell(info))
        info=struct('ImageDescription',info); 
    end
    
    tf = [];
    try 
%	try
        tf = Tiff(file_to_write,'w');
%	catch                                      	% added try-catch, but this is not optimal in terms of computation time.  by RH
%	file_to_write = fn;
%	tf = Tiff(file_to_write,'w');
%	end
        
        tagstruct.Photometric=Tiff.Photometric.MinIsBlack;
        tagstruct.Compression=Tiff.Compression.None;
        tagstruct.SamplesPerPixel=1;
        tagstruct.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
        tagstruct.ImageLength=size(stack,1);
        tagstruct.ImageWidth=size(stack,2);
            
        for ii=1:size(stack,3)
            if(numel(info)>=ii && isfield(info,'ImageDescription'))
                tagstruct.ImageDescription=info(ii).ImageDescription;
            elseif(numel(info)==1 && isfield(info,'ImageDescription'))
                tagstruct.ImageDescription=info(1).ImageDescription;
            else
                if(isfield(tagstruct,'ImageDescription'))
                    tagstruct=rmfield(tagstruct,'ImageDescription');
                end
            end
            tf.setTag(tagstruct);
            tf.write(stack(:,:,ii));
            tf.writeDirectory;
        end
        tf.close();
        
        if(~strcmp(file_to_write,fn))
            movefile(file_to_write,fn,'f');
        end
    catch e
        if(isa(tf,'Tiff'));
            tf.close();
        end
        delete(file_to_write);
        rethrow(e);
    end
end