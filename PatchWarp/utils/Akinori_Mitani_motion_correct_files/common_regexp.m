function expression = common_regexp(opt)
    switch(opt)
        case 'tiff_ext'
            expression = '\.[tT][iI][fF][fF]?$';
        case 'mat_ext'
            expression = '\.[mM][aA][tT]$';
        otherwise
            warning('Unknown option');
            expression = '';
    end
end
            