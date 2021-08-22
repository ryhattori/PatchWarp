classdef ScanImageTiffReaderTests < matlab.unittest.TestCase    
    % Test suite for the ScanImageTiffReader
    %       
    % EXAMPLE
    %   runtests('ScanImageTiffReader');
    %
    % Copyright 2016-2018 Vidrio Technologies, LLC
    %
    % Licensed under the Apache License, Version 2.0 (the "License");
    % you may not use this file except in compliance with the License.
    % You may obtain a copy of the License at
    % 
    %     http://www.apache.org/licenses/LICENSE-2.0
    % 
    % Unless required by applicable law or agreed to in writing, software
    % distributed under the License is distributed on an "AS IS" BASIS,
    % WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    % See the License for the specific language governing permissions and
    % limitations under the License.
    
    methods(Test)
        function apiVersionStringIsWellFormed(obj)
            import ScanImageTiffReader.ScanImageTiffReader;
            import matlab.unittest.constraints.ContainsSubstring;
            % ensure we didn't have any problems with macro expansion
            obj.verifyThat(...
                ScanImageTiffReader.apiVersion(),...
                ~ContainsSubstring('GIT_TAG'));
        end


        function basicOperations(~)
            % Just making the calls to check if anything crashes/throws.
            % Currently, we rely on other testing of the API to make sure outputs
            % are what we expect.
            import ScanImageTiffReader.ScanImageTiffReader;
            reader=ScanImageTiffReader('./data/resj_00001.tif');
            reader.data();
            reader.descriptions();
            reader.metadata();
            reader.apiVersion();
        end
        
        function singleImage(obj)
            import ScanImageTiffReader.ScanImageTiffReader;
            reader=ScanImageTiffReader('./data/single_image_si.tif');
            obj.verifyEqual(numel(reader.descriptions()),1);
        end
    end
end
