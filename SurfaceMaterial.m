classdef SurfaceMaterial
    %SURFACEMATERIAL Creates material with given optical properties
    %   Creates material either front or back to the surface with given
    %   indices of refraction for optical lengths of basic d, and 
    %   additional for which chromatic aberration can be minimized e, C', 
    %   F'. In basic example d corresponds to 0.5876 um, e to 0.5461 um, 
    %   C' to 0.6438 um and F' to 0.4800 um according to SCHOTT catalogue
    
    properties(SetAccess = immutable)
        indexd
        indexC
        indexF
        indexe
    end
    
    methods
        function obj = SurfaceMaterial(nd,nC,nF,ne)
            %SURFACEMATERIAL Construct an instance of this class
            %   Set surfaces' indices for basic light length and additional
            %   for which chromatic aberration is minimized 
            obj.indexd = nd;
            obj.indexC = nC;
            obj.indexF = nF;
            obj.indexe = ne;
        end
    end
end

