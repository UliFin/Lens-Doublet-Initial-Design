classdef OpticSurface
    %OPTICSURFACE
    %   Class of optical surface between two materials: in front of 
    %   and behind the surface. Prescribes aberrations, geometrical 
    %   refraction parameters of beams from given point on optical axis 
    %   in given paraxial height, geometry of surface.

    properties
        ddnn % Helping value
        dnp % Helping value
        curID % Index of optic surface

        % Materials
        frontMaterial 
        backMaterial

        %Refraction parameters
        objectPosition
        imagePosition
        paraxialHeight
        
        %Aberrations
        abbeInvariant
        chromaticAberrationCoefficient
        sphericalAberrationCoefficient
        
        % Geometry
        curveRadius % Curvature radius
        previousSurfaceDistance
    end

    methods
        function obj = OpticSurface(ndfront, nCfront, nFfront, nefront, ...
                ndback, nCback, nFback, neback, previousParaxialHeight, previousImagePosition, varargin) % indices

            %OPTICSURFACE Construct an instance of this class
            %   Creates optic surface object with given optical and 
            %   geometry properties

            %Setting geometry(optimization) parameters to symbolical
            obj.curID = OpticSurface.getID;
            strId = num2str(obj.curID);
            obj.curveRadius = sym(strcat('curveRadius', strId));
            obj.previousSurfaceDistance = sym(strcat('previousSurfaceDistance', strId));
            
            % Creating materials 
            obj.frontMaterial = SurfaceMaterial(ndfront,nCfront,nFfront,nefront);
            obj.backMaterial = SurfaceMaterial(ndback,nCback,nFback,neback);

            % For first surface paraxial height and object position is
            % directly set
            if (obj.curID == 1)
                obj.paraxialHeight = previousParaxialHeight;
                obj.objectPosition = previousImagePosition;
            else
                % Example: object position for second surface -> p2=p11-d1
                obj.objectPosition = previousImagePosition - obj.previousSurfaceDistance;
                % Example: Paraxial height for second surface -> h2=h1*p2/p11
                obj.paraxialHeight = previousParaxialHeight*obj.objectPosition/previousImagePosition;
            end

            obj = calculateall(obj);
        end
% Getters
        function out = get.curveRadius(obj)
            out = obj.curveRadius;
        end
        function out = get.previousSurfaceDistance(obj)
            out = obj.previousSurfaceDistance;
        end
       
        function obj = calculateall(obj)
        % Setting symbolical equations to surface parameters
            obj = calculateddnn(obj);
            obj = calculatesymparameters(obj);
        end
        function obj = calculateddnn(obj)
            % ddnn - Helping variable
            % Example: for first surface-> ddnn1 = (n2ff-n2cc)/n2e-(n1ff-n1cc)/n1e
            obj.ddnn = (obj.backMaterial.indexF-obj.backMaterial.indexC)/obj.backMaterial.indexe-(obj.frontMaterial.indexF-obj.frontMaterial.indexC)/obj.frontMaterial.indexe;
        end

        function obj = calculatesymparameters(obj)
        % Setting numerical equations results to surface parameters
        % depending on curvature radii and thickness of the lens (previous
        % surface position)
       
            %Example: Image position for first surface-> p11 = n2/((n2-n1)/R1+n1/(p1));
            obj.imagePosition = obj.backMaterial.indexd/((obj.backMaterial.indexd-obj.frontMaterial.indexd)/obj.curveRadius + obj.frontMaterial.indexd/obj.objectPosition);

            %Example:  Abbe Invariant for first surface-> Q1 = n1*(1/R1-1/p1)
            obj.abbeInvariant = obj.frontMaterial.indexd*(1/obj.curveRadius-1/obj.objectPosition);

            %Example:  dnp - Helping variable for first surface-> dnp1 = 1/(n2*p11)-1/(n1*p1)
            obj.dnp = 1/(obj.backMaterial.indexd*obj.imagePosition)-1/(obj.frontMaterial.indexd*obj.objectPosition);

            %Example:  sphericalAberrationCoefficient for first surface-> S1 = h1^4*Q1^2*dnp1
            obj.sphericalAberrationCoefficient=obj.paraxialHeight^4*obj.abbeInvariant^2*obj.dnp;

            %Example:  chromaticAberrationCoefficient for first surface-> E1 = h1^2*Q1*((n2ff-n2cc)/n2e-(n1ff-n1cc)/n1e)
            obj.chromaticAberrationCoefficient = obj.paraxialHeight^2*obj.abbeInvariant*obj.ddnn;

        end

        function obj = updateValues(obj, cr, psd, previousImagePosition, previousParaxialHeight, varargin)
        % Update curvature radius and previous surface distance and
        % calculate depending on them parameters of surface
            obj.curveRadius = cr;
            obj.previousSurfaceDistance = psd;
            % obj.calculateall();
            %   NEW
            if nargin > 3
            
                % Example: object position for second surface -> p2=p11-d1
                obj.objectPosition = previousImagePosition - obj.previousSurfaceDistance;
                % Example: Paraxial height for second surface -> h2=h1*p2/p11
                obj.paraxialHeight = previousParaxialHeight*obj.objectPosition/previousImagePosition;
            else
            end
            obj = calculateall(obj);
        end
    end

    methods (Static)
        function out = getID()
            persistent id;
            if isempty(id)
                id = 1;
            end
            out = id;
            id = out + 1;
        end
    end
end


