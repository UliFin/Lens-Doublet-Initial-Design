classdef LensSystemDesign
    %LENSSYSTEMDESIGN Lens systems class
    %   Creates optical system with four optical surfaces -  two lenses in
    %   air with prescribed parameters - optical, mechanical. 

    properties

        % System parameters 
        EFL; % System effective focal length
        H1; % Half of effective aperture
        Lmax;   % Maximum length of the system - vertex to vertex of first
        % and last optical surface
        magnification = sym('magnification');
        systemSphericalAberrationCoefficient = sym('systemSphericalAberrationCoefficient');
        systemSphericalAberration = sym('systemSphericalAberration');
        systemChromaticAberrationCoefficient = sym('systemChromaticAberrationCoefficient');
        systemChromaticAberration = sym('systemChromaticAberration');
        systemCalculatedEFL
        % Optic Surfaces
        optSurf1 
        optSurf2 
        optSurf3 
        optSurf4 

    end

    methods
        function obj = LensSystemDesign(systemEFL,multiplier,halfEfectiveAperture,refIndAir, ...
                refIndGlass1,refIndGlass2,refIndAirF,refIndAirC, ...
                refIndAire,refIndGlass1F,refIndGlass1C,refIndGlass1e, ...
                refIndGlass2F,refIndGlass2C,refIndGlass2e,firstParaxialHeight, ...
                maximumSysLength)
            %LENSSYSTEMDESIGN Construct an instance of this class
            %   Creates lens system object consisting of four optic
            %   surfaces - two lenses, prescribes symbolical parameters of
            %   the system depending on optical surfaces parameters.

            %Setting system parameters
            obj.EFL = systemEFL;
            obj.Lmax = maximumSysLength;
            obj.H1 = halfEfectiveAperture;
            firstObjectPosition = obj.EFL*multiplier;

            % Create optical surfaces of two lenses
            obj.optSurf1 = OpticSurface(refIndAir, refIndAirC, refIndAirF, refIndAire, refIndGlass1,refIndGlass1C, refIndGlass1F, refIndGlass1e, firstParaxialHeight, firstObjectPosition,1);
            obj.optSurf2 = OpticSurface(refIndGlass1, refIndGlass1C,refIndGlass1F, refIndGlass1e, refIndAir, refIndAirC, refIndAirF, refIndAire, obj.optSurf1.paraxialHeight, obj.optSurf1.imagePosition);
            obj.optSurf3 = OpticSurface(refIndAir, refIndAirC, refIndAirF, refIndAire, refIndGlass2, refIndGlass2C, refIndGlass2F, refIndGlass2e, obj.optSurf2.paraxialHeight, obj.optSurf2.imagePosition);
            obj.optSurf4 = OpticSurface(refIndGlass2, refIndGlass2C, refIndGlass2F, refIndGlass2e, refIndAir, refIndAirC, refIndAirF, refIndAire, obj.optSurf3.paraxialHeight, obj.optSurf3.imagePosition);

            obj = obj.CalculateSystemValues();
        end

        function obj = CalculateSystemValues(obj)
            obj.magnification = obj.optSurf1.frontMaterial.indexd/obj.optSurf4.backMaterial.indexd*obj.optSurf1.paraxialHeight/obj.optSurf4.paraxialHeight*obj.optSurf4.imagePosition/obj.optSurf1.objectPosition;
            obj.systemSphericalAberrationCoefficient = obj.optSurf1.sphericalAberrationCoefficient+obj.optSurf2.sphericalAberrationCoefficient+obj.optSurf3.sphericalAberrationCoefficient+obj.optSurf4.sphericalAberrationCoefficient;
            obj.systemChromaticAberrationCoefficient = obj.optSurf1.chromaticAberrationCoefficient+obj.optSurf2.chromaticAberrationCoefficient+obj.optSurf3.chromaticAberrationCoefficient+obj.optSurf4.chromaticAberrationCoefficient;
            obj.systemSphericalAberration = -0.5*(obj.H1*obj.magnification*obj.optSurf1.objectPosition)^2*obj.systemSphericalAberrationCoefficient;
            obj.systemChromaticAberration = -(obj.magnification*obj.optSurf1.objectPosition)^2*obj.systemChromaticAberrationCoefficient;
        end

        function out = GAfunction(obj, A, b, lb, ub, gaOptions, ksph, kchr,varargin)
            % Genetic algorithm optimization. For more information see
            % Matlab documentation for ga function. Optimization parameters
            % are optical surface curvature radii and lens thicknesses
            % resp. distances between lenses.
     if nargin > 7
     elseif nargin > 6
     kchr=1;
     else
     ksph=1;kchr=1;
     end

            % number of variables - optimization parameters
            numbvars = 7;
            %Linear constraints Ax=b, x - parameters vector in alphabetical
            %order. Here [curRad1 curRad2 curRad3 curRad4 thick2 thick3 thick4]
            Aeq=[];
            beq=[];

            % Preparing merit (aim) function for minimization - Here
            % minimizing spherical and chromatic aberrations for light
            % length (d) and (e), (F'), (C') accordingly.
            % (d)=0.5876µm, (F')=0.480µm, (C')=0.6438µm, (e)=0.5461µm
            fun = matlabFunction(ksph*abs(obj.systemSphericalAberration) + kchr*abs(obj.systemChromaticAberration));

            out = ga(@(x)fun(x(1),x(2),x(3),x(4),x(5),x(6),x(7)), ...
                numbvars,A,b,Aeq,beq,lb,ub, ...
                @(x)ConstraintFunctions.simple_constraint(x,obj.H1,obj.EFL,obj.optSurf1.backMaterial.indexd,obj.optSurf3.backMaterial.indexd, obj.optSurf1.frontMaterial.indexd), ...
                gaOptions);
        end

        function obj = inserValues(obj, values)
            
            % Renew values of curvature radius and previous surface 
            % distance* to newly found ones, from results of optimization.
            % *Note: First surface has no previous surface, so it is 0 mm
            % Values is output vector of ga optimization, so it is
            % corresponding to x and is a vector with optimization
            % parameters in alphabetical order. 
            % Here [curRad1 curRad2 curRad3 curRad4 thick2 thick3 thick4]
            R1=values(1); % First lens first surface curvature
            R2=values(2); % First lens second surface curvature
            R3=values(3); % Second lens first surface curvature
            R4=values(4); % Second lens second surface curvature
            d1=values(5); % Thickness of first lens on optic axis
            d=values(6); % Thickness of second lens on optic axis
            d2=values(7);  % Distance between lenses inbetween on optic axis

            obj.optSurf1 = obj.optSurf1.updateValues(R1, 0);
            obj.optSurf2 = obj.optSurf2.updateValues(R2, d1, obj.optSurf1.imagePosition, obj.optSurf1.paraxialHeight);
            obj.optSurf3 = obj.optSurf3.updateValues(R3, d, obj.optSurf2.imagePosition, obj.optSurf2.paraxialHeight);
            obj.optSurf4 = obj.optSurf4.updateValues(R4, d2, obj.optSurf3.imagePosition, obj.optSurf3.paraxialHeight);
            obj = obj.CalculateSystemValues();


            n2=obj.optSurf1.backMaterial.indexd;
            n1=obj.optSurf1.frontMaterial.indexd;
            n3=obj.optSurf3.backMaterial.indexd;

% Lens maker equation for thick lens, Focal length of first lens  
            f11=1/((n2/n1-1)*(1/R1-1/R2)+(n2/n1-1)^2*d1/(n2*R1*R2)); 
            
            % Lens maker equation for thick lens, Focal length of second lens
            f22=1/((n3/n1-1)*(1/R3-1/R4)+(n3/n1-1)^2*d2/(n3*R3*R4));

            % Two lenses system eqivalent focal length
obj.systemCalculatedEFL = 1/(1/f11+1/f22-d/(f11*f22*n1)...
                -(n2-n1)*d1/(n2*R1*f22)...
                +(n3-n1)*d2/(n3*R4*f11));
            % 
            % 
            % % Lens maker equation for thick lens, Focal length of second lens
            % f11=1/((obj.optSurf1.backMaterial.indexd/obj.optSurf1.frontMaterial.indexd-1)*(1/values(1)-1/(values(2)))+(obj.optSurf1.backMaterial.indexd/obj.optSurf1.frontMaterial.indexd-1)^2*values(5)/(obj.optSurf1.backMaterial.indexd*values(1)*values(2)));          
            % % Lens maker equation for thick lens, Focal length of second lens
            % f22=1/((obj.optSurf3.backMaterial.indexd/obj.optSurf3.frontMaterial.indexd-1)*(1/values(3)-1/(values(4)))+(obj.optSurf3.backMaterial.indexd/obj.optSurf3.frontMaterial.indexd-1)^2*values(6)/(obj.optSurf3.backMaterial.indexd*values(3)*values(4)));
            % 
            % obj.systemCalculatedEFL = 1/(1/f11+1/f22-values(6)/(f11*f22*obj.optSurf1.frontMaterial.indexd)...
            %     -(obj.optSurf1.backMaterial.indexd-obj.optSurf1.frontMaterial.indexd)*values(5)/(obj.optSurf1.backMaterial.indexd*values(1)*f22)...
            %     +(obj.optSurf3.backMaterial.indexd-obj.optSurf3.frontMaterial.indexd)*values(6)/(obj.optSurf3.backMaterial.indexd*values(4)*f11));

        end

        function visualiseSystem(obj)
            %% VISUALIZATION

            R1 = obj.optSurf1.curveRadius;
            R2 = obj.optSurf2.curveRadius;
            R3 = obj.optSurf3.curveRadius;
            R4 = obj.optSurf4.curveRadius;
            d1 = obj.optSurf2.previousSurfaceDistance;
            d2 = obj.optSurf4.previousSurfaceDistance;
            d = obj.optSurf3.previousSurfaceDistance;
            H1 = obj.H1;

            XC1=-d1-d/2+R1; % Center of first curvature
            XC2=-d/2+R2; % Center of second curvature
            XC3=d/2+R3; % Center of third curvature
            XC4=d/2+d2+R4; % Center of fourth curvature
            
            % Intersections of circles
            XCC1=circcirc(XC1,0,abs(R1),XC2,0,abs(R2));
            XCC2=circcirc(XC3,0,abs(R3),XC4,0,abs(R4));
            xcc1=XCC1(1);
            xcc2=XCC1(1);
            xcc3=XCC2(1);
            xcc4=XCC2(1);
            
            % Intersection with mechanical aperture - horizontal line
            XV1 = linecirc(0,0,XC1,0,abs(R1));
            XV2 = linecirc(0,0,XC2,0,abs(R2));
            XV3 = linecirc(0,0,XC3,0,abs(R3));
            XV4= linecirc(0,0,XC4,0,abs(R4));

            % Finding the limits for graphical image
            if R1>0
                xv1=min(XV1(1), XV1(2));
            else
                xv1=max(XV1(1), XV1(2));
            end
            if R2>0
                xv2=min(XV2(1), XV2(2));
            else
                xv2=max(XV2(1), XV2(2));
            end
            if R3>0
                xv3=min(XV3(1), XV3(2));
            else
                xv3=max(XV3(1), XV3(2));
            end
            if R4>0
                xv4=min(XV4(1), XV4(2));
            else
                xv4=max(XV4(1), XV4(2));
            end
            if isnan(xcc1)&&R1<0
                xcc1=max(linecirc(0,H1,XC1,0,abs(R1)));
            elseif isnan(xcc1)
                xcc1=min(linecirc(0,H1,XC1,0,abs(R1)));
            end
            if isnan(xcc2)&&R2<0
                xcc2=max(linecirc(0,H1,XC2,0,abs(R2)));
            elseif isnan(xcc2)
                xcc2=min(linecirc(0,H1,XC2,0,abs(R2)));
            end
            if isnan(xcc3)&&R3<0
                xcc3=max(linecirc(0,H1,XC3,0,abs(R3)));
            elseif isnan(xcc3)
                xcc3=min(linecirc(0,H1,XC3,0,abs(R3)));
            end
            if isnan(xcc4)&&R4<0
                xcc4=max(linecirc(0,H1,XC4,0,abs(R4)));
            elseif isnan(xcc4)
                xcc4=min(linecirc(0,H1,XC4,0,abs(R4)));
            end

            % Graphics of curvatures inbetween the limits vertical and
            % horizontal
            % First lens
            fplot(@(x) sqrt(R1^2-(x-XC1)^2),[min(xcc1,xv1) max(xcc1,xv1)],'r')
            hold on
            fplot(@(x) -sqrt(R1^2-(x-XC1)^2),[min(xcc1,xv1) max(xcc1,xv1)],'r')
            hold on
            fplot(@(x) sqrt(R2^2-(x-XC2)^2),[min(xcc2,xv2) max(xcc2,xv2)],'r')
            hold on
            fplot(@(x) -sqrt(R2^2-(x-XC2)^2),[min(xcc2,xv2) max(xcc2,xv2)],'r')
            hold on
            % Second lens
            fplot(@(x) sqrt(R3^2-(x-XC3)^2),[min(xcc3,xv3) max(xcc3,xv3)],'b')
            hold on
            fplot(@(x) -sqrt(R3^2-(x-XC3)^2),[min(xcc3,xv3) max(xcc3,xv3)],'b')
            hold on
            fplot(@(x) sqrt(R4^2-(x-XC4)^2),[min(xcc4,xv4) max(xcc4,xv4)],'b')
            hold on
            fplot(@(x) -sqrt(R4^2-(x-XC4)^2),[min(xcc4,xv4) max(xcc4,xv4)],'b')
            hold on
            MatrixMin=[xcc1,xcc2,xcc3,xcc4,xv4,xv3,xv2,xv1];
            % Optical Axis
            fplot(@(x) 0*x,[min(MatrixMin) max(MatrixMin)],'--ok')
            hold on
            % Vertical limits of aperture
            fplot(@(x) 10+0*x,[min(MatrixMin) max(MatrixMin)],'--k')
            hold on
            fplot(@(x) -10+0*x,[min(MatrixMin) max(MatrixMin)],'--k')
            hold on
            txt = '\uparrow Optical Axis';
            text(0,0-1,txt)
            txt1 = '\uparrow First Mechanical Aperture';
            text(0,H1-1,txt1)
            hold on
            grid on
            pbaspect([1 1 1])
            xlim([min(MatrixMin)-abs(min(MatrixMin))/5 max(MatrixMin)+max(MatrixMin)/5])
            ylim([-H1-H1/2 H1+H1/2])
            annotation('textbox', [0.83, 0.1, 0.1, 0.1], 'String', "d2 = " + sprintf('%.2f', d2) +" mm")
            annotation('textbox', [0.83, 0.2, 0.1, 0.1], 'String', "d = " + sprintf('%.2f', d)+" mm")
            annotation('textbox', [0.83, 0.3, 0.1, 0.1], 'String', "d1 = " + sprintf('%.2f', d1)+" mm")
            hold off
        end



        function out = GAfunction1(obj, A, b, lb, ub, gaOptions, ksph, kchr,varargin)
            % Genetic algorithm optimization. For more information see
            % Matlab documentation for ga function. Optimization parameters
            % are optical surface curvature radii and lens thicknesses
            % resp. distances between lenses.
     if nargin > 7
     elseif nargin > 6
     kchr=1;
     else
     ksph=1;kchr=1;
     end

            % number of variables - optimization parameters
            numbvars = 7;
            %Linear constraints Ax=b, x - parameters vector in alphabetical
            %order. Here [curRad1 curRad2 curRad3 curRad4 thick2 thick3 thick4]
            Aeq=[];
            beq=[];

            % Preparing merit (aim) function for minimization - Here
            % minimizing spherical and chromatic aberrations for light
            % length (d) and (e), (F'), (C') accordingly.
            % (d)=0.5876µm, (F')=0.480µm, (C')=0.6438µm, (e)=0.5461µm
            fun = matlabFunction(ksph*abs(obj.systemSphericalAberration) + kchr*abs(obj.systemChromaticAberration));

            out = ga(@(x)fun(x(1),x(2),x(3),x(4),x(5),x(6),x(7)), ...
                numbvars,A,b,Aeq,beq,lb,ub, ...
                @(x)ConstraintFunctions1.simple_constraint(x,obj.H1,obj.EFL,obj.optSurf1.backMaterial.indexd,obj.optSurf3.backMaterial.indexd, obj.optSurf1.frontMaterial.indexd), ...
                gaOptions);
        end





    end
end

