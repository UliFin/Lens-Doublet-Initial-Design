classdef ConstraintFunctions
    %CONSTRAINTFUNCTIONS Creating custom constraint function for ga
    %   optimization function
    %   For more information see GA constraints and options in 
    %   Matlab documentation 

    methods (Static)
        function [C, Ceq] = simple_constraint(x,H1,EFL,n2,n3,n1)
            %simple_constraint Custom constraint function for ga
            %   optimization. More information in Matlab documentation for
            %   ga function. Here, mechanical constraint of lens edge 
            %   thickness, lenses edges not collapsing* into each other and 
            %   effective focal length of the system is prescribed.
            %   Note: Lenses vertexes not collapsing into each other is
            %   prescribed by ga linear constraints by Ax<=b outside of the
            %   scope of this class.

            %Inputs             
            R1=x(1); % First lens first surface curvature
            R2=x(2); % First lens second surface curvature
            R3=x(3); % Second lens first surface curvature
            R4=x(4); % Second lens second surface curvature
            d1=x(5); % Thickness of first lens on optic axis
            d=x(6); % Thickness of second lens on optic axis
            d2=x(7);  % Distance between lenses inbetween on optic axis

            XC1=-d1-d/2+R1; % Center of first curvature
            XC2=-d/2+R2; % Center of second curvature
            XC3=d/2+R3; % Center of third curvature
            XC4=d/2+d2+R4; % Center of fourth curvature

            % Intersection of horizontal line on aperture height H1 with
            % curvatures
            if R1<0
                xh11 = max(linecirc(0,H1,XC1,0,abs(R1)));
            else
                xh11 = min(linecirc(0,H1,XC1,0,abs(R1)));
            end
            if R2<0
                xh12 = max(linecirc(0,H1,XC2,0,abs(R2)));
            else
                xh12 = min(linecirc(0,H1,XC2,0,abs(R2)));
            end
            if R3<0
                xh13 = max(linecirc(0,H1,XC3,0,abs(R3)));
            else
                xh13 = min(linecirc(0,H1,XC3,0,abs(R3)));
            end
            if R4<0
                xh14 = max(linecirc(0,H1,XC4,0,abs(R4)));
            else
                xh14 = min(linecirc(0,H1,XC4,0,abs(R4)));
            end

            %  According to most workshops possibilities and glass 
            %  technology:
            %  Minimum edge thickness of the lens 2 mm for both of lenses,
            %  minimum distance between edges is 1 mm
            %  C<0
            C=[2-(xh12-xh11); 2-(xh14-xh13);1-(xh13-xh12)];
           
            % Lens maker equation for thick lens, Focal length of first lens  
            f11=1/((n2/n1-1)*(1/R1-1/R2)+(n2/n1-1)^2*d1/(n2*R1*R2)); 
            
            % Lens maker equation for thick lens, Focal length of second lens
            f22=1/((n3/n1-1)*(1/R3-1/R4)+(n3/n1-1)^2*d2/(n3*R3*R4));

            % Two lenses system eqivalent focal length should be equal to requested
            % effective focal length of the system
            % Ceq=0
            % Ceq=[abs(EFL-f11*f22/(f11+f22-x(7)))];
            Ceq=[abs(-1/EFL ...
                +1/f11+1/f22-d/(f11*f22*n1)...
                -(n2-n1)*d1/(n2*R1*f22)...
                +(n3-n1)*d2/(n3*R4*f11))];
        end
    end
end

