classdef Mohrdapp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        UIAxes2                       matlab.ui.control.UIAxes
        UIAxes                        matlab.ui.control.UIAxes
        RotationControlPanel          matlab.ui.container.Panel
        RotZKnobLabel                 matlab.ui.control.Label
        RotZKnob                      matlab.ui.control.Knob
        RotYKnobLabel                 matlab.ui.control.Label
        RotYKnob                      matlab.ui.control.Knob
        RotXKnobLabel                 matlab.ui.control.Label
        RotXKnob                      matlab.ui.control.Knob
        UITable                       matlab.ui.control.Table
        MaximumPrincipleStressEditFieldLabel  matlab.ui.control.Label
        MaximumPrincipleStressEditField  matlab.ui.control.EditField
        AverageStressEditFieldLabel   matlab.ui.control.Label
        AverageStressEditField        matlab.ui.control.EditField
        MinimumPrincipleStressEditFieldLabel  matlab.ui.control.Label
        MinimumPrincipleStressEditField  matlab.ui.control.EditField
        MaximumShearStressEditFieldLabel  matlab.ui.control.Label
        MaximumShearStressEditField   matlab.ui.control.EditField
        EditField7_5                  matlab.ui.control.EditField
        EditField7_6                  matlab.ui.control.EditField
        EditField7_7                  matlab.ui.control.EditField
        EditField7_8                  matlab.ui.control.EditField
        DrawButton                    matlab.ui.control.Button
        ResetButton                   matlab.ui.control.Button
        InputStressTensorMPaPanel     matlab.ui.container.Panel
        EditField_8                   matlab.ui.control.NumericEditField
        XLabel                        matlab.ui.control.Label
        YLabel                        matlab.ui.control.Label
        ZLabel                        matlab.ui.control.Label
        Label_5                       matlab.ui.control.Label
        Label_6                       matlab.ui.control.Label
        EditField_9                   matlab.ui.control.NumericEditField
        EditField_10                  matlab.ui.control.NumericEditField
        EditField_11                  matlab.ui.control.NumericEditField
        EditField_12                  matlab.ui.control.NumericEditField
        EditField_13                  matlab.ui.control.NumericEditField
        XYLabel                       matlab.ui.control.Label
        XZLabel                       matlab.ui.control.Label
        YZLabel                       matlab.ui.control.Label
        CurrentStressTensorMPaPanel   matlab.ui.container.Panel
        EditField_14                  matlab.ui.control.NumericEditField
        XLabel_2                      matlab.ui.control.Label
        YLabel_2                      matlab.ui.control.Label
        ZLabel_2                      matlab.ui.control.Label
        Label_8                       matlab.ui.control.Label
        EditField_15                  matlab.ui.control.NumericEditField
        EditField_16                  matlab.ui.control.NumericEditField
        Label_9                       matlab.ui.control.Label
        EditField_17                  matlab.ui.control.NumericEditField
        EditField_18                  matlab.ui.control.NumericEditField
        EditField_19                  matlab.ui.control.NumericEditField
        XYLabel_2                     matlab.ui.control.Label
        XZLabel_2                     matlab.ui.control.Label
        YZLabel_2                     matlab.ui.control.Label
        RotationAnglesPanel           matlab.ui.container.Panel
        EditField_22                  matlab.ui.control.NumericEditField
        EditField_23                  matlab.ui.control.NumericEditField
        EditField_24                  matlab.ui.control.NumericEditField
        RotationaboutXaxisLabel       matlab.ui.control.Label
        RotationaboutYaxisLabel       matlab.ui.control.Label
        RotationaboutZaxisLabel       matlab.ui.control.Label
        AbsoluteStressEditFieldLabel  matlab.ui.control.Label
        AbsoluteStressEditField       matlab.ui.control.EditField
        EditField7_9                  matlab.ui.control.EditField
        Label_10                      matlab.ui.control.Label
    end

    
    methods (Access = private)
        
         function rotation = rotate(app,matrix,rotx,roty,rotz)
                    
                %%about x-axis
       
        w=[1      0             0           0;
           0      cosd(rotx)   -sind(rotx)      0; 
           0      sind(rotx)   cosd(rotx)       0;
           0      0              0          1];
    
                %%about y-axis
       
        r=[cosd(roty)   ,  0,    sind(roty)      0; 
            0            1 ,      0          0;
            -sind(roty)    0     cosd(roty)      0;
            0            0        0          1];    
                %%about z-axis
  
        t=[cosd(rotz)  -sind(rotz)  0  0;
           sind(rotz)  cosd(rotz)   0  0;
            0       0           1  0;
            0       0           0  1];
       
            rotation =w*r*t*matrix; 
        end
        
        function cube_3D(app,rotx,roty,rotz)
            
         cla(app.UIAxes,'reset');
     
        sigma_x = -(app.EditField_8.Value);
        sigma_y = -(app.EditField_9.Value);
        sigma_z = -(app.EditField_13.Value);
         
        tau_xy = -(app.EditField_10.Value);
        tau_yz = -(app.EditField_11.Value);
        tau_zx = -(app.EditField_12.Value);
       
        sigma_x_dash=-((sigma_x+sigma_y)/2 + ((sigma_x-sigma_y)/2)*cos(2*rotx)+tau_xy*sin(2*rotx)); 
        sigma_y_dash=-((sigma_y+sigma_z)/2 + ((sigma_y-sigma_z)/2)*cos(2*roty)+tau_yz*sin(2*roty));
        sigma_z_dash=-((sigma_z+sigma_x)/2 + ((sigma_z-sigma_x)/2)*cos(2*rotz)+tau_zx*sin(2*rotz));
       
        tau_xy_dash=-((-(sigma_x-sigma_y)/2)*sin(2*rotx)+tau_xy*cos(2*rotx));
        tau_yz_dash=-((-(sigma_y-sigma_z)/2)*sin(2*roty)+tau_yz*cos(2*roty));
        tau_zx_dash=-((-(sigma_z-sigma_x)/2)*sin(2*rotz)+tau_zx*cos(2*rotz));

        tau_yx_dash=tau_xy_dash;
        tau_zy_dash=tau_yz_dash;
        tau_xz_dash=tau_zx_dash;   

        
         view(app.UIAxes,3);

        c=[-.5 -.5 -.5 ;
           -.5 .5 -.5 ;
           .5 .5 -.5 ;
           .5 -.5 -.5 ;
           -.5 -.5 .5 ;
           -.5 .5 .5 ;
           .5 .5 .5 ;
           .5 -.5 .5 ];
        c=c';
        a=[c; ones(1,8)];
        
        faces = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
        
        r=[0   0    0   0    0  -0.4;
           0   0    0  -0.4  0   0  ;
           0.5 0.4  0.5 0    0.5 0  ];
        mt1=[r; 1 1 1 1 1 1];
        rt1 = rotate(app,mt1,rotx,roty,rotz);
        
        q=[ 0    0     0   0     0   -0.4;
           -0.5 -0.4  -0.5 0    -0.5  0  ;
            0    0     0   0.4   0    0  ];
        mt2=[q ; 1 1 1 1 1 1];
        rt2 = rotate(app,mt2,rotx,roty,rotz);
        
        k=[-0.5 -0.4 -0.5  0   -0.5  0  ;
            0    0    0    0    0   -0.4;
            0    0    0    0.4  0    0  ];
        mt3=[k ; 1 1 1 1 1 1];
        rt3 = rotate(app,mt3,rotx,roty,rotz);
      
        axis(app.UIAxes, 'equal');
        set(app.UIAxes, 'View', [45 45]);%axis equal and look good

        quiver3(app.UIAxes,rt1(1,1),rt1(2,1),rt1(3,1),rt1(1,2),rt1(2,2),rt1(3,2),'linewidth',2,'color','g');
        hold (app.UIAxes,"on");
        quiver3(app.UIAxes,rt1(1,3),rt1(2,3),rt1(3,3),rt1(1,4),rt1(2,4),rt1(3,4),'linewidth',2,'color','b');
        quiver3(app.UIAxes,rt1(1,5),rt1(2,5),rt1(3,5),rt1(1,6),rt1(2,6),rt1(3,6),'linewidth',2,'color','y');
        
        quiver3(app.UIAxes,rt2(1,1),rt2(2,1),rt2(3,1),rt2(1,2),rt2(2,2),rt2(3,2),'linewidth',2,'color','b');
        quiver3(app.UIAxes,rt2(1,3),rt2(2,3),rt2(3,3),rt2(1,4),rt2(2,4),rt2(3,4),'linewidth',2,'color','g');
        quiver3(app.UIAxes,rt2(1,5),rt2(2,5),rt2(3,5),rt2(1,6),rt2(2,6),rt2(3,6),'linewidth',2,'color','y');
        
        quiver3(app.UIAxes,rt3(1,1),rt3(2,1),rt3(3,1),rt3(1,2),rt3(2,2),rt3(3,2),'linewidth',2,'color','y');
        quiver3(app.UIAxes,rt3(1,3),rt3(2,3),rt3(3,3),rt3(1,4),rt3(2,4),rt3(3,4),'linewidth',2,'color','g');
        quiver3(app.UIAxes,rt3(1,5),rt3(2,5),rt3(3,5),rt3(1,6),rt3(2,6),rt3(3,6),'linewidth',2,'color','b');
        
        quiver3(app.UIAxes,1,1,-1,-2,0,0,'linewidth',2,'color','y');
        quiver3(app.UIAxes,1,1,-1,0,-2,0,'linewidth',2,'color','b');
        quiver3(app.UIAxes,1,1,-1,0,0,2,'linewidth',2,'color','g');
        %patch('Vertices',a,'Faces',faces,'FaceColor','r','facealpha',0.3);
        
        h = rotate(app,a,rotx,roty,rotz);
        h=h';
        results=h(:,(1:3));
        hold (app.UIAxes,"on");
        patch(app.UIAxes,'Vertices',results,'Faces',faces,'FaceColor','m','facealpha',1);
        
        p1=[-0.8 0.1 0;
            0.1 -0.8 -0.1;
            0 0 0.8];
        p2=[p1;1 1 1];
        p3=rotate(app,p2,rotx,roty,rotz);
        
        text(app.UIAxes,p3(1,1),p3(2,1),p3(3,1),sprintf('%0.2f',sigma_x_dash),'color','k','fontsize',11);
        text(app.UIAxes,p3(1,2),p3(2,2),p3(3,2),sprintf('%0.2f',sigma_y_dash),'color','k','fontsize',11);
        text(app.UIAxes,p3(1,3),p3(2,3),p3(3,3),sprintf('%0.2f',sigma_z_dash),'color','k','fontsize',11);
        
        p4=[-0.5  -0.1 -0.35;
            -0.35 -0.5 -0.1;
             0.1   0.4  0.5];
        p5=[p4;1 1 1];
        p6=rotate(app,p5,rotx,roty,rotz);
      
        text(app.UIAxes,p6(1,1),p6(2,1),p6(3,1),sprintf('%0.2f',tau_xy_dash),'color','k','fontsize',11);
        text(app.UIAxes,p6(1,2),p6(2,2),p6(3,2),sprintf('%0.2f',tau_yz_dash),'color','k','fontsize',11);
        text(app.UIAxes,p6(1,3),p6(2,3),p6(3,3),sprintf('%0.2f',tau_zx_dash),'color','k','fontsize',11);
        
        p7=[-0.4  0.1  -0.5;
            -0.5 -0.35  0.1;
            -0.1  0.5   0.4];
        p8=[p7;1 1 1];
        p9=rotate(app,p8,rotx,roty,rotz);
   
        text(app.UIAxes,p9(1,1),p9(2,1),p9(3,1),sprintf('%0.2f',tau_yx_dash),'color','k','fontsize',11);
        text(app.UIAxes,p9(1,2),p9(2,2),p9(3,2),sprintf('%0.2f',tau_zy_dash),'color','k','fontsize',11);
        text(app.UIAxes,p9(1,3),p9(2,3),p9(3,3),sprintf('%0.2f',tau_xz_dash),'color','k','fontsize',11);
        
        xlim(app.UIAxes,[-1 1]);
        ylim(app.UIAxes,[-1 1]);
        zlim(app.UIAxes,[-1 1]);
         hold (app.UIAxes,"on");

        grid (app.UIAxes,"on");
         
        end
        
        
        function  Tensor(app,rotx,roty,rotz)
       
        sigma_x = -(app.EditField_8.Value);
        sigma_y = -(app.EditField_9.Value);
        sigma_z = -(app.EditField_13.Value);
         
        tau_xy = -(app.EditField_10.Value);
        tau_yz = -(app.EditField_11.Value);
        tau_zx = -(app.EditField_12.Value);
       
        sigma_x_dash=-((sigma_x+sigma_y)/2 + ((sigma_x-sigma_y)/2)*cos(2*rotx)+tau_xy*sin(2*rotx)); 
        sigma_y_dash=-((sigma_y+sigma_z)/2 + ((sigma_y-sigma_z)/2)*cos(2*roty)+tau_yz*sin(2*roty));
        sigma_z_dash=-((sigma_z+sigma_x)/2 + ((sigma_z-sigma_x)/2)*cos(2*rotz)+tau_zx*sin(2*rotz));
       
        tau_xy_dash=-((-(sigma_x-sigma_y)/2)*sin(2*rotx)+tau_xy*cos(2*rotx));
        tau_yz_dash=-((-(sigma_y-sigma_z)/2)*sin(2*roty)+tau_yz*cos(2*roty));
        tau_zx_dash=-((-(sigma_z-sigma_x)/2)*sin(2*rotz)+tau_zx*cos(2*rotz));

        tau_yx_dash=tau_xy_dash;
        tau_zy_dash=tau_yz_dash;
        tau_xz_dash=tau_zx_dash;    
       
        app.EditField_14.Value=double(sigma_x_dash);
        app.EditField_15.Value=double(sigma_y_dash);
        app.EditField_16.Value=double(sigma_z_dash);
       
        app.EditField_17.Value=double(tau_xy_dash);
        app.EditField_18.Value=double(tau_yz_dash);
        app.EditField_19.Value=double(tau_zx_dash);
        
        end
        
        function restart(app)
       
        cla(app.UIAxes,'reset');
        cla(app.UIAxes2,'reset');
    
        Tensor(app,0,0,0);
        
        app.MaximumPrincipleStressEditField.Value= num2str(0);
        app.MinimumPrincipleStressEditField.Value= num2str(0);
        app.MaximumShearStressEditField.Value=num2str(0);
        app.AbsoluteStressEditField.Value=num2str(0);
        app.AverageStressEditField.Value=num2str(0);
        
        app.EditField_8.Value=0;
        app.EditField_9.Value=0;
        app.EditField_13.Value=0;
         
        app.EditField_10.Value=0;
        app.EditField_11.Value=0;
        app.EditField_12.Value=0;
        
        app.EditField_14.Value=double(0);
        app.EditField_15.Value=double(0);
        app.EditField_16.Value=double(0);
       
        app.EditField_17.Value=double(0);
        app.EditField_18.Value=double(0);
        app.EditField_19.Value=double(0);
        
        app.EditField_23.Value=0;
        app.EditField_24.Value=0;
        app.EditField_22.Value=0;
        
        app.RotZKnob.Value=0;
        app.RotYKnob.Value=0;
        app.RotXKnob.Value=0;
        
        end
        
        
        function Mohr_3d_2d(app)
           
        cla(app.UIAxes2, "reset");
        sigma_x = -(app.EditField_8.Value);
        sigma_y = -(app.EditField_9.Value);
        sigma_z = -(app.EditField_13.Value);
         
        tau_xy = -(app.EditField_10.Value);
        tau_yz = -(app.EditField_11.Value);
        tau_zx = -(app.EditField_12.Value);
        
        tau_yx=tau_xy;
        tau_zy=tau_yz;
        tau_xz=tau_zx;
         
        % Coefficients for Mohr's Ciecle
        c3 = 1;
        c2 = sigma_x+sigma_y+sigma_z;
        c1 = sigma_x*sigma_y+sigma_y*sigma_z+sigma_z*sigma_x -tau_xy^2 -tau_yz^2 -tau_zx^2;
        c0 = sigma_x*sigma_y*sigma_z + 2*tau_xy*tau_yz*tau_zx -(sigma_x*tau_yz^2+sigma_y*tau_zx^2+sigma_z*tau_xy^2);
         
        % Principal stresses
        normal_stresses = roots([c3 c2 c1 c0]);
        A = sort(normal_stresses,'descend');
         
        sigma_1 = A(1);
        sigma_2 = A(2);
        sigma_3 = A(3);
        sigma = [sigma_1;sigma_2;sigma_3];
         
        sigma_max = max(sigma);
        sigma_min = min(sigma);
         
        tau_1 = (sigma_1-sigma_3)/2;
        tau_2 = (sigma_1-sigma_2)/2;
        tau_3 = (sigma_2-sigma_3)/2;
        tau = [tau_1;tau_2;tau_3];
         
        tau_max = max(tau);
      
        absolute_stress=abs((sigma_max-sigma_min)/2);
       
        app.MaximumPrincipleStressEditField.Value= num2str(sigma_max);
        app.MinimumPrincipleStressEditField.Value= num2str(sigma_min);
        app.MaximumShearStressEditField.Value=num2str(tau_max);
       
        app.AbsoluteStressEditField.Value=num2str(absolute_stress);
        
        % Plotting
        theta = 0:0.01:2*pi;
         
        C1 = [(sigma_1+sigma_3)/2 0];
        C2 = [(sigma_1+sigma_2)/2 0];
        C3 = [(sigma_2+sigma_3)/2 0];
        
        %cirlce1=[C1(1)+tau_1*cos(theta') C1(2)+tau_1*sin(theta')];
        %cirlce2=[C2(1)+tau_2*cos(theta') C2(2)+tau_2*sin(theta')];
        %cirlce3=[C3(1)+tau_3*cos(theta') C3(2)+tau_3*sin(theta')];
         
        if sigma_z==0 && tau_yz==0 && tau_xz==0
           % sigma_max_2D=(sigma_x+sigma_y)/2+(((sigma_x-sigma_y)/2)^(2) + (tau_xy)^(2))^(1/2)
            
          %  sigma_max_2D = max(sigma_max_2D);
          %  sigma_min_2D = min(MAt_MAX_MIN_2D);
            
         %   app.MaximumPrincipleStressEditField.Value= num2str(sigma_max_2D);
          %  app.MinimumPrincipleStressEditField.Value= num2str(sigma_min_2D);
             
%             
%              Cemter_2d=(sigma_x+sigma_y)/2;
%              raduis_2d=(((sigma_x+sigma_y)/2)^(2)+(tau_xy)^(2))^(1/2);
%              sigma_max_2D = Cemter_2d+raduis_2d;
%              sigma_min_2D = Cemter_2d-raduis_2d;
             
            plot(app.UIAxes2,C1(1)+tau_1*cos(theta),C1(2)+tau_1*sin(theta),'m');
            axis(app.UIAxes2,"equal")
            grid (app.UIAxes2,"on")
            fill(app.UIAxes2,C1(1)+tau_1*cos(theta),C1(2)+tau_1*sin(theta),"m");
            
            if sigma_1>0
               
                text(app.UIAxes2,sigma_1*1.01,0,'\sigma_1','fontsize',15);
            else
               
                text(app.UIAxes2,sigma_1*0.95,0,'\sigma_1','fontsize',15);
            end
            
            if sigma_2>0
                
                text(app.UIAxes2,-tau_1*1.0,0,'\sigma_2','fontsize',15);
            else
               
                text(app.UIAxes2,-tau_1*0.99,0,'\sigma_2','fontsize',15);
            end
            
            text(app.UIAxes2,C1(1),tau_1*0.9,'\tau_1','fontsize',15);
            
            %average for the 2d 
            average_stress_2D=(sigma_x+sigma_y)/2;
            app.AverageStressEditField.Value=num2str(average_stress_2D);
        
        elseif (sigma_x~=0 && sigma_y~=0 && sigma_z~=0) && (tau_xy~=0 || tau_yz~=0 || tau_xz~=0)
           
             sigma_1 = A(1);
             sigma_2 = A(2);
             sigma_3 = A(3);
             sigma = [sigma_1;sigma_2;sigma_3];
         
             sigma_max = max(sigma);
             sigma_min = min(sigma);
             app.MaximumPrincipleStressEditField.Value= num2str(sigma_max);
             app.MinimumPrincipleStressEditField.Value= num2str(sigma_min);
          
            plot(app.UIAxes2,C1(1)+tau_1*cos(theta),C1(2)+tau_1*sin(theta),'m');
            hold (app.UIAxes2,"on")
            plot(app.UIAxes2,C2(1)+tau_2*cos(theta),C2(2)+tau_2*sin(theta),'g');
            hold (app.UIAxes2,"on") 
            plot(app.UIAxes2,C3(1)+tau_3*cos(theta),C3(2)+tau_3*sin(theta),'r');
            hold (app.UIAxes2,"on")
            axis(app.UIAxes2,"equal")
            grid (app.UIAxes2,"on")
        
            if sigma_1>0
                text(app.UIAxes2,sigma_1*1.01,0,'\sigma_1','fontsize',15);
            else
                text(app.UIAxes2,sigma_1*0.95,0,'\sigma_1','fontsize',15);
            end
         
            if sigma_2>0
                text(app.UIAxes2,sigma_2*1.0,0,'\sigma_2','fontsize',15);
            else
                text(app.UIAxes2,sigma_2*0.99,0,'\sigma_2','fontsize',15);
            end
         
            if sigma_3>0
                text(app.UIAxes2,sigma_3*1.1,0,'\sigma_3','fontsize',15);
            else
                text(app.UIAxes2,sigma_3*0.99,0,'\sigma_3','fontsize',15);
            end
            
            text(app.UIAxes2,C1(1),tau_1*0.9,'\tau_1','fontsize',15);
            text(app.UIAxes2,C2(1),tau_2*0.9,'\tau_2','fontsize',15);
            text(app.UIAxes2,C3(1),tau_3*0.9,'\tau_3','fontsize',15);
            
             %avergae for the 3d 
            average_stress_3D=(sigma_x+sigma_y+sigma_z)/3;
            app.AverageStressEditField.Value=num2str(average_stress_3D);
        
        elseif (tau_xy~=0 && tau_yz~=0 && tau_xz~=0) && (sigma_x~=0 || sigma_y~=0 || sigma_z~=0)
           
            sigma_1 = A(1);
             sigma_2 = A(2);
             sigma_3 = A(3);
             sigma = [sigma_1;sigma_2;sigma_3];
         
             sigma_max = max(sigma);
             sigma_min = min(sigma);
             app.MaximumPrincipleStressEditField.Value= num2str(sigma_max);
             app.MinimumPrincipleStressEditField.Value= num2str(sigma_min);
            
            plot(app.UIAxes2,C1(1)+tau_1*cos(theta),C1(2)+tau_1*sin(theta),'m');
            hold (app.UIAxes2,"on")
            plot(app.UIAxes2,C2(1)+tau_2*cos(theta),C2(2)+tau_2*sin(theta),'g');
            hold (app.UIAxes2,"on") 
            plot(app.UIAxes2,C3(1)+tau_3*cos(theta),C3(2)+tau_3*sin(theta),'r');
            hold (app.UIAxes2,"on")
            axis(app.UIAxes2,"equal")
            grid (app.UIAxes2,"on")
             
            
             
            
            if sigma_1>0
                text(app.UIAxes2,sigma_1*1.01,0,'\sigma_1','fontsize',15);
            else
                text(app.UIAxes2,sigma_1*0.95,0,'\sigma_1','fontsize',15);
            end
         
            if sigma_2>0
                text(app.UIAxes2,sigma_2*1.0,0,'\sigma_2','fontsize',15);
            else
                text(app.UIAxes2,sigma_2*0.99,0,'\sigma_2','fontsize',15);
            end
         
            if sigma_3>0
                text(app.UIAxes2,sigma_3*1.1,0,'\sigma_3','fontsize',15);
            else
                text(app.UIAxes2,sigma_3*0.99,0,'\sigma_3','fontsize',15);
            end
            
            text(app.UIAxes2,C1(1),tau_1*0.9,'\tau_1','fontsize',15);
            text(app.UIAxes2,C2(1),tau_2*0.9,'\tau_2','fontsize',15);
            text(app.UIAxes2,C3(1),tau_3*0.9,'\tau_3','fontsize',15);
            
            %avergae for the 3d 
            average_stress_3D=(sigma_x+sigma_y+sigma_z)/3;
            app.AverageStressEditField.Value=num2str(average_stress_3D);
             
        else
           
            sigma_1 = A(1);
            sigma_2 = A(2);
            sigma_3 = A(3);
            sigma = [sigma_1;sigma_2;sigma_3];
         
            sigma_max = max(sigma);
            sigma_min = min(sigma);
            app.MaximumPrincipleStressEditField.Value= num2str(sigma_max);
            app.MinimumPrincipleStressEditField.Value= num2str(sigma_min);
            
            plot(app.UIAxes2,C1(1)+tau_1*cos(theta),C1(2)+tau_1*sin(theta),'m');
            hold (app.UIAxes2,"on")
            plot(app.UIAxes2,C2(1)+tau_2*cos(theta),C2(2)+tau_2*sin(theta),'g');
            hold (app.UIAxes2,"on") 
            plot(app.UIAxes2,C3(1)+tau_3*cos(theta),C3(2)+tau_3*sin(theta),'r');
            hold (app.UIAxes2,"on")
            axis(app.UIAxes2,"equal")
            grid (app.UIAxes2,"on")
           
         
            if sigma_1>0
                text(app.UIAxes2,sigma_1*1.01,0,'\sigma_1','fontsize',15);
            else
                text(app.UIAxes2,sigma_1*0.95,0,'\sigma_1','fontsize',15);
            end
        
            if sigma_2>0
                text(app.UIAxes2,sigma_2*1.0,0,'\sigma_2','fontsize',15);
            else
                text(app.UIAxes2,sigma_2*0.99,0,'\sigma_2','fontsize',15);
            end
        
            if sigma_3>0
                text(app.UIAxes2,sigma_3*1.1,0,'\sigma_3','fontsize',15);
            else
                text(app.UIAxes2,sigma_3*0.99,0,'\sigma_3','fontsize',15);
            end
        
            text(app.UIAxes2,C1(1),tau_1*0.9,'\tau_1','fontsize',15);
            text(app.UIAxes2,C2(1),tau_2*0.9,'\tau_2','fontsize',15);
            text(app.UIAxes2,C3(1),tau_3*0.9,'\tau_3','fontsize',15);
           
            %avergae for the 3d 
            average_stress_3D=(sigma_x+sigma_y+sigma_z)/3;
            app.AverageStressEditField.Value=num2str(average_stress_3D);
       
        end
      
        xlabel(app.UIAxes2,'Normal Stress,\sigma (MPa)','fontsize',15);
        ylabel(app.UIAxes2,'Shear Stress,\tau (MPa)','fontsize',15);
        title(app.UIAxes2,'3D Mohr Circle','fontsize',15);

        
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: DrawButton
        function DrawButtonPushed(app, event)
        cla(app.UIAxes2, "reset");
                
       Mohr_3d_2d(app);
         
              
        %plotting the stress cube 
     
        Tensor(app,0,0,0);

        
       
        
        end

        % Callback function
        function RotXKnobValueChanging(app, event)
       
       
        end

        % Value changing function: RotYKnob
        function RotYKnobValueChanging(app, event)
       changingValue = event.Value;      
        app.EditField_24.Value=changingValue;
       
       
        rotx=app.RotXKnob.Value;     
        rotz=app.RotZKnob.Value;
       
        cube_3D(app,rotx,changingValue,rotz);
        Tensor(app,rotx,changingValue,rotz);
        
       
      
        end

        % Value changing function: RotZKnob
        function RotZKnobValueChanging(app, event)
        changingValue = event.Value;      
        app.EditField_23.Value=changingValue;
     
        rotx=app.RotXKnob.Value;
        roty=app.RotYKnob.Value;
       
        cube_3D(app,rotx,roty,changingValue);
        Tensor(app,rotx,roty,changingValue);
        
          
        end

        % Value changing function: RotXKnob
        function RotXKnobValueChanging2(app, event)
       changingValue = event.Value;             
       app.EditField_22.Value=changingValue; 
            
      
      
        roty=app.RotYKnob.Value;
        rotz=app.RotZKnob.Value;
       
        cube_3D(app,changingValue,roty,rotz);
        Tensor(app,changingValue,roty,rotz);
        
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            restart(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.9412 0.9412 0.9412];
            app.UIFigure.Position = [100 100 1205 771];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, {'Mohr''s Circle'; ''})
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.PlotBoxAspectRatio = [2.0573476702509 1 1];
            app.UIAxes2.Position = [1 1 565 341];

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, {'Stress Cube'; ''})
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.DataAspectRatio = [1 1 1];
            app.UIAxes.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes.View = [45 45];
            app.UIAxes.CameraTarget = [0.5 0.5 0.5];
            app.UIAxes.Position = [628 4 578 350];

            % Create RotationControlPanel
            app.RotationControlPanel = uipanel(app.UIFigure);
            app.RotationControlPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.RotationControlPanel.Title = '                                  Rotation Control';
            app.RotationControlPanel.BackgroundColor = [0.902 0.902 0.902];
            app.RotationControlPanel.Position = [857 516 349 256];

            % Create RotZKnobLabel
            app.RotZKnobLabel = uilabel(app.RotationControlPanel);
            app.RotZKnobLabel.HorizontalAlignment = 'center';
            app.RotZKnobLabel.FontWeight = 'bold';
            app.RotZKnobLabel.FontColor = [0.4941 0.1843 0.5569];
            app.RotZKnobLabel.Position = [258 8 36 22];
            app.RotZKnobLabel.Text = {'Rot Z'; ''};

            % Create RotZKnob
            app.RotZKnob = uiknob(app.RotationControlPanel, 'continuous');
            app.RotZKnob.Limits = [0 360];
            app.RotZKnob.ValueChangingFcn = createCallbackFcn(app, @RotZKnobValueChanging, true);
            app.RotZKnob.FontColor = [0.4941 0.1843 0.5569];
            app.RotZKnob.Position = [250 63 53 53];

            % Create RotYKnobLabel
            app.RotYKnobLabel = uilabel(app.RotationControlPanel);
            app.RotYKnobLabel.HorizontalAlignment = 'center';
            app.RotYKnobLabel.FontWeight = 'bold';
            app.RotYKnobLabel.FontColor = [0.4941 0.1843 0.5569];
            app.RotYKnobLabel.Position = [42 3 37 22];
            app.RotYKnobLabel.Text = {'Rot Y'; ''};

            % Create RotYKnob
            app.RotYKnob = uiknob(app.RotationControlPanel, 'continuous');
            app.RotYKnob.Limits = [0 360];
            app.RotYKnob.MajorTicks = [-180 -135 -90 -45 0 45 90 135 180 225 270 315 360];
            app.RotYKnob.ValueChangingFcn = createCallbackFcn(app, @RotYKnobValueChanging, true);
            app.RotYKnob.FontColor = [0.4941 0.1843 0.5569];
            app.RotYKnob.Position = [33 59 56 56];

            % Create RotXKnobLabel
            app.RotXKnobLabel = uilabel(app.RotationControlPanel);
            app.RotXKnobLabel.HorizontalAlignment = 'center';
            app.RotXKnobLabel.FontWeight = 'bold';
            app.RotXKnobLabel.FontColor = [0.4941 0.1843 0.5569];
            app.RotXKnobLabel.Position = [152 96 37 22];
            app.RotXKnobLabel.Text = {'Rot X'; ''};

            % Create RotXKnob
            app.RotXKnob = uiknob(app.RotationControlPanel, 'continuous');
            app.RotXKnob.Limits = [0 360];
            app.RotXKnob.ValueChangingFcn = createCallbackFcn(app, @RotXKnobValueChanging2, true);
            app.RotXKnob.FontColor = [0.4941 0.1843 0.5569];
            app.RotXKnob.Position = [141 152 58 58];

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.ColumnName = {'Parameter'; 'Value'; 'Unit'};
            app.UITable.RowName = {};
            app.UITable.ForegroundColor = [0 0.4471 0.7412];
            app.UITable.FontSize = 14;
            app.UITable.Position = [321 545 528 227];

            % Create MaximumPrincipleStressEditFieldLabel
            app.MaximumPrincipleStressEditFieldLabel = uilabel(app.UIFigure);
            app.MaximumPrincipleStressEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumPrincipleStressEditFieldLabel.FontSize = 14;
            app.MaximumPrincipleStressEditFieldLabel.FontWeight = 'bold';
            app.MaximumPrincipleStressEditFieldLabel.FontColor = [0.0745 0.6235 1];
            app.MaximumPrincipleStressEditFieldLabel.Position = [322 706 177 22];
            app.MaximumPrincipleStressEditFieldLabel.Text = {'Maximum Principle Stress'; ''};

            % Create MaximumPrincipleStressEditField
            app.MaximumPrincipleStressEditField = uieditfield(app.UIFigure, 'text');
            app.MaximumPrincipleStressEditField.FontSize = 14;
            app.MaximumPrincipleStressEditField.Position = [537 706 100 22];

            % Create AverageStressEditFieldLabel
            app.AverageStressEditFieldLabel = uilabel(app.UIFigure);
            app.AverageStressEditFieldLabel.HorizontalAlignment = 'center';
            app.AverageStressEditFieldLabel.FontSize = 14;
            app.AverageStressEditFieldLabel.FontWeight = 'bold';
            app.AverageStressEditFieldLabel.FontColor = [0.0745 0.6235 1];
            app.AverageStressEditFieldLabel.Position = [322 558 107 22];
            app.AverageStressEditFieldLabel.Text = 'Average Stress';

            % Create AverageStressEditField
            app.AverageStressEditField = uieditfield(app.UIFigure, 'text');
            app.AverageStressEditField.Position = [536 556 102 22];

            % Create MinimumPrincipleStressEditFieldLabel
            app.MinimumPrincipleStressEditFieldLabel = uilabel(app.UIFigure);
            app.MinimumPrincipleStressEditFieldLabel.HorizontalAlignment = 'right';
            app.MinimumPrincipleStressEditFieldLabel.FontSize = 14;
            app.MinimumPrincipleStressEditFieldLabel.FontWeight = 'bold';
            app.MinimumPrincipleStressEditFieldLabel.FontColor = [0.0745 0.6235 1];
            app.MinimumPrincipleStressEditFieldLabel.Position = [322 628 177 22];
            app.MinimumPrincipleStressEditFieldLabel.Text = {'Minimum Principle Stress'; ''};

            % Create MinimumPrincipleStressEditField
            app.MinimumPrincipleStressEditField = uieditfield(app.UIFigure, 'text');
            app.MinimumPrincipleStressEditField.Position = [537 628 100 22];

            % Create MaximumShearStressEditFieldLabel
            app.MaximumShearStressEditFieldLabel = uilabel(app.UIFigure);
            app.MaximumShearStressEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumShearStressEditFieldLabel.FontSize = 14;
            app.MaximumShearStressEditFieldLabel.FontWeight = 'bold';
            app.MaximumShearStressEditFieldLabel.FontColor = [0.0745 0.6235 1];
            app.MaximumShearStressEditFieldLabel.Position = [322 667 160 22];
            app.MaximumShearStressEditFieldLabel.Text = 'Maximum Shear Stress';

            % Create MaximumShearStressEditField
            app.MaximumShearStressEditField = uieditfield(app.UIFigure, 'text');
            app.MaximumShearStressEditField.Position = [537 669 100 22];

            % Create EditField7_5
            app.EditField7_5 = uieditfield(app.UIFigure, 'text');
            app.EditField7_5.FontWeight = 'bold';
            app.EditField7_5.Position = [710 591 100 22];
            app.EditField7_5.Value = 'MPa';

            % Create EditField7_6
            app.EditField7_6 = uieditfield(app.UIFigure, 'text');
            app.EditField7_6.FontWeight = 'bold';
            app.EditField7_6.Position = [710 628 100 22];
            app.EditField7_6.Value = 'Mpa';

            % Create EditField7_7
            app.EditField7_7 = uieditfield(app.UIFigure, 'text');
            app.EditField7_7.FontWeight = 'bold';
            app.EditField7_7.Position = [710 669 100 22];
            app.EditField7_7.Value = 'MPa';

            % Create EditField7_8
            app.EditField7_8 = uieditfield(app.UIFigure, 'text');
            app.EditField7_8.FontWeight = 'bold';
            app.EditField7_8.Position = [710 706 100 22];
            app.EditField7_8.Value = ' MPa';

            % Create DrawButton
            app.DrawButton = uibutton(app.UIFigure, 'push');
            app.DrawButton.ButtonPushedFcn = createCallbackFcn(app, @DrawButtonPushed, true);
            app.DrawButton.BackgroundColor = [0.651 0.651 0.651];
            app.DrawButton.FontName = 'Arial';
            app.DrawButton.FontSize = 20;
            app.DrawButton.FontWeight = 'bold';
            app.DrawButton.FontColor = [0 0.4471 0.7412];
            app.DrawButton.Position = [448 457 118 53];
            app.DrawButton.Text = 'Draw';

            % Create ResetButton
            app.ResetButton = uibutton(app.UIFigure, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.BackgroundColor = [0.651 0.651 0.651];
            app.ResetButton.FontName = 'Arial Black';
            app.ResetButton.FontSize = 20;
            app.ResetButton.FontWeight = 'bold';
            app.ResetButton.FontColor = [0.7176 0.2745 1];
            app.ResetButton.Position = [636 455 114 56];
            app.ResetButton.Text = 'Reset';

            % Create InputStressTensorMPaPanel
            app.InputStressTensorMPaPanel = uipanel(app.UIFigure);
            app.InputStressTensorMPaPanel.ForegroundColor = [1 0 0];
            app.InputStressTensorMPaPanel.Title = '              Input Stress Tensor(MPa)';
            app.InputStressTensorMPaPanel.BackgroundColor = [0.8 0.8 0.8];
            app.InputStressTensorMPaPanel.FontSize = 16;
            app.InputStressTensorMPaPanel.Position = [1 540 314 232];

            % Create EditField_8
            app.EditField_8 = uieditfield(app.InputStressTensorMPaPanel, 'numeric');
            app.EditField_8.Position = [69 118 51 35];

            % Create XLabel
            app.XLabel = uilabel(app.InputStressTensorMPaPanel);
            app.XLabel.HorizontalAlignment = 'center';
            app.XLabel.FontSize = 15;
            app.XLabel.FontWeight = 'bold';
            app.XLabel.FontColor = [1 0 0];
            app.XLabel.Position = [76 180 25 22];
            app.XLabel.Text = 'X';

            % Create YLabel
            app.YLabel = uilabel(app.InputStressTensorMPaPanel);
            app.YLabel.HorizontalAlignment = 'center';
            app.YLabel.FontSize = 15;
            app.YLabel.FontWeight = 'bold';
            app.YLabel.FontColor = [1 0 0];
            app.YLabel.Position = [171 180 25 22];
            app.YLabel.Text = 'Y';

            % Create ZLabel
            app.ZLabel = uilabel(app.InputStressTensorMPaPanel);
            app.ZLabel.HorizontalAlignment = 'center';
            app.ZLabel.FontSize = 15;
            app.ZLabel.FontWeight = 'bold';
            app.ZLabel.FontColor = [1 0 0];
            app.ZLabel.Position = [251 180 25 22];
            app.ZLabel.Text = 'Z';

            % Create Label_5
            app.Label_5 = uilabel(app.InputStressTensorMPaPanel);
            app.Label_5.FontSize = 30;
            app.Label_5.FontWeight = 'bold';
            app.Label_5.FontColor = [1 0 0];
            app.Label_5.Position = [16 117 26 36];
            app.Label_5.Text = 'ÿ';

            % Create Label_6
            app.Label_6 = uilabel(app.InputStressTensorMPaPanel);
            app.Label_6.FontName = 'Calibri';
            app.Label_6.FontSize = 30;
            app.Label_6.FontWeight = 'bold';
            app.Label_6.FontColor = [1 0 0];
            app.Label_6.Position = [17 21 25 40];
            app.Label_6.Text = 'ÿ';

            % Create EditField_9
            app.EditField_9 = uieditfield(app.InputStressTensorMPaPanel, 'numeric');
            app.EditField_9.Position = [157 118 51 35];

            % Create EditField_10
            app.EditField_10 = uieditfield(app.InputStressTensorMPaPanel, 'numeric');
            app.EditField_10.Position = [69 24 51 35];

            % Create EditField_11
            app.EditField_11 = uieditfield(app.InputStressTensorMPaPanel, 'numeric');
            app.EditField_11.Position = [157 24 51 35];

            % Create EditField_12
            app.EditField_12 = uieditfield(app.InputStressTensorMPaPanel, 'numeric');
            app.EditField_12.Position = [237 24 51 35];

            % Create EditField_13
            app.EditField_13 = uieditfield(app.InputStressTensorMPaPanel, 'numeric');
            app.EditField_13.Position = [237 118 51 35];

            % Create XYLabel
            app.XYLabel = uilabel(app.InputStressTensorMPaPanel);
            app.XYLabel.HorizontalAlignment = 'center';
            app.XYLabel.FontSize = 15;
            app.XYLabel.FontWeight = 'bold';
            app.XYLabel.FontColor = [1 0 0];
            app.XYLabel.Position = [76.5 72 26 22];
            app.XYLabel.Text = 'XY';

            % Create XZLabel
            app.XZLabel = uilabel(app.InputStressTensorMPaPanel);
            app.XZLabel.HorizontalAlignment = 'center';
            app.XZLabel.FontSize = 15;
            app.XZLabel.FontWeight = 'bold';
            app.XZLabel.FontColor = [1 0 0];
            app.XZLabel.Position = [170 72 25 22];
            app.XZLabel.Text = 'XZ';

            % Create YZLabel
            app.YZLabel = uilabel(app.InputStressTensorMPaPanel);
            app.YZLabel.HorizontalAlignment = 'center';
            app.YZLabel.FontSize = 15;
            app.YZLabel.FontWeight = 'bold';
            app.YZLabel.FontColor = [1 0 0];
            app.YZLabel.Position = [250 72 25 22];
            app.YZLabel.Text = 'YZ';

            % Create CurrentStressTensorMPaPanel
            app.CurrentStressTensorMPaPanel = uipanel(app.UIFigure);
            app.CurrentStressTensorMPaPanel.ForegroundColor = [0 0 1];
            app.CurrentStressTensorMPaPanel.TitlePosition = 'centertop';
            app.CurrentStressTensorMPaPanel.Title = '  Current Stress Tensor(MPa)';
            app.CurrentStressTensorMPaPanel.BackgroundColor = [0.8 0.8 0.8];
            app.CurrentStressTensorMPaPanel.FontSize = 16;
            app.CurrentStressTensorMPaPanel.Position = [1 341 314 213];

            % Create EditField_14
            app.EditField_14 = uieditfield(app.CurrentStressTensorMPaPanel, 'numeric');
            app.EditField_14.Position = [64 101 51 35];

            % Create XLabel_2
            app.XLabel_2 = uilabel(app.CurrentStressTensorMPaPanel);
            app.XLabel_2.HorizontalAlignment = 'center';
            app.XLabel_2.FontSize = 15;
            app.XLabel_2.FontWeight = 'bold';
            app.XLabel_2.FontColor = [0 0 1];
            app.XLabel_2.Position = [76 156 25 22];
            app.XLabel_2.Text = 'X';

            % Create YLabel_2
            app.YLabel_2 = uilabel(app.CurrentStressTensorMPaPanel);
            app.YLabel_2.HorizontalAlignment = 'center';
            app.YLabel_2.FontSize = 15;
            app.YLabel_2.FontWeight = 'bold';
            app.YLabel_2.FontColor = [0 0 1];
            app.YLabel_2.Position = [171 156 25 22];
            app.YLabel_2.Text = 'Y';

            % Create ZLabel_2
            app.ZLabel_2 = uilabel(app.CurrentStressTensorMPaPanel);
            app.ZLabel_2.HorizontalAlignment = 'center';
            app.ZLabel_2.FontSize = 15;
            app.ZLabel_2.FontWeight = 'bold';
            app.ZLabel_2.FontColor = [0 0 1];
            app.ZLabel_2.Position = [251 156 25 22];
            app.ZLabel_2.Text = 'Z';

            % Create Label_8
            app.Label_8 = uilabel(app.CurrentStressTensorMPaPanel);
            app.Label_8.FontSize = 30;
            app.Label_8.FontWeight = 'bold';
            app.Label_8.FontColor = [0 0 1];
            app.Label_8.Position = [16 99 26 36];
            app.Label_8.Text = 'ÿ';

            % Create EditField_15
            app.EditField_15 = uieditfield(app.CurrentStressTensorMPaPanel, 'numeric');
            app.EditField_15.Position = [157 101 51 35];

            % Create EditField_16
            app.EditField_16 = uieditfield(app.CurrentStressTensorMPaPanel, 'numeric');
            app.EditField_16.Position = [237 101 51 35];

            % Create Label_9
            app.Label_9 = uilabel(app.CurrentStressTensorMPaPanel);
            app.Label_9.FontName = 'Calibri';
            app.Label_9.FontSize = 30;
            app.Label_9.FontWeight = 'bold';
            app.Label_9.FontColor = [0 0 1];
            app.Label_9.Position = [16 30 25 40];
            app.Label_9.Text = 'ÿ';

            % Create EditField_17
            app.EditField_17 = uieditfield(app.CurrentStressTensorMPaPanel, 'numeric');
            app.EditField_17.Position = [66 33 51 35];

            % Create EditField_18
            app.EditField_18 = uieditfield(app.CurrentStressTensorMPaPanel, 'numeric');
            app.EditField_18.Position = [157 33 51 35];

            % Create EditField_19
            app.EditField_19 = uieditfield(app.CurrentStressTensorMPaPanel, 'numeric');
            app.EditField_19.Position = [238 33 51 35];

            % Create XYLabel_2
            app.XYLabel_2 = uilabel(app.CurrentStressTensorMPaPanel);
            app.XYLabel_2.HorizontalAlignment = 'center';
            app.XYLabel_2.FontSize = 15;
            app.XYLabel_2.FontWeight = 'bold';
            app.XYLabel_2.FontColor = [0 0 1];
            app.XYLabel_2.Position = [77 74 26 22];
            app.XYLabel_2.Text = 'XY';

            % Create XZLabel_2
            app.XZLabel_2 = uilabel(app.CurrentStressTensorMPaPanel);
            app.XZLabel_2.HorizontalAlignment = 'center';
            app.XZLabel_2.FontSize = 15;
            app.XZLabel_2.FontWeight = 'bold';
            app.XZLabel_2.FontColor = [0 0 1];
            app.XZLabel_2.Position = [171 74 25 22];
            app.XZLabel_2.Text = 'XZ';

            % Create YZLabel_2
            app.YZLabel_2 = uilabel(app.CurrentStressTensorMPaPanel);
            app.YZLabel_2.HorizontalAlignment = 'center';
            app.YZLabel_2.FontSize = 15;
            app.YZLabel_2.FontWeight = 'bold';
            app.YZLabel_2.FontColor = [0 0 1];
            app.YZLabel_2.Position = [250 74 25 22];
            app.YZLabel_2.Text = 'YZ';

            % Create RotationAnglesPanel
            app.RotationAnglesPanel = uipanel(app.UIFigure);
            app.RotationAnglesPanel.ForegroundColor = [0.6353 0.0784 0.1843];
            app.RotationAnglesPanel.TitlePosition = 'centertop';
            app.RotationAnglesPanel.Title = '    Rotation Angles';
            app.RotationAnglesPanel.BackgroundColor = [0.902 0.902 0.902];
            app.RotationAnglesPanel.Position = [857 341 349 179];

            % Create EditField_22
            app.EditField_22 = uieditfield(app.RotationAnglesPanel, 'numeric');
            app.EditField_22.FontSize = 16;
            app.EditField_22.Position = [239 110 94 35];

            % Create EditField_23
            app.EditField_23 = uieditfield(app.RotationAnglesPanel, 'numeric');
            app.EditField_23.FontSize = 16;
            app.EditField_23.Position = [239 12 95 35];

            % Create EditField_24
            app.EditField_24 = uieditfield(app.RotationAnglesPanel, 'numeric');
            app.EditField_24.FontSize = 16;
            app.EditField_24.Position = [239 59 95 35];

            % Create RotationaboutXaxisLabel
            app.RotationaboutXaxisLabel = uilabel(app.RotationAnglesPanel);
            app.RotationaboutXaxisLabel.BackgroundColor = [0 1 1];
            app.RotationaboutXaxisLabel.FontSize = 16;
            app.RotationaboutXaxisLabel.Position = [32 116 158 22];
            app.RotationaboutXaxisLabel.Text = 'Rotation about X-axis';

            % Create RotationaboutYaxisLabel
            app.RotationaboutYaxisLabel = uilabel(app.RotationAnglesPanel);
            app.RotationaboutYaxisLabel.BackgroundColor = [0.0588 1 1];
            app.RotationaboutYaxisLabel.FontSize = 16;
            app.RotationaboutYaxisLabel.Position = [32 65 157 22];
            app.RotationaboutYaxisLabel.Text = 'Rotation about Y-axis';

            % Create RotationaboutZaxisLabel
            app.RotationaboutZaxisLabel = uilabel(app.RotationAnglesPanel);
            app.RotationaboutZaxisLabel.BackgroundColor = [0.0588 1 1];
            app.RotationaboutZaxisLabel.FontSize = 16;
            app.RotationaboutZaxisLabel.Position = [32 18 158 22];
            app.RotationaboutZaxisLabel.Text = 'Rotation about Z-axis';

            % Create AbsoluteStressEditFieldLabel
            app.AbsoluteStressEditFieldLabel = uilabel(app.UIFigure);
            app.AbsoluteStressEditFieldLabel.HorizontalAlignment = 'right';
            app.AbsoluteStressEditFieldLabel.FontSize = 14;
            app.AbsoluteStressEditFieldLabel.FontWeight = 'bold';
            app.AbsoluteStressEditFieldLabel.FontColor = [0.0745 0.6235 1];
            app.AbsoluteStressEditFieldLabel.Position = [321 595 112 22];
            app.AbsoluteStressEditFieldLabel.Text = 'Absolute Stress';

            % Create AbsoluteStressEditField
            app.AbsoluteStressEditField = uieditfield(app.UIFigure, 'text');
            app.AbsoluteStressEditField.Position = [537 591 100 22];

            % Create EditField7_9
            app.EditField7_9 = uieditfield(app.UIFigure, 'text');
            app.EditField7_9.FontWeight = 'bold';
            app.EditField7_9.Position = [710 553 100 22];
            app.EditField7_9.Value = 'MPa';

            % Create Label_10
            app.Label_10 = uilabel(app.UIFigure);
            app.Label_10.BackgroundColor = [0.902 0.902 0.902];
            app.Label_10.VerticalAlignment = 'bottom';
            app.Label_10.FontName = 'Arial';
            app.Label_10.FontSize = 14;
            app.Label_10.Position = [354 366 466 48];
            app.Label_10.Text = {'-In case of using Mohr''s circle in 2D you must use it in one plane XY-plane'; '-Taking 6 values of user , while getting 6 current values on stress while '; '  rotating th Knobs on right side.   '};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mohrdapp_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end