% Taken from Dr. Bilal Siddiqui's code, modified to 
%   be a useful utility function
% https://www.mathworks.com/matlabcentral/fileexchange/49941-3d-mohr-s-circle

% Some modifications made by me (Lucas), and are generally noted

function [tau_max, sigma_1, sigma_2, sigma_3] = mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
    % 3D-Stresses using Mohr's Method
    % by Dr. Bilal Siddiqui
    % Assistant Professor, Mechanical Engineering
    % DHA Suffa University
    %Inputs
    % sigma_x=stress_mat(1, 1);
    % sigma_y=stress_mat(2, 2);
    % sigma_z=stress_mat(3, 3);
    % tau_xy=stress_mat(2, 1);
    % tau_yz=stress_mat(3,2);
    % tau_xz=stress_mat(3,1);

    % TODO -- check! -- LJ
    sigma_x=-sigma_x
    sigma_y=-sigma_y
    sigma_z=-sigma_z
    tau_xy=tau_xy
    tau_yz=tau_yz
    tau_xz=tau_xz

    %coefficients of Mohr's circle
    a=1;
    b=sigma_x+sigma_y+sigma_z;
    c=sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z...
        -tau_xy^2-tau_yz^2-tau_xz^2;
    d=sigma_x*sigma_y*sigma_z+2*tau_xy*tau_yz*tau_xz...
        -sigma_x*tau_yz^2-sigma_y*tau_xz^2-sigma_z*tau_xy^2;
    %Find Principle stresses
    normal_stresses=roots([a b c d]);
    %arrange in descending order
    ordered_normal_stresses=sort(normal_stresses,'descend');
    sigma_1=ordered_normal_stresses(1);
    sigma_2=ordered_normal_stresses(2);
    sigma_3=ordered_normal_stresses(3);
    tau_1  =(sigma_1-sigma_3)/2;
    tau_2  =(sigma_1-sigma_2)/2;
    tau_3  =(sigma_2-sigma_3)/2;

    % Added by Lucas -- LJ
    fprintf("sigma_1 = %.5f\n", sigma_1)
    fprintf("sigma_2 = %.5f\n", sigma_2)
    fprintf("sigma_3 = %.5f\n", sigma_3)
    fprintf("tau_1 = %.5f\n", tau_1)
    fprintf("tau_2 = %.5f\n", tau_2)
    fprintf("tau_3 = %.5f\n", tau_3)


    %Plotting the 3-D Mohr's Circle
    angles=0:0.01:2*pi;
    center1=[(sigma_1+sigma_3)/2 0];
    center2=[(sigma_1+sigma_2)/2 0];
    center3=[(sigma_2+sigma_3)/2 0];
    cirlce1=[center1(1)+tau_1*cos(angles') center1(2)+tau_1*sin(angles')];
    cirlce2=[center2(1)+tau_2*cos(angles') center2(2)+tau_2*sin(angles')];
    cirlce3=[center3(1)+tau_3*cos(angles') center3(2)+tau_3*sin(angles')];
    plot(cirlce1(:,1),cirlce1(:,2),'b',cirlce2(:,1),cirlce2(:,2),'g',...
        cirlce3(:,1),cirlce3(:,2),'r');axis equal;grid on;
    %ANNOTATIONS
    if sigma_1>0
        text(sigma_1*1.01,0,'\sigma_1','fontsize',15);
    else
        text(sigma_1*0.95,0,'\sigma_1','fontsize',15)
    end
    if sigma_2>0
        text(sigma_2*1.1,0,'\sigma_2','fontsize',15)
    else
        text(sigma_2*0.99,0,'\sigma_2','fontsize',15)
    end
    if sigma_3>0
        text(sigma_3*1.1,0,'\sigma_3','fontsize',15);
    else
        text(sigma_3*0.99,0,'\sigma_3','fontsize',15);
    end
    text(center1(1),tau_1*0.9,'\tau_1','fontsize',15)
    text(center2(1),tau_2*0.9,'\tau_2','fontsize',15)
    text(center3(1),tau_3*0.9,'\tau_3','fontsize',15);
    xlabel('Normal Stress, \sigma','fontsize',15);
    ylabel('Shear Stess, \tau','fontsize',15)

    % Another addition by me -- LJ
    tau_max = tau_1
end
