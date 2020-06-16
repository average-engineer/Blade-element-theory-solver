%% ARAD 10 + MH 32 8.7 Prop
clear all
close all
clc
%Re = 1lakh        
%extracting aerodynamic file from excel file
    section = xlsread('QProp results Final',13,'A8:A28');
    pitch = xlsread('QProp results Final',13,'H8:H28');%pitch in radians
    chord_length = xlsread('QProp results Final',13,'B8:B28');
    lift_coeff1 = xlsread('airfoil data Re 1 lakh',3,'B24:B84');
    drag_coeff1 = xlsread('airfoil data Re 1 lakh',3,'C24:C84');
    lift_coeff2 = xlsread('airfoil data Re 1 lakh',9,'B24:B84');
    drag_coeff2 = xlsread('airfoil data Re 1 lakh',9,'C24:C84');
    %non truncated data taken
    alpha1 = xlsread('airfoil data Re 1 lakh',3,'D24:D84');%AoA values in radian
    %non truncated data taken
    alpha2 = xlsread('airfoil data Re 1 lakh',9,'D24:D84');%AoA values in radian
    %avg_chord = xlsread(FileName,1,'E20');  
    %diameter of propeller
    R = 0.1397;%[m]
    di = 2*R;
    %hub radius
    r_h = 0.006;%[m]
    %number of blades
    B = str2num(input('Enter the number of blades ','s')); 
    %RPM is 5000 rev per minute
    RPM = 5000;
    omega = (RPM*2*pi)/60;%[rad/s]
    rho = 1.225;
    %visc_k = 1.46*(10^(-5));
    %v_total = Re*visc_k/avg_chord;
    %vi_max = sqrt((v_total).^(2) - (((di/2)*0.85)*RPM)^(2));
    vi_max = str2num(input('What is the maximum velocity of the drone? ','s'))
    vi = [1:1:vi_max];
%starting the iterative process of calculation of thrust and torque
%running a loop over each velocity value of the aircraft (V infinity)
for itv=1:length(vi)
    
    %thrust
    T = 0;
    %torque
    Q = 0;
    
    vel = vi(itv);
    
    %running a loop over each section of the blade
    %for a aircraft velocity, each section will have a different velocity
    for i = 1:length(section)
        %defining the pitch angle (beta) in radians
        %we are going to use the experimentally calculated values of the
        %pitch angles 
        pitch_beta(i) = pitch(i) %+20*pi/180; %adding 15.6 degrees to the pitch angle 
        
        %assuming the values of the axial and radial inflow factors
        %axial inflow factor
        a = 1.0;
        %radial inflow factor
        ao = 1.0;
        
        %variable for checking if the error in the values of the inflow
        %factors are below the prescribed tolerance of 10^-5
        %testvalue will be equal to 0 as long as the error tolerance is
        %more than 10^-5
        %when error tolerance is less than 10^-5 then testvalue variable is
        %1
        testvalue = 0;
        
        %variable for limiting the maximum number of iterations
        check = 1;
        
        %while loop used to iterate the value of a and ao
        %we use the interpolation function to relate/curve fit the
        %experimental values of AoA and CL and CD and using that curve we
        %find the new values of alpha and cl and cd
        %we will use a spline to describe the relationship between cl,cd
        %and AoA
        
        while(~testvalue) %~testvalue means that testvalue is not equal to 1
            %radial velocity of each section
            vrad(i) = (1-ao)*omega*section(i);
            %axial velocity of each section
            vax(i) = (1+a)*vel;
            
            %resultant velocity aka relative wind velocity at each section
            v(i) = sqrt(((vrad(i))^2) + ((vax(i))^2));
            
            %the flow angle (pitch angle) for each section
            phi(i) = atan(vax(i)/vrad(i));
            
            % angle of attack = pitch angle - flow angle
            aoa(i) = pitch_beta(i)-phi(i); %the angle of attack obtained is the new angle of attack
            %we need to find the Cl and Cd corresponding to this angle of
            %attack
            
            %interploation to get the spline function defining relation between
            %experimental alphas and Cl and Cd
            %extrapval is the upper bound 
            %coefficient of lift for each section corresponding to the AoA
            %of each section
            %if the AoA calculated is higher than the largest value of the
            %experimental AoA then lift coefficient is set to a specified
            %maximum value
            if i<=18 %ARAD 10
                extrapval = 1.5641; %maximum value of lift coefficient
                Clint(i) = interp1(alpha1,lift_coeff1,aoa(i),'spline',extrapval);
            
                %similarly for drag coefficient
                extrapval = 0.05467;
                Cdint(i) = interp1(alpha1,drag_coeff1,aoa(i),'spline',extrapval);
                
            %if the new calculated AoA has a value less than the lowest
            %experimental AoA then the coefficients are set as the values
            %corresponding to the lowest experimental AoA
                   if aoa(i)<alpha1(1)
                    Cdint(i) = drag_coeff1(1);
                    Clint(i) = lift_coeff1(1);
                   end
            
            else %MH 32 8.7%
                extrapval = 1.1497; %maximum value
                Clint(i) = interp1(alpha2,lift_coeff2,aoa(i),'spline',extrapval);
            
                %similarly for drag coefficient
                extrapval = 0.06909;
                Cdint(i) = interp1(alpha2,drag_coeff2,aoa(i),'spline',extrapval);
                
            %if the new calculated AoA has a value less than the lowest
            %experimental AoA then the coefficients are set as the values
            %corresponding to the lowest experimental AoA
                
                   if aoa(i)<alpha2(1)
                    Cdint(i) = drag_coeff2(1);
                    Clint(i) = lift_coeff2(1);
                   end
            end   
            %elemental thrust for each section for all B blades
            dT(i) = (0.5*rho*(v(i)*v(i)))*chord_length(i)*B*(Clint(i)*cos(phi(i)) - Cdint(i)*sin(phi(i)));
            %elemental resisting torque for each section for all B blades
            dQ(i) = (0.5*rho*(v(i)*v(i)))*chord_length(i)*section(i)*B*(Clint(i)*sin(phi(i)) + Cdint(i)*cos(phi(i)));
            %using the relation between elemental thrust and torque and the
            %inflow axial and radial parameters to find the new parameters
            %after finding the new parameters we will find the average
            %between the new and the old (assumed) values and find the
            %difference between the average and the old values.
            %the differfence is the error toelrance which needs to be less
            %than 10^-20
            %new axial inflow factor
            a_new = dT(i)/(4*pi*section(i)*rho*(vel*vel)*(1+a)); 
            %new radial inflow parameter
            ao_new = dQ(i)/(4*pi*((section(i))^3)*rho*vel*(1+a)*omega);
            %average of the new and old values
            %loop calcualtion done
            %check if next while loop needs to be run 
            %axial inflow factor
            a_avg = (a_new+a)/2;
            %radial inflow factor
            ao_avg = (ao_new+ao)/2;
            %the values of the inflow factors are to be converged
            %error tolerance if less than 10^-20 then loop has to stop
            if abs(a_avg - a)<10^(-20) && abs(ao_avg-ao)<10^(-20)
                testvalue = 1; %while loop will stop
            end
            
            %in any case to reduce the computational time we can stop the
            %while loop after 500 iterations
            
            check = check+1;
            if check>500
                testvalue = 1;
            end
            
            %the average inflow factors now become the old factor values in
            %the next while loop
            a = a_avg;
            ao = ao_avg;
        end %while loop ended
        
        %the values of a and ao are converged numerically in the above
        %while loop and using the converged values of the factors we have
        %calculated the elemental torque and thrust
        
        %for the total thust and torque over the blade, we have to
        %integrate
        
        %elemental step for integration = radius of the propeller/number of
        %sections
        dr = (R-r_h)/length(section);
        
        %total thrust over the blades for a particular itv velocity
        T = T + (dT(i)*dr);
        %0.005 is the elementals step for a small propeller
        
        %total resisting torque over the blades for itv velocity
        Q = Q + (dQ(i)*dr);
        
        
    end %for loop for each section ended
    
    %for calculating the torque and thurst coefficient of the whole
    %propeller, we need the total torque and thrust over the propeller
    
    %the coefficients calculated will be for each velocity itv
    
    %thrust for each velocity
    thrust(itv) = T;
    
    %resisting torque for each velocity
    torque(itv) = Q;
    
    %converting RPM to rps
    rps = RPM/60;
    
    %torque coefficient for each velocity
    CQ(itv) = torque(itv)/(rho*(rps*rps)*(di^5));
    
    %thrust coefficient for each velocity
    CT(itv) = thrust(itv)/(rho*(rps*rps)*(di^4));
    
    %power coefficient
    CP(itv) = 2*pi*CQ(itv);
    
   
    
    %advance ratio for each velocity
    J(itv) = vel/(rps*di);
    
    %power generated by the motor for each velocity
    power(itv) = Q*omega;
    
    %propeller efficiency according for a particualr speed
    eff_prop1(itv) = (J(itv)*CT(itv))/(2*pi*CQ(itv));
    
    %efficiency by a slightly different formula
    eff_prop2(itv) = thrust(itv)*vel/power(itv);
    
    
    %if any of CT or CP is negative then efficiency will be zero
    if(CT(itv)<0||CP(itv)<0)
        eff_prop1(itv) = 0;
        eff_prop2(itv) = 0;
    end
    
end %for loop for each velocity is ended

%plotting

figure(1)
plot(J,eff_prop1,J,eff_prop2,'-*')
xlabel('Advance Ratio')
ylabel('Efficiencies')  


figure(2)
hold on
aa = plot(J,CQ,'color','b')
bb = plot(J,CT,'color','r')
cc = plot(J,CP,'color','g')
legend([aa,bb,cc],'Torque Coeff', 'Thrust Coeff', 'Power coeff')
xlabel('Advance Ratio')
ylabel('Coefficients')

figure(3)
hold on
xx = plot(vi,torque)
yy = plot(vi,thrust)
xlabel('Aircraft Velocity')
ylabel('Torque and Thrust')
%axis([0 140 0 1200])
legend([xx,yy], 'Torque', 'Thrust')
 

%% APC 11x4.7 Prop
clear all
close all
clc
%Re = 1lakh        
%extracting aerodynamic file from excel file
    section = xlsread('APC Props',1,'G2:G19');
    pitch = xlsread('APC Props',1,'D2:D19');%pitch in radians
    chord_length = xlsread('APC Props',1,'E2:E19');
    lift_coeff1 = xlsread('airfoil data Re 1 lakh',1,'B2:B67');
    drag_coeff1 = xlsread('airfoil data Re 1 lakh',1,'C2:C67');
    lift_coeff2 = xlsread('airfoil data Re 1 lakh',2,'B2:B114');
    drag_coeff2 = xlsread('airfoil data Re 1 lakh',2,'C2:C114');
    %non truncated data taken
    alpha1 = xlsread('airfoil data Re 1 lakh',1,'D2:D67');%AoA values in radian
    %non truncated data taken
    alpha2 = xlsread('airfoil data Re 1 lakh',2,'D2:D114');%AoA values in radian
    %avg_chord = xlsread(FileName,1,'E20');  
    %diameter of propeller
    R = 0.1397;%[m]
    di = 2*R;
    %hub radius
    r_h = 0.006;%[m]
    %number of blades
    B = str2double(input('Enter the number of blades ','s')); 
    %RPM is 5000 rev per minute
    RPM = 5000;
    omega = (RPM*2*pi)/60;%[rad/s]
    rho = 1.225;
    %visc_k = 1.46*(10^(-5));
    %v_total = Re*visc_k/avg_chord;
    %vi_max = sqrt((v_total).^(2) - (((di/2)*0.85)*RPM)^(2));
    vi_max = str2double(input('What is the maximum velocity of the drone? ','s'));
    vi = [1:1:vi_max];
%starting the iterative process of calculation of thrust and torque
%running a loop over each velocity value of the aircraft (V infinity)
for itv=1:length(vi)
    
    %thrust
    T = 0;
    %torque
    Q = 0;
    
    vel = vi(itv);
    
    %running a loop over each section of the blade
    %for a aircraft velocity, each section will have a different velocity
    for i = 1:length(section)
        %defining the pitch angle (beta) in radians
        %we are going to use the experimentally calculated values of the
        %pitch angles 
        pitch_beta(i) = pitch(i); 
        
        %assuming the values of the axial and radial inflow factors
        %axial inflow factor
        a = 0.1;
        %radial inflow factor
        ao = 0.01;
        
        %variable for checking if the error in the values of the inflow
        %factors are below the prescribed tolerance of 10^-5
        %testvalue will be equal to 0 as long as the error tolerance is
        %more than 10^-5
        %when error tolerance is less than 10^-5 then testvalue variable is
        %1
        testvalue = 0;
        
        %variable for limiting the maximum number of iterations
        check = 1;
        
        %while loop used to iterate the value of a and ao
        %we use the interpolation function to relate/curve fit the
        %experimental values of AoA and CL and CD and using that curve we
        %find the new values of alpha and cl and cd
        %we will use a spline to describe the relationship between cl,cd
        %and AoA
        
        while(~testvalue) %~testvalue means that testvalue is not equal to 1
            %radial velocity of each section
            vrad(i) = (1-ao)*omega*section(i);
            %axial velocity of each section
            vax(i) = (1+a)*vel;
            
            %resultant velocity aka relative wind velocity at each section
            v(i) = sqrt(((vrad(i))^2) + ((vax(i))^2));
            
            %the flow angle (pitch angle) for each section
            phi(i) = atan(vax(i)/vrad(i));
            
            % angle of attack = pitch angle - flow angle
            aoa(i) = pitch_beta(i)-phi(i); %the angle of attack obtained is the new angle of attack
            %we need to find the Cl and Cd corresponding to this angle of
            %attack
            
            %interploation to get the spline function defining relation between
            %experimental alphas and Cl and Cd
            %extrapval is the upper bound 
            %coefficient of lift for each section corresponding to the AoA
            %of each section
            %if the AoA calculated is higher than the largest value of the
            %experimental AoA then lift coefficient is set to a specified
            %maximum value
            if i<=15 %E63
                extrapval = 1.4681; %maximum value of lift coefficient
                Clint(i) = interp1(alpha1,lift_coeff1,aoa(i),'spline',extrapval);
            
                %similarly for drag coefficient
                extrapval = 0.10888;
                Cdint(i) = interp1(alpha1,drag_coeff1,aoa(i),'spline',extrapval);
                
            %if the new calculated AoA has a value less than the lowest
            %experimental AoA then the coefficients are set as the values
            %corresponding to the lowest experimental AoA
                   if aoa(i)<alpha1(1)
                    Cdint(i) = drag_coeff1(1);
                    Clint(i) = lift_coeff1(1);
                   end
            
            else %CLARK Y
                extrapval = 1.3694; %maximum value
                Clint(i) = interp1(alpha2,lift_coeff2,aoa(i),'spline',extrapval);
            
                %similarly for drag coefficient
                extrapval = 0.22395;
                Cdint(i) = interp1(alpha2,drag_coeff2,aoa(i),'spline',extrapval);
                
            %if the new calculated AoA has a value less than the lowest
            %experimental AoA then the coefficients are set as the values
            %corresponding to the lowest experimental AoA
                
                   if aoa(i)<alpha2(1)
                    Cdint(i) = drag_coeff2(1);
                    Clint(i) = lift_coeff2(1);
                   end
            end   
            %elemental thrust for each section for all B blades
            dT(i) = (0.5*rho*(v(i)*v(i)))*chord_length(i)*B*(Clint(i)*cos(phi(i)) - Cdint(i)*sin(phi(i)));
            %elemental resisting torque for each section for all B blades
            dQ(i) = (0.5*rho*(v(i)*v(i)))*chord_length(i)*section(i)*B*(Clint(i)*sin(phi(i)) + Cdint(i)*cos(phi(i)));
            %using the relation between elemental thrust and torque and the
            %inflow axial and radial parameters to find the new parameters
            %after finding the new parameters we will find the average
            %between the new and the old (assumed) values and find the
            %difference between the average and the old values.
            %the differfence is the error toelrance which needs to be less
            %than 10^-5
            %new axial inflow factor
            a_new = dT(i)/(4*pi*section(i)*rho*(vel*vel)*(1+a)); 
            %new radial inflow parameter
            ao_new = dQ(i)/(4*pi*((section(i))^3)*rho*vel*(1+a)*omega);
            %average of the new and old values
            %loop calcualtion done
            %check if next while loop needs to be run 
            %axial inflow factor
            a_avg = (a_new+a)/2;
            %radial inflow factor
            ao_avg = (ao_new+ao)/2;
            %the values of the inflow factors are to be converged
            %error tolerance if less than 10^-7 then loop has to stop
            if abs(a_avg - a)<10^(-20) && abs(ao_avg-ao)<10^(-20)
                testvalue = 1; %while loop will stop
            end
            
            %in any case to reduce the computational time we can stop the
            %while loop after 500 iterations
            
            check = check+1;
            if check>500
                testvalue = 1;
            end
            
            %the average inflow factors now become the old factor values in
            %the next while loop
            a = a_avg;
            ao = ao_avg;
        end %while loop ended
        
        %the values of a and ao are converged numerically in the above
        %while loop and using the converged values of the factors we have
        %calculated the elemental torque and thrust
        
        %for the total thust and torque over the blade, we have to
        %integrate
        
        %elemental step for integration = radius of the propeller/number of
        %sections
        dr = (R-r_h)/length(section);
        
        %total thrust over the blades for a particular itv velocity
        T = T + (dT(i)*dr);
        %0.005 is the elementals step for a small propeller
        
        %total resisting torque over the blades for itv velocity
        Q = Q + (dQ(i)*dr);
        
        
    end %for loop for each section ended
    
    %for calculating the torque and thurst coefficient of the whole
    %propeller, we need the total torque and thrust over the propeller
    
    %the coefficients calculated will be for each velocity itv
    
    %thrust for each velocity
    thrust(itv) = T;
    
    %resisting torque for each velocity
    torque(itv) = Q;
    
    %converting RPM to rps
    rps = RPM/60;
    
    %torque coefficient for each velocity
    CQ(itv) = torque(itv)/(rho*(rps*rps)*(di^5));
    
    %thrust coefficient for each velocity
    CT(itv) = thrust(itv)/(rho*(rps*rps)*(di^4));
    
    %power coefficient
    CP(itv) = 2*pi*CQ(itv);
    
   
    
    %advance ratio for each velocity
    J(itv) = vel/(rps*di);
    
    %power generated by the motor for each velocity
    power(itv) = Q*omega;
    
    %propeller efficiency according for a particualr speed
    eff_prop1(itv) = (J(itv)*CT(itv))/(2*pi*CQ(itv));
    
    %efficiency by a slightly different formula
    eff_prop2(itv) = thrust(itv)*vel/power(itv);
    
    
    %if any of CT or CP is negative then efficiency will be zero
    if(CT(itv)<0||CP(itv)<0)
        eff_prop1(itv) = 0;
        eff_prop2(itv) = 0;
    end
    
end %for loop for each velocity is ended

%plotting

figure(1)
plot(J,eff_prop1,J,eff_prop2,'-*')
xlabel('Advance Ratio')
ylabel('Efficiencies')  


figure(2)
hold on
aa = plot(J,CQ,'color','b')
bb = plot(J,CT,'color','r')
cc = plot(J,CP,'color','g')
legend([aa,bb,cc],'Torque Coeff', 'Thrust Coeff', 'Power coeff')
xlabel('Advance Ratio')
ylabel('Coefficients')

figure(3)
hold on
xx = plot(vi,torque)
yy = plot(vi,thrust)
xlabel('Aircraft Velocity')
ylabel('Torque and Thrust')
%axis([0 140 0 1200])
legend([xx,yy], 'Torque', 'Thrust')
 
