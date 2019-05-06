# Quadrotor Simulation Environment - Gabriel Barsi Haberfeld
# University of Illinois at Urbana-Champaign, 2019
# Heavily based on original MATLAB version by William Selby:
# https://www.wilselby.com/research/arducopter/simulation-environment/


# This script implements a basic PID position control and associated plots
# and animation. The code follows sequentially and everyrhing lives inside the
# parameter struct QuadParam. It initializes with default values.


using Parameters
using BenchmarkTools
using Plots

# If the struct can be non-mutable, restart julia to redefine
@with_kw mutable struct QuadParam @deftype Float64

        init = 0;     # used in initilization
        Ts = .01;     # Sampling time (100 Hz)
        sim_time = 10; # Simulation time (seconds)
        counter = 1;   # the counter that holds the time value

        # Plotting Variables
        t_plot::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}} = 0:Ts:sim_time-Ts;       # the time values

        # Environmental Parametes
        g = 9.81;     # Gravity (m/s^2)

        # Quadrotor Physical Parameters
        m = 1.4;      # Quadrotor mass (kg)
        l = .56;     # Distance from the center of mass to the each motor (m)
        t = .02;   #Thickness of the quadrotor's arms for drawing purposes (m)
        rot_rad = .1;   #Radius of the propellor (m)
        Kd = 1.3858e-6;    # Drag torque coeffecient (kg-m^2)

        Kdx = 0.16481;    # Translational drag force coeffecient (kg/s)
        Kdy = 0.31892;    # Translational drag force coeffecient (kg/s)
        Kdz = 1.1E-6;    # Translational drag force coeffecient (kg/s)

        Jx = .05;     # Moment of inertia about X axis (kg-m^2)
        Jy = .05;     # Moment of inertia about Y axis (kg-m^2)
        Jz = .24;    # Moment of inertia about Z axis (kg-m^2)

        # Quadrotor Sensor Paramaters
        GPS_freq = (1/Ts)/1;
        X_error = .01;  #+/- m
        Y_error = .01;  #+/- m
        Z_error = .01;  #+/- m
        phi_error = .01;
        theta_error = .01;
        psi_error = .01;

        x_acc_bias = 0.16594;  # m/s^2
        x_acc_sd = 0.0093907;
        y_acc_bias = 0.31691;  # m/s^2
        y_acc_sd = 0.011045;
        z_acc_bias = -8.6759;  # m/s^2
        z_acc_sd = 0.016189;

        x_gyro_bias = 0.00053417;  # rad/s
        x_gyro_sd = 0.00066675;
        y_gyro_bias = -0.0011035;  # rad/s
        y_gyro_sd = 0.00053642;
        z_gyro_bias = 0.00020838;  # rad/s
        z_gyro_sd = 0.0004403;

        ground_truth = 1;  # Use perfect sensor measurements
        sensor_unfiltered = 0; # Use sensor errors, no filter
        sensor_kf = 0;     # Use sensor error, Kalman Filter

        # Motor Parameters
        KT = 1.3328e-5;    # Thrust force coeffecient (kg-m)
        Jp = 0.044;     # Moment of Intertia of the rotor (kg-m^2)
        max_motor_speed = 925; # motors upper limit (rad/s)
        min_motor_speed = 0; #-1*((400)^2); # motors lower limit (can't spin in reverse)

        Obar = 0;     # sum of motor speeds (O1-O2+O3-O4, N-m)
        O1 = 0;       # Front motor speed (raidans/s)
        O2 = 0;       # Right motor speed (raidans/s)
        O3 = 0;       # Rear motor speed (raidans/s)
        O4 = 0;       # Left motor speed (raidans/s)

        #Translational Positions
        X = 0;        # Initial position in X direction GF (m)
        Y = 0;        # Initial position in Y direction GF (m)
        Z = 0;        # Initial position in Z direction GF (m)
        X_BF = 0;     # Initial position in X direction BF (m)
        Y_BF = 0;     # Initial position in Y direction BF (m)
        Z_BF = 0;     # Initial position in the Z direction BF (m)

        #Translational Velocities
        X_dot = 0;    # Initial velocity in X direction GF (m/s)
        Y_dot = 0;    # Initial velocity in Y direction GF (m/s)
        Z_dot = 0;    # Initial velocity in Z direction GF (m/s)
        X_dot_BF = 0;    # Initial velocity in X direction BF (m/s)
        Y_dot_BF = 0;    # Initial velocity in Y direction BF (m/s)
        Z_dot_BF = 0;    # Initial velocity in Y direction BF (m/s)

        #Angular Positions
        phi = 0;      # Initial phi value (rotation about X GF, roll,  radians)
        theta = 0;    # Initial theta value (rotation about Y GF, pitch, radians)
        psi = 0;      # Initial psi value (rotation about Z GF, yaw, radians)

        #Angular Velocities
        p = 0;        # Initial p value (angular rate rotation about X BF, radians/s)
        q = 0;        # Initial q value (angular rate rotation about Y BF, radians/s)
        r = 0;        # Initial r value (angular rate rotation about Z BF, radians/s)

        # Desired variables
        X_des_GF = 5;         # desired value of X in Global frame
        Y_des_GF = 5;         # desired value of Y in Global frame
        Z_des_GF = -0.5;      # desired value of Z in Global frame
        X_des = 0;            # desired value of X in Body frame
        Y_des = 0;            # desired value of Y in Body frame
        Z_des = 0;            # desired value of Z in Body frame

        phi_des = 0;          # desired value of phi (radians)
        theta_des = 0;        # desired value of theta (radians)
        psi_des = 0;          # desired value of psi (radians)

        p_des = 0;
        q_des = 0;
        r_des = 0;

        # Measured variables
        X_meas = 0;
        Y_meas = 0;
        Z_meas = 0;
        phi_meas = 0;
        theta_meas = 0;
        psi_meas = 0;

        # Disturbance Variables
        Z_dis = 0;            # Disturbance in Z direction
        X_dis = 0;            # Disturbance in X direction
        Y_dis = 0;            # Ddisturbance in Y direction
        phi_dis = 0;            # Disturbance in Yaw direction
        theta_dis = 0;            # Disturbance in Pitch direction
        psi_dis = 0;            # Disturbance in Roll direction

        # Control Inputs
        U1 = 0;       # Total thrust (N)
        U2 = 0;       # Torque about X axis BF (N-m)
        U3 = 0;       # Torque about Y axis BF (N-m)
        U4 = 0;       # Torque about Z axis BF (N-m)

        # Control Limits (update values)
        U1_max = 43.5;   # KT*4*max_motor_speed^2
        U1_min = 0;      #
        U2_max = 6.25;  # KT*l*max_motor_speed^2
        U2_min = -6.25; # KT*l*max_motor_speed^2
        U3_max = 6.25;  # KT*l*max_motor_speed^2
        U3_min = -6.25; # KT*l*max_motor_speed^2
        U4_max = 2.25;  # Kd*2*max_motor_speed^2
        U4_min = -2.25; # Kd*2*max_motor_speed^2

        # PID errors

        x_error_sum = 0;
        y_error_sum = 0
        z_error_sum = 0;
        phi_error_sum = 0;
        theta_error_sum = 0;
        psi_error_sum = 0;
        p_error_sum = 0;
        q_error_sum = 0;
        r_error_sum = 0;

        # PID parameters
        X_KP = .35;          # KP value in X position control
        X_KI = .25;            # KI value in X position control
        X_KD = -.35;         # KD value in X position control
        X_KI_lim = .25;         # Error to start calculating integral term

        Y_KP = .35;          # KP value in Y position control
        Y_KI = .25;            # KI value in Y position control
        Y_KD = -.35;         # KD value in Y position control
        Y_KI_lim = .25;         # Error to start calculating integral term

        Z_KP = 5/1.7;    # KP value in altitude control
        Z_KI = 0*3;    # KI value in altitude control
        Z_KD = -10/1.980;  # KD value in altitude control
        Z_KI_lim = .25;         # Error to start calculating integral term

        phi_KP = 4.5;      # KP value in roll control 2
        phi_KI = 0;       # KI value in roll control   1
        phi_KD = 0;     # KD value in roll control  -.5
        phi_max = pi/8;   # Maximum roll angle commanded
        phi_KI_lim = 2*(2*pi/360);  # Error to start calculating integral

        theta_KP = 4.5;    # KP value in pitch control 2
        theta_KI = 0;     # KI value in pitch control 1
        theta_KD = 0;   # KD value in pitch control -.5
        theta_max = pi/8; # Maximum pitch angle commanded
        theta_KI_lim = 2*(2*pi/360);  # Error to start calculating integral

        psi_KP = 10;     # KP value in yaw control
        psi_KI = 0;     # KI value in yaw control .75
        psi_KD = 0;     # KD value in yaw control -.5
        psi_KI_lim = 8*(2*pi/360);  # Error to start calculating integral

        p_KP = 2.7;    # KP value in pitch control 2
        p_KI = 1;     # KI value in pitch control
        p_KD = -.01;   # KD value in pitch control -.5
        p_max = 500*(2*pi/360); # Maximum pitch angle commanded
        p_KI_lim = 10*(2*pi/360);  # Error to start calculating integral

        q_KP = 2.7;    # KP value in pitch control
        q_KI = 1;     # KI value in pitch control
        q_KD = -.01;   # KD value in pitch control -.5
        q_max = 500*(2*pi/360); # Maximum pitch angle commanded
        q_KI_lim = 10*(2*pi/360);  # Error to start calculating integral

        r_KP = 2.7;    # KP value in pitch control
        r_KI = 1;     # KI value in pitch control
        r_KD = -.01;   # KD value in pitch control
        r_max = 50*(2*pi/360); # Maximum pitch angle commanded
        r_KI_lim = 10*(2*pi/360);  # Error to start calculating integral


        #dynamics
        X_ddot = (-cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)*U1-Kdx*X_dot)/m;
        Y_ddot = (-cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)*U1-Kdy*Y_dot)/m;
        Z_ddot = (-cos(phi)*cos(theta)*U1-Kdz*Z_dot)/m+g;

        p_dot = (q*r*(Jy - Jz) - Jp*p*Obar + l*U2)/Jx;
        q_dot = (p*r*(Jz - Jx) + Jp*q*Obar + l*U3)/Jy;
        r_dot = (p*q*(Jx - Jy) + U4)/Jz;

        phi_dot   = p + sin(phi)*tan(theta)*q + cos(phi)*tan(theta)*r;
        theta_dot = cos(phi)*q - sin(phi)*r;
        psi_dot   = sin(phi)/cos(theta)*q + cos(phi)/cos(theta)*r;
end

# Nonlinear Dynamics
function QuadDynamicsUpdate(Quad)

        ## Update Accelerations

        Quad.X_ddot = (-(cos(Quad.phi)*sin(Quad.theta)*cos(Quad.psi)+sin(Quad.phi)*sin(Quad.psi))*Quad.U1-Quad.Kdx*Quad.X_dot)/Quad.m;
        Quad.Y_ddot = (-(cos(Quad.phi)*sin(Quad.psi)*sin(Quad.theta)-cos(Quad.psi)*sin(Quad.phi))*Quad.U1-Quad.Kdy*Quad.Y_dot)/Quad.m;
        Quad.Z_ddot = (-(cos(Quad.phi)*cos(Quad.theta))*Quad.U1-Quad.Kdz*Quad.Z_dot)/Quad.m+Quad.g;

        Quad.p_dot = (Quad.q*Quad.r*(Quad.Jy - Quad.Jz) - Quad.Jp*Quad.p*Quad.Obar + Quad.l*Quad.U2)/Quad.Jx;
        Quad.q_dot = (Quad.p*Quad.r*(Quad.Jz - Quad.Jx) + Quad.Jp*Quad.q*Quad.Obar + Quad.l*Quad.U3)/Quad.Jy;
        Quad.r_dot = (Quad.p*Quad.q*(Quad.Jx - Quad.Jy) + Quad.U4)/Quad.Jz;

        Quad.phi_dot   = Quad.p + sin(Quad.phi)*tan(Quad.theta)*Quad.q + cos(Quad.phi)*tan(Quad.theta)*Quad.r;
        Quad.theta_dot = cos(Quad.phi)*Quad.q - sin(Quad.phi)*Quad.r;
        Quad.psi_dot   = sin(Quad.phi)/cos(Quad.theta)*Quad.q + cos(Quad.phi)/cos(Quad.theta)*Quad.r;

        ## Disturbance model

        Quad.X_ddot = Quad.X_ddot + Quad.X_dis/Quad.m;
        Quad.Y_ddot = Quad.Y_ddot + Quad.Y_dis/Quad.m;
        Quad.Z_ddot = Quad.Z_ddot + Quad.Z_dis/Quad.m;
        Quad.phi_dot = Quad.phi_dot + Quad.phi_dis/Quad.Jx*Quad.Ts;
        Quad.theta_dot = Quad.theta_dot + Quad.theta_dis/Quad.Jy*Quad.Ts;
        Quad.psi_dot = Quad.psi_dot + Quad.psi_dis/Quad.Jz*Quad.Ts;

        ## Update Velocities and Positions

        # Calculating the Z velocity & position
        Quad.Z_dot = Quad.Z_ddot*Quad.Ts + Quad.Z_dot;
        Quad.Z = Quad.Z_dot*Quad.Ts + Quad.Z;

        # Calculating the X velocity & position
        Quad.X_dot = Quad.X_ddot*Quad.Ts + Quad.X_dot;
        Quad.X = Quad.X_dot*Quad.Ts + Quad.X;

        # Calculating the Y velocity & position
        Quad.Y_dot = Quad.Y_ddot*Quad.Ts + Quad.Y_dot;
        Quad.Y = Quad.Y_dot*Quad.Ts + Quad.Y;

        # Calculating p,q,r
        Quad.p = Quad.p_dot*Quad.Ts+Quad.p;
        Quad.q = Quad.q_dot*Quad.Ts+Quad.q;
        Quad.r = Quad.r_dot*Quad.Ts+Quad.r;

        # Calculating angular velocity and position
        Quad.phi = Quad.phi_dot*Quad.Ts + Quad.phi;
        Quad.theta = Quad.theta_dot*Quad.Ts+Quad.theta;
        Quad.psi = Quad.psi_dot*Quad.Ts+Quad.psi;

        return Quad
end

function QuadSensorMeasurements(Quad)

        if Quad.ground_truth == 1
            noise = 0;
        else
            noise = 1;
        end

        ## VICON
        Quad.X_meas = Quad.X + noise*(randn()*Quad.X_error);
        Quad.Y_meas = Quad.Y + noise*(randn()*Quad.Y_error);
        Quad.Z_meas = Quad.Z + noise*(randn()*Quad.Z_error);
        Quad.phi_meas = Quad.phi + noise*(randn()*Quad.phi_error);
        Quad.theta_meas = Quad.theta + noise*(randn()*Quad.theta_error);
        Quad.psi_meas = Quad.psi + noise*(randn()*Quad.psi_error);

        ## IMU Measurements

        Quad.X_ddot = Quad.X_ddot + noise*(Quad.x_acc_bias + Quad.x_acc_sd*randn());
        Quad.Y_ddot = Quad.Y_ddot + noise*(Quad.y_acc_bias + Quad.y_acc_sd*randn());
        Quad.Z_ddot = Quad.Z_ddot + noise*(Quad.z_acc_bias + Quad.z_acc_sd*randn());

        Quad.p = Quad.p + noise*(Quad.x_gyro_bias + Quad.x_gyro_sd*randn());
        Quad.q = Quad.q + noise*(Quad.y_gyro_bias + Quad.y_gyro_sd*randn());
        Quad.r = Quad.r + noise*(Quad.z_gyro_bias + Quad.z_gyro_sd*randn());

        return Quad
end

function QuadPositionPID(Quad)

        x = Quad.X_meas;
        y = Quad.Y_meas;
        z = Quad.Z_meas;
        phi = Quad.phi_meas;
        theta = Quad.theta_meas;
        psi = Quad.psi_meas;

        # Rotate Desired Position from GF to BF (Z axis rotation only maybe)
        pts = [Quad.X_des_GF Quad.Y_des_GF Quad.Z_des_GF]'
        pts = rotateGFtoBF(pts,1*phi,1*theta,psi);
        Quad.X_des = pts[1]
        Quad.Y_des = pts[2]
        Quad.Z_des = pts[3]

        # Rotate Current Position from GF to BF
        pts = [x y z]'
        pts = rotateGFtoBF(pts,phi,theta,psi);
        Quad.X_BF = pts[1]
        Quad.Y_BF = pts[2]
        Quad.Z_BF = pts[3]

        # Rotate Current Velocity from GF to BF
        pts = [Quad.X_dot Quad.Y_dot Quad.Z_dot]'
        pts = rotateGFtoBF(pts,phi,theta,psi);
        Quad.X_dot_BF = pts[1]
        Quad.Y_dot_BF = pts[2]
        Quad.Z_dot_BF = pts[3]

        # X Position PID controller
        x_error = Quad.X_des - Quad.X_BF;
        if(abs(x_error) < Quad.X_KI_lim) #only use KI for fine tracking
            Quad.x_error_sum = Quad.x_error_sum + x_error;
        end
        cp = Quad.X_KP*x_error;    #Proportional term
        ci = Quad.X_KI*Quad.Ts*Quad.x_error_sum;
        ci = min(Quad.theta_max, max(-Quad.theta_max, ci));    #Saturate ci
        cd = Quad.X_KD*Quad.X_dot_BF;                    #Derivative term
        Quad.theta_des =  - (cp + ci + cd);   #Theta and X inversely related
        Quad.theta_des = min(Quad.theta_max, max(-Quad.theta_max, Quad.theta_des));


        # Y Position PID controller
        y_error = Quad.Y_des - Quad.Y_BF;
        if(abs(y_error) < Quad.Y_KI_lim)
            Quad.y_error_sum = Quad.y_error_sum + y_error;
        end
        cp = Quad.Y_KP*y_error;    #Proportional term
        ci = Quad.Y_KI*Quad.Ts*Quad.y_error_sum;
        ci = min(Quad.phi_max, max(-Quad.phi_max, ci));    #Saturate ci
        cd = Quad.Y_KD*Quad.Y_dot_BF;                      #Derivative term
        Quad.phi_des = cp + ci + cd;
        Quad.phi_des = min(Quad.phi_max, max(-Quad.phi_max, Quad.phi_des));


        ## Z Position PID Controller/Altitude Controller
        z_error = Quad.Z_des_GF-Quad.Z_BF;
        if(abs(z_error) < Quad.Z_KI_lim)
            Quad.z_error_sum = Quad.z_error_sum + z_error;
        end
        cp = Quad.Z_KP*z_error;         #Proportional term
        ci = Quad.Z_KI*Quad.Ts*Quad.z_error_sum; #Integral term
        ci = min(Quad.U1_max, max(Quad.U1_min, ci));    #Saturate ci
        cd = Quad.Z_KD*Quad.Z_dot;                  #Derivative term
        Quad.U1 = -(cp + ci + cd)/(cos(theta)*cos(phi)) + (Quad.m * Quad.g)/(cos(theta)*cos(phi));   #Negative since Thurst and Z inversely related
        Quad.U1 = min(Quad.U1_max, max(Quad.U1_min, Quad.U1));

        return Quad
end

function QuadAttitudePID(Quad)

        phi = Quad.phi_meas;
        theta = Quad.theta_meas;
        psi = Quad.psi_meas;

        # Roll PID Controller
        phi_error = Quad.phi_des - phi;
        if(abs(phi_error) < Quad.phi_KI_lim)
            Quad.phi_error_sum = Quad.phi_error_sum + phi_error;
        end
        cp = Quad.phi_KP*phi_error;
        ci = Quad.phi_KI*Quad.Ts*Quad.phi_error_sum;
        ci = min(Quad.p_max, max(-Quad.p_max, ci));
        cd = Quad.phi_KD*Quad.p;
        Quad.p_des = cp + ci + cd;
        Quad.p_des = min(Quad.p_max, max(-Quad.p_max, Quad.p_des));

        # Pitch PID Controller
        theta_error = Quad.theta_des - theta;
        if(abs(theta_error) < Quad.theta_KI_lim)
            Quad.theta_error_sum = Quad.theta_error_sum + theta_error;
        end
        cp = Quad.theta_KP*theta_error;
        ci = Quad.theta_KI*Quad.Ts*Quad.theta_error_sum;
        ci = min(Quad.q_max, max(-Quad.q_max, ci));
        cd = Quad.theta_KD*Quad.q;
        Quad.q_des = cp + ci + cd;
        Quad.q_des = min(Quad.q_max, max(-Quad.q_max, Quad.q_des));


        # Yaw PID Controller
        psi_error = Quad.psi_des - psi;
        if(abs(psi_error) < Quad.psi_KI_lim)
            Quad.psi_error_sum = Quad.psi_error_sum + psi_error;
        end
        cp = Quad.psi_KP*psi_error;
        ci = Quad.psi_KI*Quad.Ts*Quad.psi_error_sum;
        ci = min(Quad.r_max, max(-Quad.r_max, ci));
        cd = Quad.psi_KD*Quad.r;
        Quad.r_des = cp + ci + cd;
        Quad.r_des = min(Quad.r_max, max(-Quad.r_max, Quad.r_des));

        return Quad
end

function QuadRatePID(Quad)

        p = Quad.p;
        q = Quad.q;
        r = Quad.r;

        ## Angular Rate Controller

        # Roll PID Controller
        p_error = Quad.p_des - p;
        if(abs(p_error) < Quad.p_KI_lim)
            Quad.p_error_sum = Quad.p_error_sum + p_error;
        end
        cp = Quad.p_KP*p_error;
        ci = Quad.p_KI*Quad.Ts*Quad.p_error_sum;
        ci = min(Quad.U2_max, max(Quad.U2_min, ci));
        cd = Quad.p_KD*Quad.p_dot;
        Quad.U2 = cp + ci + cd;
        Quad.U2 = min(Quad.U2_max, max(Quad.U2_min, Quad.U2));

        # Pitch PID Controller
        q_error = Quad.q_des - q;
        if(abs(q_error) < Quad.q_KI_lim)
            Quad.q_error_sum = Quad.q_error_sum + q_error;
        end
        cp = Quad.q_KP*q_error;
        ci = Quad.q_KI*Quad.Ts*Quad.q_error_sum;
        ci = min(Quad.U3_max, max(Quad.U3_min, ci));
        cd = Quad.q_KD*Quad.q_dot;
        Quad.U3 = cp + ci + cd;
        Quad.U3 = min(Quad.U3_max, max(Quad.U3_min, Quad.U3));

        # Yaw PID Controller
        r_error = Quad.r_des - r;
        if(abs(r_error) < Quad.r_KI_lim)
            Quad.r_error_sum = Quad.r_error_sum + r_error;
        end
        cp = Quad.r_KP*r_error;
        ci = Quad.r_KI*Quad.Ts*Quad.r_error_sum;
        ci = min(Quad.U4_max, max(Quad.U4_min, ci));
        cd = Quad.r_KD*Quad.r_dot;
        Quad.U4 = cp + ci + cd;
        Quad.U4 = min(Quad.U4_max, max(Quad.U4_min, Quad.U4));


        return Quad
end

function rotateGFtoBF(pts,phi,theta,psi)

        # define rotation matrix
          R_roll = [
                  1     0           0;
                  0     cos(phi)    -sin(phi);
                  0     sin(phi)    cos(phi)];
          R_pitch = [
                  cos(theta)    0       sin(theta);
                  0             1       0;
                  -sin(theta)   0       cos(theta)];
          R_yaw = [
                  cos(psi)      -sin(psi)       0;
                  sin(psi)      cos(psi)        0;
                  0             0               1];

          R = R_yaw*R_pitch*R_roll;
          pts = pts'*R;

        return pts
end

function QuadMotorSpeeds(Quad)

        # Calculate motor speeds (rad/s)^2
        w1 = Quad.U1/(4*Quad.KT) + Quad.U3/(2*Quad.KT*Quad.l) + Quad.U4/(4*Quad.Kd);
        w2 = Quad.U1/(4*Quad.KT) - Quad.U2/(2*Quad.KT*Quad.l) - Quad.U4/(4*Quad.Kd);
        w3 = Quad.U1/(4*Quad.KT) - Quad.U3/(2*Quad.KT*Quad.l) + Quad.U4/(4*Quad.Kd);
        w4 = Quad.U1/(4*Quad.KT) + Quad.U2/(2*Quad.KT*Quad.l) - Quad.U4/(4*Quad.Kd);

        # Apply realistic motor speed limits
        if w1 > Quad.max_motor_speed^2
            w1 = Quad.max_motor_speed^2;
        end
        if w1 < Quad.min_motor_speed^2
            w1 = Quad.min_motor_speed^2;
        end

        if w2 > Quad.max_motor_speed^2
            w2 = Quad.max_motor_speed^2;
        end
        if w2 < Quad.min_motor_speed^2
            w2 = Quad.min_motor_speed^2;
        end

        if w3 > Quad.max_motor_speed^2
            w3 = Quad.max_motor_speed^2;
        end
        if w3 < Quad.min_motor_speed^2
            w3 = Quad.min_motor_speed^2;
        end

        if w4 > Quad.max_motor_speed^2
            w4 = Quad.max_motor_speed^2;
        end
        if w4 < Quad.min_motor_speed^2
            w4 = Quad.min_motor_speed^2;
        end

        Quad.O1 = sqrt(w1);    # Front M
        Quad.O2 = sqrt(w2);    # Right M
        Quad.O3 = sqrt(w3);    # Rear M
        Quad.O4 = sqrt(w4);    # Left M


        ## Re-compute traditional control inputs

        Quad.U1 = Quad.KT*(Quad.O1^2 + Quad.O2^2 + Quad.O3^2 + Quad.O4^2);
        Quad.U2 = Quad.KT*Quad.l*(Quad.O4^2 - Quad.O2^2);
        Quad.U3 = Quad.KT*Quad.l*(Quad.O1^2 - Quad.O3^2);
        Quad.U4 = Quad.Kd*(Quad.O1^2 + Quad.O3^2 - Quad.O2^2 - Quad.O4^2);
        Quad.Obar = (Quad.O1 - Quad.O2 + Quad.O3 - Quad.O4);

        return Quad
end

function QuadTorqueDisturbances(Quad,t)

        Tx = 1.0*randn() + 0.5*cos(t*5)
        Ty = 1.0*randn() + 0.5*sin(t*3)
        Tz = 1.0*randn() + 0.5*sin(t*1)

        Quad.U2 += Tx
        Quad.U3 += Ty
        Quad.U4 += Tz

        return Quad
end

function QuadSim(Xd,Yd,Zd,ground_truth)

        Quad = QuadParam()
        Quad.X_des_GF = Xd
        Quad.Y_des_GF = Yd
        Quad.Z_des_GF = Zd
        Quad.ground_truth = ground_truth
        Quad.sim_time = 10
        Quad.t_plot::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}} = 0:Quad.Ts:Quad.sim_time-Quad.Ts;
        l = size(Quad.t_plot,1) + 1
        U1 = zeros(l)
        U2 = zeros(l)
        U3 = zeros(l)
        U4 = zeros(l)
        X = zeros(l)
        Y = zeros(l)
        Z = zeros(l)
        phi = zeros(l)
        theta = zeros(l)
        psi = zeros(l)

        t_plot = 0:Quad.Ts:Quad.sim_time
        c = 1

        for t = 0:Quad.Ts:Quad.sim_time
                QuadSensorMeasurements(Quad)
                QuadPositionPID(Quad)
                QuadAttitudePID(Quad)
                QuadRatePID(Quad)
                QuadTorqueDisturbances(Quad,t)
                #QuadMotorSpeeds(Quad) # this is broken
                QuadDynamicsUpdate(Quad)
                U1[c] = Quad.U1
                U2[c] = Quad.U2
                U3[c] = Quad.U3
                U4[c] = Quad.U4
                X[c] = Quad.X
                Y[c] = Quad.Y
                Z[c] = -Quad.Z # flip for plotting (NED)
                phi[c] = Quad.phi
                theta[c] = Quad.theta
                psi[c] = Quad.psi
                c = c + 1
        end

        return Quad, t_plot, U1, U2, U3, U4, X, Y, Z, phi, theta, psi
end

function QuadPlot(j, plotlim, qb)
        for i = 1:size(qb,1)
                qb[i,:] = rotateGFtoBF(qb[i,:], phi[j], theta[j], psi[j])
        end
        qbx = qb[:,1] .+ X[j]
        qby = qb[:,2] .+ Y[j]
        qbz = qb[:,3] .+ Z[j]
        plot(qbx,qby,qbz,lw=2 )
        plot!(X[1:j], Y[1:j], Z[1:j],lw=2)
        plot!(X[1:j], Y[1:j], 0 .*Z[1:j],lw=2, color=:black)
        scatter!(qbx[1:4],qby[1:4],qbz[1:4], markersize=10, markercolor=:blue, markeralpha=0.5)
        scatter!(qbx[3:4],qby[3:4],qbz[3:4], markersize=10, markercolor=:red, markeralpha=0.5)
        plot!(aspect_ratio=:equal, xlim=(-plotlim, plotlim), ylim=(-plotlim, plotlim), zlim=(-0.1*plotlim, plotlim), legend=:false)
        plot!(xlabel = "x", ylabel = "y")
        plot!(size = (1064, 768))
end


Quad, t_plot, U1, U2, U3, U4, X, Y, Z, phi, theta, psi  = QuadSim(-4,-3,-2,0)
@show Quad.X Quad.Y Quad.Z

l = Quad.l


frames = 30
animtime = Quad.sim_time
numframes = animtime*frames
skip = Int64(div(size(X,1),numframes))
plotlim = 5

if true

        anim = @animate for j = 1:skip:size(X,1)
                qb = [  l l 0 ;
                        -l l 0 ;
                        -l -l 0 ;
                        l -l 0 ;
                        l l 0 ]
                QuadPlot(j, plotlim, qb)
                progress = div((j+2)*100,size(X,1))
                @show progress
        end
        gif(anim, "quad.gif", fps = frames)

end
p1 = plot(t_plot,U1,label="F")
p1 = plot!(t_plot,U2,label="Tx")
p1 = plot!(t_plot,U3,label="Ty")
p1 = plot!(t_plot,U4,label="Tz")

p2 = plot(t_plot,X,label="X")
p2 = plot!(t_plot,Y,label="Y")
p2 = plot!(t_plot,Z,label="Z")

p3 = plot(t_plot, phi, label = "roll")
p3 = plot!(t_plot, theta, label = "pitch")
p3 = plot!(t_plot, psi, label = "yaw")

plot(p1, p2, p3, layout=(3,1))
