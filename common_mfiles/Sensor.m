% ------------------------------------------%
% ---------- Class def - sensor ------------%
% ------------------------------------------%
classdef Sensor
    properties
        name = 'DefaultName'
        x = 0 
        y = 0
        z = 0
        AOA = 0
        AOAMean = 0
        AOAVariance = 0
        TOA = 0
        TDOA = 0
        TOAMean = 0
        TOAVariance = 0
        position = [0; 0; 0]
        c = 299792458;% speed of light [m/s]
    end

    methods
       
        function obj = Sensor(name, x, y, z, AOA, AOAMean, AOAVariance, TOA, TDOA, TOAMean, TOAVariance)
            if nargin > 0
                obj.name = name;
                obj.x = x;
                obj.y = y;
                obj.z = z;
                obj.AOA = AOA;
                obj.AOAMean = AOAMean;
                obj.AOAVariance = AOAVariance;
                obj.TOA = TOA;
                obj.TDOA = TDOA;
                obj.TOAMean = TOAMean;
                obj.TOAVariance = TOAVariance;
                obj.position = [x; y; z];
            end
        end
        % ------------------------------------%
        % --- Getters for all properties: ----%
        % ------------------------------------%
        function value = get.x(obj)
            value = obj.x;
        end

        function value = get.y(obj)
            value = obj.y;
        end
        
        function value = get.z(obj)
            value = obj.z;
        end
        
        function value = get.AOA(obj)
            value = obj.AOA;
        end
        
        function value = get.AOAMean(obj)
            value = obj.AOAMean;
        end

        function value = get.AOAVariance(obj)
            value = obj.AOAVariance;
        end

        function value = get.TOAMean(obj)
            value = obj.TOAMean;
        end

        function value = get.TOAVariance(obj)
            value = obj.TOAVariance;
        end
        
        function value = get.position(obj)
            value = [obj.x; obj.y; obj.z];
        end
        
        % ------------------------------------%
        
        % ------------------------------------%
        % --- Setters for all properties: ----%
        % ------------------------------------%
        function obj = set.x(obj, value)
            obj.x = value;
        end
        
        function obj = set.y(obj, value)
            obj.y = value;
        end
        function obj = set.z(obj, value)
            obj.z = value;
        end
        
        function obj = set.AOA(obj, value)
            obj.AOA = value;
        end
        
        function obj = set.AOAMean(obj, value)
            obj.AOAMean = value;
        end
        
        function obj = set.AOAVariance(obj, value)
            obj.AOAVariance = value;
        end
        
        function obj = set.TOAMean(obj, value)
            obj.TOAMean = value;
        end
        
        function obj = set.TOAVariance(obj, value)
            obj.TOAVariance = value;
        end
        % ------------------------------------%
        
        
        % Calculate TOA of actual object
        function obj = calcTOA(obj, emitter)
            % Calculate time of arrival
            noise= mvnrnd(obj.TOAMean,obj.TOAVariance);
            rec = [obj.x obj.y obj.z];
            obj.TOA = norm(emitter - rec)/obj.c + noise;
        end
        
        % Calculate TDOA of actual object
        function TDOA = calcTDOA(obj, emitter, receiver_ref)
            % receiver_ref =  receiver against which everything counts (reference)
            %              - it is also object
            % Positioning from Time Difference Of Arrival - (2.6) equation
            %     \deltaT and h function
            % ... parameterized by the receivers (antennas) and the speed of light c
            noise= (obj.TOAVariance*randn(1) + obj.TOAMean) - (receiver_ref.TOAVariance*randn(1) + receiver_ref.TOAMean); % v_{ij} = v_i - v_j
            rec = [obj.x obj.y obj.z]';
            ref = [receiver_ref.x receiver_ref.y receiver_ref.z]';
            Meas = norm(emitter - rec)/obj.c - norm(emitter - ref)/obj.c;
            obj.TDOA = Meas + noise;
            TDOA = obj.TDOA;
        end
        
        % TDOA without noise
        function Meas = calcTDOAMeas(obj, emitter, receiver_ref)
            % receiver_ref =  receiver against which everything counts (reference)
            %              - it is also object
            % ... parameterized by the receivers (antennas) and the speed of light c
            rec = [obj.x obj.y obj.z]';
            ref = [receiver_ref.x receiver_ref.y receiver_ref.z]';
            Meas = norm(emitter - rec)/obj.c - norm(emitter - ref)/obj.c;
        end
        
        function [az, el] = calcAOA(obj, emitter)
            % AOA azimut and elevation calculation 
            rec = [obj.x obj.y obj.z]'; % receiver = actual sensor
            noise= obj.AOAVariance*randn(1) + obj.AOAMean;
            
            
            % Azimut: atan( (Py_e - Py_r) / (Px_e - Px_r) ) + noise
            az = atan2( (rec(2) - emitter(2)), (rec(1) - emitter(1)) ) + noise;
            % Normalize azimut to [0 2pi]
            az = wrapTo2Pi(az);
            
            % Elevation: atan( (Pz_e - Pz_r) / sqrt( (Py_e - Py_r)^2 + (Px_e - Px_r)^2 ) + noise
            el = atan2( (rec(3) - emitter(3)), sqrt( (rec(2) - emitter(2))^2 + (rec(1) - emitter(1))^2 ) ) + noise;
            % Normalize elevation to [0 2pi]
            el = wrapTo2Pi(el);
        end
        
        function [az, el] = calcAOAMeas(obj, emitter)
            % AOA azimut and elevation calculation 
            rec = [obj.x obj.y obj.z]'; % receiver = actual sensor
            
            % Azimut: atan( (Py_e - Py_r) / (Px_e - Px_r) ) + noise
            az = atan2( (rec(2) - emitter(2)), (rec(1) - emitter(1)) );
            % Normalize azimut to [0 2pi]
            az = wrapTo2Pi(az);
            
            % Elevation: atan( (Pz_e - Pz_r) / sqrt( (Py_e - Py_r)^2 + (Px_e - Px_r)^2 ) + noise
            el = atan2( (rec(3) - emitter(3)), sqrt( (rec(2) - emitter(2))^2 + (rec(1) - emitter(1))^2 ) );
            % Normalize elevation to [0 2pi]
            el = wrapTo2Pi(el);
        end
        
        % Calculate Jacobian of TDOA with respect to sensor's position
        function J = calcTDOAJacobian(obj, emitter, receiver_ref)
            rec = [obj.x; obj.y; obj.z];
            ref = [receiver_ref.x; receiver_ref.y; receiver_ref.z];
            J = (emitter - rec)'/norm(emitter - rec) - (emitter - ref)'/norm(emitter - ref); % Derrivate measurement eq
            J./sqrt(obj.TOAVariance^2);
            J = J/obj.c;
        end
        
        
        % Calculate the Jacobian matrix for AOA measurements
        function J = calcAOAJacobian(obj, emitter)       
            rec = [obj.x; obj.y; obj.z]; % Receiver's position
            
            % Partial derivatives of azimuth with respect to position
            % Azimuth = atan(DeltaY/DeltaX)
            diff = rec - emitter; % [deltaX, deltaY, deltaZ]
            dAz_dx = diff(2) / ((diff(1)^2 + diff(2)^2));
            % It is analytical for y 
            dAz_dy = - diff(1) / (diff(1)^2 + diff(2)^2);
            % Z axis is not included in calculation
            dAz_dz = 0;
            
            % Partial derivatives of elevation with respect to position
            % atan( deltaZ/sqrt(deltaY^2 + deltaX^2))
            r = (diff(1)^2 + diff(2)^2 + diff(3)^2); 
            dEl_dx = (diff(1) * diff(3)) / (r * sqrt(diff(1)^2 + diff(2)^2));
            dEl_dy = (diff(2) * diff(3)) / (r * sqrt(diff(1)^2 + diff(2)^2));
            dEl_dz = - (diff(1)^2 + diff(2)^2) / (r * sqrt(diff(1)^2 + diff(2)^2));
            
            % Assemble Jacobian
            J = [dAz_dx, dAz_dy, dAz_dz;
                 dEl_dx, dEl_dy, dEl_dz];
        end
        
    end
end
