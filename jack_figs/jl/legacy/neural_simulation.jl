using OrdinaryDiffEq, DSP

# Heaviside function
heaviside(x) = x >= 0 ? 1.0 : 0.0

# ODE system function
function neural_network!(du, u, p, t)
    # Unpack state variables
    V0, V1, V2, V3, V4 = u[1:5]
    x0, x1, x2, x3, x4 = u[6:10]
    Ca0, Ca1, Ca2, Ca3, Ca4 = u[11:15]
    h0, h1, h2, h3, h4 = u[16:20]
    n0, n1, n2, n3, n4 = u[21:25]
    y1, y2 = u[26:27]
    s0, s1, s2, s3, s4 = u[28:32]
    s12, s21, s34, s43 = u[33:36]
    sm = u[37]
    
    # Unpack parameters
    gm, gelec, g0, g41, g32, g14, g23, g34, g43, g21, g12 = p[1:11]
    alpha1, beta1, alpham, betam, alphax, betax, alphai, betai, alpha2, beta2, alpha3, beta3 = p[12:23]
    scale1, scalem, scalei, scale2, scale3 = p[24:28]
    Ca_shift0, Ca_shift1, Ca_shift2, Ca_shift3, Ca_shift4, x_shift, Vhh = p[29:35]
    Iapp, t1, t2 = p[36:38]
    
    # Convert time from seconds to msec for compatibility
    tt = t * 1000.0
    
    scalex2 = 1.0
    
    # Voltage equations
    du[1] = Iapp*heaviside(tt-t1)*heaviside(t2+1250-tt)+4*((0.1*(50-(127*V0/105+8265/105))/(exp((50 - (127*V0/105 +8265/105))/10) - 1))/((0.1*(50 - (127*V0/105 + 8265/105))/(exp((50 - (127*V0/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V0/105 + 8265/105))/18))))^3*h0*(30 - V0) + 0.3*n0^4*(-75 - V0)+0.01*x0*(30-V0) +0.03*Ca0/(0.5 + Ca0)*(-75 - V0)+0.003*(-40 - V0)
    
    du[2] = 4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) +0.03*Ca1/(0.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) -g41*(V1-30)*s4/scalex2 - g12*(V1+80)*s21/scale2 + gelec*(V4-V1)
    
    du[3] = 4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) +0.03*Ca2/(0.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) -g32*(V2-30)*s3/scalex2 - g12*(V2+80)*s12/scale2 + gelec*(V3-V2)
    
    du[4] = 4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) +0.03*Ca3/(0.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3) -g23*(V3+80)*s2/scalei - g34*(V3+80)*s43/scale3 + gelec*(V2-V3) - g0*(V3-30)*s0/scale1
    
    du[5] = 4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) +0.03*Ca4/(0.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4) -g14*(V4+80)*s1/scalei - g34*(V4+80)*s34/scale3 + gelec*(V1-V4) - g0*(V4-30)*s0/scale1
    
    # Calcium dynamics
    du[11] = 0.000012*(0.0085*x0*(140-V0+Ca_shift0)-Ca0)
    du[12] = 0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1)
    du[13] = 0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2)
    du[14] = 0.0003*(0.0085*x3*(140-V3+Ca_shift3)-Ca3)
    du[15] = 0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4)
    
    # x dynamics
    du[6] = (((1/(exp(0.15*(-V0-50+x_shift))+1))-x0)/100)
    du[7] = (((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100)
    du[8] = (((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100)
    du[9] = (((1/(exp(0.15*(-V3-50+x_shift))+1))-x3)/100)
    du[10] = (((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100)
    
    # h dynamics
    du[16] = (((1-h0)*(0.07*exp((25 - (127*V0/105 + 8265/105))/20))-h0*(1.0/(1 + exp((55 - (127*V0/105 + 8265/105))/10))))/12.5)
    du[17] = (((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5)
    du[18] = (((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5)
    du[19] = (((1-h3)*(0.07*exp((25 - (127*V3/105 + 8265/105))/20))-h3*(1.0/(1 + exp((55 - (127*V3/105 + 8265/105))/10))))/12.5)
    du[20] = (((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5)
    
    # n dynamics
    du[21] = (((1-n0)*(0.01*(55 - (127*V0/105 + 8265/105))/(exp((55 - (127*V0/105 + 8265/105))/10) - 1))-n0*(0.125*exp((45 - (127*V0/105 + 8265/105))/80)))/12.5)
    du[22] = (((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5)
    du[23] = (((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5)
    du[24] = (((1-n3)*(0.01*(55 - (127*V3/105 + 8265/105))/(exp((55 - (127*V3/105 + 8265/105))/10) - 1))-n3*(0.125*exp((45 - (127*V3/105 + 8265/105))/80)))/12.5)
    du[25] = (((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5)
    
    # y dynamics
    du[26] = 0.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2)))
    du[27] = 0.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2)))
    
    # Synaptic dynamics
    du[28] = alpha1*(1-s0)/(1+exp(-20*(V0+20)))-beta1*s0
    du[29] = alphai*s1*(1-s1)/(1+exp(-20*(V1+20)))-betai*(s1-0.0001)
    du[30] = alphai*s2*(1-s2)/(1+exp(-20*(V2+20)))-betai*(s2-0.0001)
    du[31] = alphax*(1+gm*sm)/scalem*s3*(1-s3)/(1+exp(-20*(V3+20)))-betax*(s3-0.0001)
    du[32] = alphax*(1+gm*sm)/scalem*s4*(1-s4)/(1+exp(-20*(V4+20)))-betax*(s4-0.0001)
    
    du[33] = alpha2*(1-s12)/(1+exp(-20*(V1+20)))-beta2*s12
    du[34] = alpha2*(1-s21)/(1+exp(-20*(V2+20)))-beta2*s21
    du[35] = alpha3*(1-s34)/(1+exp(-20*(V3+20)))-beta3*s34
    du[36] = alpha3*(1-s43)/(1+exp(-20*(V4+20)))-beta3*s43
    
    du[37] = alpham*sm*(1-sm)/(1+exp(-20*(V0+20)))-betam*(sm-0.0001)
end

# Function to generate band-limited noise (simplified version)
function bandlimnoise(cutoff_freq, timestep, total_time)
    nt = Int(total_time / timestep)
    # Generate white noise
    white_noise = randn(nt)
    
    # Simple low-pass filter approximation
    # For more accurate filtering, you could use DSP.jl
    filtered_noise = copy(white_noise)
    alpha_filter = cutoff_freq * timestep * 2π
    for i in 2:nt
        filtered_noise[i] = alpha_filter * white_noise[i] + (1 - alpha_filter) * filtered_noise[i-1]
    end
    
    return filtered_noise, nt
end

# Function to find peaks (simplified version of MATLAB's findpeaks)
function findpeaks(data; min_height=-Inf, min_distance=1)
    peaks = Int[]
    peak_values = Float64[]
    
    for i in 2:(length(data)-1)
        if data[i] > data[i-1] && data[i] > data[i+1] && data[i] >= min_height
            # Check minimum distance constraint
            if isempty(peaks) || (i - peaks[end]) >= min_distance
                push!(peaks, i)
                push!(peak_values, data[i])
            end
        end
    end
    
    return peak_values, peaks
end

# Main simulation function
function run_neural_simulation(gm; tmax=200.0, tsamp=2.0)
    # Simulation parameters
    gelec = 0.002
    
    # Excitation from Si1
    alpha1 = 0.01
    beta1 = 0.002
    scale1 = alpha1/(alpha1+beta1)
    g0 = 0.002
    
    # Excitation from Si3s
    alpham = 0.005
    betam = 0.0001
    scalem = (alpham-betam)/alpham
    
    alphax = 0.0001
    betax = 0.001
    scalex = (alphax-betax)/alphax
    
    g41 = 0.001
    g32 = g41
    
    # Inhibition from Si2
    g14 = 0.01
    g23 = g14
    alphai = 0.015
    betai = 0.0005
    scalei = (alphai-betai)/alphai
    
    # Inhibition between Si3s
    alpha3 = 0.01
    beta3 = 0.002
    scale3 = alpha3/(alpha3+beta3)
    
    g34 = 0.005
    g43 = g34
    
    # Inhibition between Si2s
    alpha2 = 0.01
    beta2 = 0.01
    scale2 = alpha2/(alpha2+beta2)
    
    g21 = 0.01
    g12 = g21
    
    # Parameters
    Ca_shift0 = -10.0
    Ca_shift4 = -25.0
    Ca_shift1 = -25.0
    Ca_shift3 = Ca_shift4
    Ca_shift2 = Ca_shift1
    x_shift = -4.0
    
    Iapp = 0.0
    t1 = 1*1000.0 # this time is in msec
    t2 = 35*1000.0 # this time is in msec
    
    # H-current      
    gh = 0.0005
    Vhh = -53.0
    
    # Set up initial conditions
    u0 = zeros(37)
    u0[1:5] = [-44.0, -44.0, -44.0, -44.0, -44.0]  # V0, V1, V2, V3, V4
    u0[6:10] = [0.0, 0.6, 0.6, 0.5, 0.6]  # x0, x1, x2, x3, x4
    u0[11:15] = [0.3, 1.0, 1.0, 1.1, 1.0]  # Ca0, Ca1, Ca2, Ca3, Ca4
    u0[16:20] = [0.0, 0.0, 0.0, 0.0, 0.0]  # h0, h1, h2, h3, h4
    u0[21:25] = [0.0, 0.0, 0.0, 0.0, 0.0]  # n0, n1, n2, n3, n4
    u0[26:27] = [0.0, 0.0]  # y1, y2
    u0[28:32] = [0.0, 0.0, 0.0, 0.0, 0.0]  # s0, s1, s2, s3, s4
    u0[33:36] = [0.0, 0.0, 0.1, 0.1]  # s12, s21, s34, s43
    u0[37] = 0.99  # sm
    
    # Set up parameters
    p = [gm, gelec, g0, g41, g32, g14, g23, g34, g43, g21, g12,
         alpha1, beta1, alpham, betam, alphax, betax, alphai, betai, alpha2, beta2, alpha3, beta3,
         scale1, scalem, scalei, scale2, scale3,
         Ca_shift0, Ca_shift1, Ca_shift2, Ca_shift3, Ca_shift4, x_shift, Vhh,
         Iapp, t1, t2]
    
    # Time span in seconds
    tspan = (0.0, tmax)
    
    # Create and solve ODE problem
    prob = ODEProblem(neural_network!, u0, tspan, p)
    sol = solve(prob, RK4(), saveat=tsamp, reltol=1e-6, abstol=1e-6)
    
    # Extract solution
    time = sol.t
    vv0 = [sol.u[i][1] for i in 1:length(sol.u)]
    vv1 = [sol.u[i][2] for i in 1:length(sol.u)]
    vv2 = [sol.u[i][3] for i in 1:length(sol.u)]
    vv3 = [sol.u[i][4] for i in 1:length(sol.u)]
    vv4 = [sol.u[i][5] for i in 1:length(sol.u)]
    
    return time, vv0, vv1, vv2, vv3, vv4, sol
end

# Burst detection function
function detect_bursts(time, voltage; burst_threshold=-20.0, min_interburst_interval=500.0, 
                      analysis_start_time=50.0, analysis_end_time=190.0, tsamp=2.0)
    
    # Find peaks above threshold
    min_peak_distance = Int(round(min_interburst_interval/tsamp))
    peaks, peak_locs = findpeaks(voltage, min_height=burst_threshold, min_distance=min_peak_distance)
    
    if length(peak_locs) > 1
        # Find peaks in analysis window
        analysis_indices = findall(x -> x >= analysis_start_time && x <= analysis_end_time, time[peak_locs])
        num_bursts = length(analysis_indices)
        analysis_duration = analysis_end_time - analysis_start_time
        
        burst_frequency = num_bursts / analysis_duration
        return burst_frequency, num_bursts, peak_locs, peaks
    else
        return 0.0, 0, peak_locs, peaks
    end
end
