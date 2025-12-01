using LinearAlgebra
using Random
using Distributions
using Plots
using Statistics
using Printf



"""
Parámetros físicos de la partícula y el medio
"""
struct ParticleParameters
    radius::Float64              # Radio de la partícula (m)
    mass::Float64               # Masa de la partícula (kg)
    magnetic_moment::Vector{Float64}  # Momento magnético dipolar (A·m²)
    temperature::Float64        # Temperatura (K)
    viscosity::Float64          # Viscosidad del fluido (Pa·s)
    D_t::Float64               # Coeficiente de difusión traslacional (m²/s)
    D_r::Float64               # Coeficiente de difusión rotacional (rad²/s)
end

"""
Constructor de parámetros con cálculo automático de coeficientes de difusión
"""
function ParticleParameters(radius, temperature, magnetic_moment_magnitude)
    # Constantes físicas
    k_B = 1.380649e-23  # Constante de Boltzmann (J/K)
    
    # Viscosidad del agua a 300 K
    viscosity = 8.5e-4  # Pa·s
    
    # Densidad aproximada de partícula magnética (Fe₃O₄)
    density = 5200.0  # kg/m³
    
    # Masa de la partícula
    volume = (4/3) * π * radius^3
    mass = density * volume
    
    # Coeficientes de difusión (Stokes-Einstein)
    D_t = k_B * temperature / (6 * π * viscosity * radius)
    D_r = k_B * temperature / (8 * π * viscosity * radius^3)
    
    # Momento magnético en dirección z del marco del cuerpo
    magnetic_moment = [0.0, 0.0, magnetic_moment_magnitude]
    
    return ParticleParameters(radius, mass, magnetic_moment, temperature, 
                            viscosity, D_t, D_r)
end

"""
Estado de la partícula (posición y orientación)
"""
mutable struct ParticleState
    position::Vector{Float64}      # Posición (x, y, z)
    orientation::Vector{Float64}   # Cuaternión (q0, q1, q2, q3)
    velocity::Vector{Float64}      # Velocidad (opcional, para análisis)
    angular_velocity::Vector{Float64}  # Velocidad angular (opcional)
end

"""
Inicializar estado de la partícula
"""
function initialize_particle()
    position = [0.0, 0.0, 0.0]
    # Cuaternión identidad (sin rotación inicial)
    orientation = [1.0, 0.0, 0.0, 0.0]
    velocity = [0.0, 0.0, 0.0]
    angular_velocity = [0.0, 0.0, 0.0]
    
    return ParticleState(position, orientation, velocity, angular_velocity)
end

"""
Normalizar cuaternión
"""
function normalize_quaternion!(q::Vector{Float64})
    norm_q = sqrt(sum(q.^2))
    q ./= norm_q
end

"""
Convertir cuaternión a matriz de rotación
"""
function quaternion_to_rotation_matrix(q::Vector{Float64})
    q0, q1, q2, q3 = q
    
    R = [
        1 - 2(q2^2 + q3^2)    2(q1*q2 - q0*q3)      2(q1*q3 + q0*q2);
        2(q1*q2 + q0*q3)      1 - 2(q1^2 + q3^2)    2(q2*q3 - q0*q1);
        2(q1*q3 - q0*q2)      2(q2*q3 + q0*q1)      1 - 2(q1^2 + q2^2)
    ]
    
    return R
end

"""
Calcular energía potencial magnética: U = -μ·B
"""
function magnetic_energy(state::ParticleState, params::ParticleParameters, 
                        B_field::Vector{Float64})
    # Rotar momento magnético al sistema de laboratorio
    R = quaternion_to_rotation_matrix(state.orientation)
    μ_lab = R * params.magnetic_moment
    
    # Energía: U = -μ·B
    U = -dot(μ_lab, B_field)
    
    return U
end

"""
Calcular torque magnético: τ = μ × B
"""
function magnetic_torque(state::ParticleState, params::ParticleParameters,
                        B_field::Vector{Float64})
    # Rotar momento magnético al sistema de laboratorio
    R = quaternion_to_rotation_matrix(state.orientation)
    μ_lab = R * params.magnetic_moment
    
    # Torque: τ = μ × B
    τ = cross(μ_lab, B_field)
    
    return τ
end

"""
Paso de integración para dinámica browniana traslacional
Ecuación de Langevin: m*dv/dt = -ξ*v + F + F_B
Aproximación de alto amortiguamiento (régimen browniano)
"""
function translational_step!(state::ParticleState, params::ParticleParameters,
                           B_field::Vector{Float64}, dt::Float64, rng::AbstractRNG)
    # Coeficiente de fricción
    ξ = 6 * π * params.viscosity * params.radius
    
    # Fuerza determinística (en este caso, solo gravedad o fuerza magnética si hay gradiente)
    # Para campo uniforme, no hay fuerza neta, solo torque
    F_det = [0.0, 0.0, 0.0]
    
    # Ruido térmico traslacional (distribución gaussiana)
    σ_t = sqrt(2 * params.D_t * dt)
    ΔW_t = σ_t * randn(rng, 3)
    
    # Actualización de posición (régimen sobreamortiguado)
    # dr = (F_det/ξ)*dt + ΔW_t
    state.position .+= (F_det / ξ) * dt .+ ΔW_t
end

"""
Paso de integración para dinámica browniana rotacional
"""
function rotational_step!(state::ParticleState, params::ParticleParameters,
                         B_field::Vector{Float64}, dt::Float64, rng::AbstractRNG)
    # Calcular torque magnético
    τ = magnetic_torque(state, params, B_field)
    
    # Coeficiente de fricción rotacional
    ξ_r = 8 * π * params.viscosity * params.radius^3
    
    # Ruido térmico rotacional
    σ_r = sqrt(2 * params.D_r * dt)
    ΔW_r = σ_r * randn(rng, 3)
    
    # Velocidad angular en régimen sobreamortiguado
    # ω = τ/ξ_r + ruido
    ω = (τ / ξ_r) .+ ΔW_r / dt
    
    # Actualizar cuaternión usando velocidad angular
    # dq/dt = (1/2) * Ω * q, donde Ω es la matriz de ω
    q = state.orientation
    Ω = [
        0      -ω[1]  -ω[2]  -ω[3];
        ω[1]    0      ω[3]  -ω[2];
        ω[2]   -ω[3]   0      ω[1];
        ω[3]    ω[2]  -ω[1]   0
    ]
    
    dq = 0.5 * Ω * q * dt
    state.orientation .+= dq
    
    # Renormalizar cuaternión
    normalize_quaternion!(state.orientation)
    
    # Guardar velocidad angular para análisis
    state.angular_velocity .= ω
end

"""
Simulación completa de dinámica browniana
"""
function simulate_brownian_dynamics(params::ParticleParameters, 
                                   B_field::Vector{Float64},
                                   t_total::Float64, dt::Float64;
                                   save_interval::Int=10)
    # Inicializar
    rng = MersenneTwister(123)
    state = initialize_particle()
    
    # Arrays para guardar resultados
    n_steps = Int(round(t_total / dt))
    n_saved = Int(round(n_steps / save_interval))
    
    times = zeros(n_saved)
    positions = zeros(n_saved, 3)
    orientations = zeros(n_saved, 4)
    energies = zeros(n_saved)
    
    # Simulación
    save_counter = 1
    for step in 1:n_steps
        t = step * dt
        
        # Actualizar estado
        translational_step!(state, params, B_field, dt, rng)
        rotational_step!(state, params, B_field, dt, rng)
        
        # Guardar datos
        if step % save_interval == 0
            times[save_counter] = t
            positions[save_counter, :] = state.position
            orientations[save_counter, :] = state.orientation
            energies[save_counter] = magnetic_energy(state, params, B_field)
            save_counter += 1
        end
    end
    
    return times, positions, orientations, energies
end

"""
Calcular Mean Square Displacement (MSD)
"""
function calculate_msd(positions::Matrix{Float64}, times::Vector{Float64})
    n = size(positions, 1)
    max_lag = div(n, 4)  # Calcular hasta 1/4 del tiempo total
    
    msd = zeros(max_lag)
    time_lags = zeros(max_lag)
    
    for lag in 1:max_lag
        displacements = 0.0
        count = 0
        
        for i in 1:(n - lag)
            dr = positions[i + lag, :] - positions[i, :]
            displacements += sum(dr.^2)
            count += 1
        end
        
        msd[lag] = displacements / count
        time_lags[lag] = times[min(1 + lag, n)] - times[1]
    end
    
    return time_lags, msd
end

"""
Generar todas las gráficas para reporte (de mec.estadistica)
"""
function generate_report_plots(times, positions, orientations, energies, params, B_field)
    plots_array = []
    
    # 1. Trayectoria 3D
    p1 = plot(positions[:, 1], positions[:, 2], positions[:, 3],
             xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
             title="Trayectoria Browniana 3D",
             legend=false, linewidth=1.5, alpha=0.7,
             camera=(45, 30))
    scatter!(p1, [positions[1, 1]], [positions[1, 2]], [positions[1, 3]], 
            marker=:circle, markersize=6, color=:green, label="Inicio")
    scatter!(p1, [positions[end, 1]], [positions[end, 2]], [positions[end, 3]], 
            marker=:star, markersize=8, color=:red, label="Fin")
    push!(plots_array, p1)
    
    # 2. Proyecciones 2D
    p2 = plot(positions[:, 1], positions[:, 2],
             xlabel="x (m)", ylabel="y (m)",
             title="Proyección XY",
             legend=false, linewidth=1.5, alpha=0.7)
    push!(plots_array, p2)
    
    # 3. Posición vs tiempo
    p3 = plot(times, positions[:, 1], label="x", linewidth=1.5)
    plot!(p3, times, positions[:, 2], label="y", linewidth=1.5)
    plot!(p3, times, positions[:, 3], label="z", linewidth=1.5)
    xlabel!(p3, "Tiempo (s)")
    ylabel!(p3, "Posición (m)")
    title!(p3, "Posición vs Tiempo")
    push!(plots_array, p3)
    
    # 4. Mean Square Displacement
    time_lags, msd = calculate_msd(positions, times)
    p4 = plot(time_lags, msd, 
             xlabel="Intervalo de tiempo (s)", ylabel="MSD (m²)",
             title="Mean Square Displacement",
             legend=false, linewidth=2, marker=:circle, markersize=3)
    
    # Ajuste lineal para verificar comportamiento difusivo
    if length(time_lags) > 10
        idx_fit = 5:min(50, length(time_lags))
        A = hcat(time_lags[idx_fit], ones(length(idx_fit)))
        coef = A \ msd[idx_fit]
        D_fitted = coef[1] / 6  # MSD = 6*D*t en 3D
        
        plot!(p4, time_lags, coef[1] * time_lags .+ coef[2], 
             label="Ajuste lineal", linestyle=:dash, linewidth=2)
        annotate!(p4, [(time_lags[end]*0.5, maximum(msd)*0.8, 
                       text(@sprintf("D_ajustado = %.2e m²/s", D_fitted), 10))])
    end
    push!(plots_array, p4)
    
    # 5. Distribución de desplazamientos
    all_displacements = []
    for lag in [10, 50, 100]
        if lag < size(positions, 1)
            for i in 1:(size(positions, 1) - lag)
                dr = norm(positions[i + lag, :] - positions[i, :])
                push!(all_displacements, dr)
            end
        end
    end
    
    p5 = histogram(all_displacements, bins=30, normalize=:pdf,
                  xlabel="Desplazamiento (m)", ylabel="Densidad de probabilidad",
                  title="Distribución de Desplazamientos",
                  legend=false, alpha=0.7)
    push!(plots_array, p5)
    
    # 6. Energía magnética vs tiempo
    p6 = plot(times, energies * 1e21,  # Convertir a 10^-21 J para mejor visualización
             xlabel="Tiempo (s)", ylabel="Energía magnética (×10⁻²¹ J)",
             title="Energía Potencial Magnética",
             legend=false, linewidth=1.5)
    hline!(p6, [mean(energies) * 1e21], linestyle=:dash, linewidth=2, color=:red, 
          label="Promedio")
    push!(plots_array, p6)
    
    # 7. Orientación del momento magnético
    μ_directions = zeros(size(orientations, 1), 3)
    for i in 1:size(orientations, 1)
        R = quaternion_to_rotation_matrix(orientations[i, :])
        μ_directions[i, :] = R * params.magnetic_moment
    end
    
    # Normalizar para visualización en esfera unitaria
    μ_norm = μ_directions ./ norm(params.magnetic_moment)
    
    p7 = plot(xlabel="μₓ/|μ|", ylabel="μᵧ/|μ|", zlabel="μᵤ/|μ|",
             title="Evolución de la Orientación del Momento Magnético",
             legend=false, camera=(45, 30))
    plot!(p7, μ_norm[:, 1], μ_norm[:, 2], μ_norm[:, 3],
         linewidth=1.5, alpha=0.7)
    scatter!(p7, [μ_norm[1, 1]], [μ_norm[1, 2]], [μ_norm[1, 3]],
            marker=:circle, markersize=6, color=:green)
    scatter!(p7, [μ_norm[end, 1]], [μ_norm[end, 2]], [μ_norm[end, 3]],
            marker=:star, markersize=8, color=:red)
    push!(plots_array, p7)
    
    return plots_array
end

"""
Imprimir información adicional 
"""
function print_system_info(params::ParticleParameters, B_field::Vector{Float64})
    println("="^70)
    println("PARÁMETROS DEL SISTEMA - PARTÍCULA COLOIDAL MAGNÉTICA")
    println("="^70)
    println()
    println("Parámetros de la partícula:")
    println("  Radio:                    ", @sprintf("%.2e m (%.1f nm)", 
            params.radius, params.radius * 1e9))
    println("  Masa:                     ", @sprintf("%.2e kg", params.mass))
    println("  Momento magnético:        ", @sprintf("%.2e A·m²", 
            norm(params.magnetic_moment)))
    println()
    println("Parámetros del medio:")
    println("  Temperatura:              ", @sprintf("%.1f K", params.temperature))
    println("  Viscosidad:               ", @sprintf("%.2e Pa·s", params.viscosity))
    println()
    println("Coeficientes de difusión:")
    println("  D_traslacional:           ", @sprintf("%.2e m²/s", params.D_t))
    println("  D_rotacional:             ", @sprintf("%.2e rad²/s", params.D_r))
    println()
    println("Campo magnético aplicado:")
    println("  B = [", @sprintf("%.2e, %.2e, %.2e", B_field...), "] T")
    println("  |B| = ", @sprintf("%.2e T", norm(B_field)))
    println()
    
    # Números adimensionales
    k_B = 1.380649e-23
    E_thermal = k_B * params.temperature
    E_magnetic = norm(params.magnetic_moment) * norm(B_field)
    
    println("Relación de energías:")
    println("  E_térmica (k_B*T):        ", @sprintf("%.2e J", E_thermal))
    println("  E_magnética (μ·B):        ", @sprintf("%.2e J", E_magnetic))
    println("  E_magnética/E_térmica:    ", @sprintf("%.2f", E_magnetic/E_thermal))
    println()
    println("="^70)
end
#

println("\n Iniciando Simulación de Dinámica Browniana de Partícula Magnética\n")

# Definir parámetros
radius = 100e-9  # 100 nm
temperature = 300.0  # K
magnetic_moment_magnitude = 1e-18  # A·m²

# Crear objeto de parámetros
params = ParticleParameters(radius, temperature, magnetic_moment_magnitude)

# Campo magnético aplicado (T)
B_field = [0.0, 0.0, 0.01]  # 10 mT en dirección z

# Imprimir información del sistema
print_system_info(params, B_field)

# Parámetros de simulación
t_total = 0.1  # Tiempo total de simulación (s)
dt = 1e-6      # Paso de tiempo (μs)
save_interval = 100  # Guardar cada N pasos

println("\n Parámetros de simulación:")
println("  Tiempo total:             ", t_total, " s")
println("  Paso de tiempo:           ", dt * 1e6, " μs")
println("  Número de pasos:          ", Int(round(t_total / dt)))
println("  Intervalo de guardado:    ", save_interval)
println()

# Ejecutar simulación
println("  Ejecutando simulación...\n")
@time times, positions, orientations, energies = simulate_brownian_dynamics(
    params, B_field, t_total, dt; save_interval=save_interval
)

println(" Simulación completada!\n")

# Generar gráficas
println(" Generando gráficas...\n")
plots_array = generate_report_plots(times, positions, orientations, energies, params, B_field)

# Crear carpeta para guardar resultados
output_dir = "resultados_browniano"
if !isdir(output_dir)
    mkdir(output_dir)
end

# Guardar gráficas individuales
println("Guardando gráficas individuales...")
savefig(plots_array[1], joinpath(output_dir, "01_trayectoria_3D.png"))
savefig(plots_array[2], joinpath(output_dir, "02_proyeccion_XY.png"))
savefig(plots_array[3], joinpath(output_dir, "03_posicion_vs_tiempo.png"))
savefig(plots_array[4], joinpath(output_dir, "04_MSD.png"))
savefig(plots_array[5], joinpath(output_dir, "05_distribucion_desplazamientos.png"))
savefig(plots_array[6], joinpath(output_dir, "06_energia_magnetica.png"))
savefig(plots_array[7], joinpath(output_dir, "07_orientacion_momento.png"))

# Guardar gráficas combinadas
println(" Guardando gráficas combinadas...")
p_combined1 = plot(plots_array[1:4]..., layout=(2, 2), size=(1200, 900))
savefig(p_combined1, joinpath(output_dir, "figuras_combinadas_1-4.png"))

p_combined2 = plot(plots_array[5:7]..., layout=(1, 3), size=(1400, 400))
savefig(p_combined2, joinpath(output_dir, "figuras_combinadas_5-7.png"))

# Mostrar gráficas en pantalla
display(p_combined1)
display(p_combined2)

println("\n Gráficas guardadas en la carpeta: ", output_dir)
println(" Análisis completado. Gráficas listas para el reporte.")
println("\n Archivos generados:")
println("  - 7 gráficas individuales (01-07)")
println("  - 2 gráficas combinadas")
println("  - Total: 9 archivos PNG")