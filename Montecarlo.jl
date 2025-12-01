
using Plots, Random, Statistics, LinearAlgebra, CSV, DataFrames, Printf

# Constantes físicas
const k_B = 1.380649e-23      # Constante de Boltzmann [J/K]
const μ_0 = 4π * 1e-7         # Permeabilidad del vacío [H/m]

# Propiedades del glóbulo rojo
struct RedBloodCell
    radius::Float64           # Radio [m]
    volume::Float64           # Volumen [m³]
    density::Float64          # Densidad [kg/m³]
    chi_magnetic::Float64     # Susceptibilidad magnética (adimensional)
end

function inicializar_globulo_rojo()
    radius = 3.8e-6           # 3.8 μm
    volume = 90e-15           # 90 femtolitros = 90×10⁻¹⁵ m³
    density = 1125.0          # kg/m³
    chi_magnetic = -9e-6      # Diamagnético
    return RedBloodCell(radius, volume, density, chi_magnetic)
end

# Parámetros de simulación
struct SimulationParams
    N_globulos::Int           # Número de glóbulos rojos
    N_equilibrio::Int         # Pasos de equilibración
    N_produccion::Int         # Pasos de producción
    H_magnitud::Float64       # Campo magnético [Tesla]
    T::Float64                # Temperatura [K]
    save_interval::Int        # Guardar cada N pasos
end


"""
Calcula el momento magnético inducido en el glóbulo rojo
μ = χ * V * H / μ₀
"""
function calcular_momento_magnetico(rbc::RedBloodCell, H::Float64)
    return abs(rbc.chi_magnetic * rbc.volume * H / μ_0)
end

"""
Energía de interacción magnética: φ(θ) = -μ⃗ · H⃗ = -μ*H*cos(θ)
θ: ángulo entre momento magnético y campo (en radianes)
"""
function energia_magnetica(μ::Float64, H::Float64, θ::Float64)
    return -μ * H * cos(θ)
end

"""
Propone un nuevo ángulo θ con un paso aleatorio
Usa condiciones de frontera periódicas correctas
"""
function proponer_nuevo_angulo(θ_actual::Float64, δθ_max::Float64=0.5)
    δθ = (2*rand() - 1) * δθ_max  # Paso aleatorio en [-δθ_max, δθ_max]
    θ_nuevo = θ_actual + δθ
    
    # Condiciones de frontera: reflexión en los extremos
    # Si θ < 0, reflejamos: θ → -θ
    # Si θ > π, reflejamos: θ → 2π - θ
    if θ_nuevo < 0
        θ_nuevo = -θ_nuevo
    elseif θ_nuevo > π
        θ_nuevo = 2π - θ_nuevo
    end
    
    # Asegurar que está en [0, π] (por si acaso)
    θ_nuevo = clamp(θ_nuevo, 0.0, π)
    
    return θ_nuevo
end


"""
Criterio de Metropolis: acepta o rechaza la nueva configuración
"""
function criterio_metropolis(ΔE::Float64, β::Float64)
    if ΔE <= 0
        return true  # Siempre acepta si disminuye energía
    else
        return rand() < exp(-β * ΔE)  # Acepta con probabilidad e^(-βΔE)
    end
end

"""
Realiza un paso de Monte Carlo para un glóbulo rojo
Retorna: (θ_nuevo, aceptado)
"""
function paso_montecarlo!(θ::Float64, μ::Float64, H::Float64, β::Float64, δθ_max::Float64=0.8)
    # Proponer nueva orientación con paso más grande
    θ_nuevo = proponer_nuevo_angulo(θ, δθ_max)
    
    # Calcular cambio de energía
    E_vieja = energia_magnetica(μ, H, θ)
    E_nueva = energia_magnetica(μ, H, θ_nuevo)
    ΔE = E_nueva - E_vieja
    
    # Aplicar criterio de Metropolis
    if criterio_metropolis(ΔE, β)
        return θ_nuevo, true
    else
        return θ, false
    end
end

"""
Simula el ensemble de N glóbulos rojos
"""
function simular_ensemble(params::SimulationParams)
    # Inicializar
    rbc = inicializar_globulo_rojo()
    μ = calcular_momento_magnetico(rbc, params.H_magnitud)
    β = 1.0 / (k_B * params.T)
    
    println("=" ^ 70)
    println("PARÁMETROS DE LA SIMULACIÓN")
    println("=" ^ 70)
    println("Glóbulos rojos: $(params.N_globulos)")
    println("Radio: $(rbc.radius*1e6) μm")
    println("Volumen: $(rbc.volume*1e15) fL")
    println("Susceptibilidad magnética χ: $(rbc.chi_magnetic)")
    println("Campo magnético H: $(params.H_magnitud) T")
    println("Temperatura T: $(params.T) K")
    println("Momento magnético μ: $(@sprintf("%.3e", μ)) A·m²")
    println("Energía térmica k_B*T: $(@sprintf("%.3e", k_B*params.T)) J")
    println("Energía magnética μ*H: $(@sprintf("%.3e", μ*params.H_magnitud)) J")
    println("Ratio μH/(k_B*T): $(@sprintf("%.3f", μ*params.H_magnitud/(k_B*params.T)))")
    println("=" ^ 70)
    
    # Estados iniciales: ángulos UNIFORMEMENTE DISTRIBUIDOS en [0, π]
    θ_globulos = π * rand(params.N_globulos)
    
    # DIAGNÓSTICO: Verificar inicialización
    println("\nDIAGNÓSTICO DE INICIALIZACIÓN:")
    println("θ mínimo inicial: $(@sprintf("%.4f", minimum(θ_globulos))) rad = $(@sprintf("%.2f", rad2deg(minimum(θ_globulos))))°")
    println("θ máximo inicial: $(@sprintf("%.4f", maximum(θ_globulos))) rad = $(@sprintf("%.2f", rad2deg(maximum(θ_globulos))))°")
    println("θ promedio inicial: $(@sprintf("%.4f", mean(θ_globulos))) rad = $(@sprintf("%.2f", rad2deg(mean(θ_globulos))))°")
    println("=" ^ 70)
    
    # Arrays para guardar resultados
    n_total = params.N_equilibrio + params.N_produccion
    n_guardados = div(n_total, params.save_interval)
    
    energias = zeros(n_guardados)
    aceptacion = zeros(n_guardados)
    configuraciones = zeros(n_guardados, params.N_globulos)
    
    # Simulación
    println("\nIniciando simulación...")
    idx_guardado = 1
    total_aceptados = 0
    total_intentos = 0
    
    for paso in 1:n_total
        # Un paso MC = intentar mover cada glóbulo una vez
        aceptados_paso = 0
        
        for i in 1:params.N_globulos
            θ_globulos[i], aceptado = paso_montecarlo!(θ_globulos[i], μ, params.H_magnitud, β)
            if aceptado
                aceptados_paso += 1
            end
        end
        
        total_aceptados += aceptados_paso
        total_intentos += params.N_globulos
        
        # Guardar datos
        if paso % params.save_interval == 0
            # Calcular energía total del sistema
            E_total = sum([energia_magnetica(μ, params.H_magnitud, θ) for θ in θ_globulos])
            
            energias[idx_guardado] = E_total
            aceptacion[idx_guardado] = aceptados_paso / params.N_globulos
            configuraciones[idx_guardado, :] = θ_globulos
            
            idx_guardado += 1
        end
        
        # Progreso
        if paso % 10000 == 0
            fase = paso <= params.N_equilibrio ? "EQUILIBRIO" : "PRODUCCIÓN"
            tasa_acept = total_aceptados / total_intentos * 100
            E_actual = sum([energia_magnetica(μ, params.H_magnitud, θ) for θ in θ_globulos])
            θ_prom = mean(θ_globulos)
            println("Paso $paso / $n_total [$fase] | Aceptación: $(@sprintf("%.1f", tasa_acept))% | E: $(@sprintf("%.3e", E_actual)) J | ⟨θ⟩: $(@sprintf("%.2f", rad2deg(θ_prom)))°")
        end
    end
    
    println("\nSimulación completada!")
    println("Tasa de aceptación global: $(@sprintf("%.1f", total_aceptados/total_intentos*100))%")
    
    return Dict(
        "energias" => energias,
        "aceptacion" => aceptacion,
        "configuraciones" => configuraciones,
        "params" => params,
        "rbc" => rbc,
        "mu" => μ,
        "theta_final" => θ_globulos
    )
end


"""
Genera todas las gráficas y animación
"""
function generar_resultados(resultados::Dict)
    params = resultados["params"]
    energias = resultados["energias"]
    configs = resultados["configuraciones"]
    μ = resultados["mu"]
    H = params.H_magnitud
    
    # Separar fase de equilibrio y producción
    n_eq = div(params.N_equilibrio, params.save_interval)
    energias_eq = energias[1:n_eq]
    energias_prod = energias[n_eq+1:end]
    configs_prod = configs[n_eq+1:end, :]
    
    println("\n" * "=" ^ 70)
    println("ANÁLISIS DE RESULTADOS")
    println("=" ^ 70)
    
    # Estadísticas
    E_min_teorico = -μ * H * params.N_globulos
    E_max_teorico = μ * H * params.N_globulos
    E_promedio = mean(energias_prod)
    E_std = std(energias_prod)
    
    println("Energía mínima teórica: $(@sprintf("%.3e", E_min_teorico)) J")
    println("Energía máxima teórica: $(@sprintf("%.3e", E_max_teorico)) J")
    println("Energía promedio (producción): $(@sprintf("%.3e", E_promedio)) J")
    println("Desviación estándar: $(@sprintf("%.3e", E_std)) J")
    println("Fracción del mínimo: $(@sprintf("%.2f", E_promedio/E_min_teorico))")
    
    # Análisis de configuraciones finales
    θ_final = vec(configs_prod[end, :])
    θ_promedio = mean(θ_final)
    println("\nÁngulo promedio final: $(@sprintf("%.2f", rad2deg(θ_promedio)))°")
    println("Glóbulos con θ < 45°: $(sum(θ_final .< π/4)) / $(params.N_globulos)")
    println("Glóbulos con θ > 135°: $(sum(θ_final .> 3π/4)) / $(params.N_globulos)")
    
    # ========== GRÁFICA 1: Curva teórica E(θ) ==========
    θ_range = range(0, π, length=1000)
    E_teorica = [energia_magnetica(μ, H, θ) for θ in θ_range]
    
    p1 = plot(rad2deg.(θ_range), E_teorica .* 1e18,
              xlabel="Ángulo θ (grados)", ylabel="Energía (×10⁻¹⁸ J)",
              title="Energía vs Orientación (1 glóbulo)",
              linewidth=3, color=:blue, legend=false,
              grid=true, size=(800, 500))
    scatter!([0, 180], [E_teorica[1], E_teorica[end]] .* 1e18,
             markersize=8, color=:red, label="Extremos")
    
    # ========== GRÁFICA 2: Serie temporal completa ==========
    pasos = (1:length(energias)) .* params.save_interval
    
    p2 = plot(pasos, energias .* 1e18,
              xlabel="Paso MC", ylabel="Energía Total (×10⁻¹⁸ J)",
              title="Evolución Temporal (N=$(params.N_globulos) glóbulos)",
              linewidth=2, color=:darkblue, legend=false,
              grid=true, size=(800, 500))
    vline!([params.N_equilibrio], linewidth=2, color=:red,
           linestyle=:dash, label="Fin equilibrio")
    
    # ========== GRÁFICA 3: Histograma de orientaciones ==========
    p3 = histogram(rad2deg.(θ_final), bins=30,
                   xlabel="Ángulo θ (grados)", ylabel="Frecuencia",
                   title="Distribución de Orientaciones (Estado Final)",
                   color=:orange, alpha=0.7, legend=false,
                   grid=true, size=(800, 500))
    
    # ========== GRÁFICA 4: Energía promedio vs paso (producción) ==========
    pasos_prod = ((n_eq+1):length(energias)) .* params.save_interval
    
    p4 = plot(pasos_prod, energias_prod .* 1e18,
              xlabel="Paso MC", ylabel="Energía Total (×10⁻¹⁸ J)",
              title="Fase de Producción",
              linewidth=2, color=:green, legend=false,
              grid=true, size=(800, 500))
    hline!([E_promedio * 1e18], linewidth=2, color=:red,
           linestyle=:dash, label="Promedio")
    
    # Guardar gráficas
    println("\n" * "=" ^ 70)
    println("Guardando resultados...")
    println("=" ^ 70)
    
    savefig(p1, "energia_vs_angulo.png")
    println("Guardado: energia_vs_angulo.png")
    
    savefig(p2, "serie_temporal_completa.png")
    println("Guardado: serie_temporal_completa.png")
    
    savefig(p3, "histograma_orientaciones.png")
    println("Guardado: histograma_orientaciones.png")
    
    savefig(p4, "fase_produccion.png")
    println("Guardado: fase_produccion.png")
    
    # Panel combinado
    p_combined = plot(p1, p2, p3, p4, layout=(2,2), size=(1600, 1200))
    savefig(p_combined, "resultados_completos.png")
    println("Guardado: resultados_completos.png")
    
    # Guardar datos numéricos
    df_energias = DataFrame(
        paso = pasos,
        energia = energias,
        fase = [i <= n_eq ? "equilibrio" : "produccion" for i in 1:length(energias)]
    )
    CSV.write("energias.csv", df_energias)
    println("Guardado: energias.csv")
    
    df_config_final = DataFrame(
        globulo_id = 1:params.N_globulos,
        theta_rad = θ_final,
        theta_deg = rad2deg.(θ_final),
        energia = [energia_magnetica(μ, H, θ) for θ in θ_final]
    )
    CSV.write("configuracion_final.csv", df_config_final)
    println("Guardado: configuracion_final.csv")
    
    # ANIMACIÓN 
    println("\nGenerando animación...")
    crear_animacion(resultados)
    
    println("\n" * "=" ^ 70)
    println("¡ANÁLISIS COMPLETO!")
    println("=" ^ 70)
end

"""
Crea animación de la evolución del sistema
"""
function crear_animacion(resultados::Dict)
    configs = resultados["configuraciones"]
    params = resultados["params"]
    μ = resultados["mu"]
    H = params.H_magnitud
    
    # Seleccionar frames (una muestra cada 10 configuraciones guardadas)
    n_frames = min(200, size(configs, 1))
    indices_frames = Int.(round.(range(1, size(configs, 1), length=n_frames)))
    
    anim = @animate for i in indices_frames
        θ_config = configs[i, :]
        E_total = sum([energia_magnetica(μ, H, θ) for θ in θ_config])
        paso = i * params.save_interval
        
        # Representación polar de las orientaciones
        scatter(rad2deg.(θ_config), ones(params.N_globulos),
                xlabel="Ángulo θ (grados)", ylabel="",
                title="Paso $paso | E = $(@sprintf("%.2e", E_total)) J",
                xlims=(0, 180), ylims=(0.5, 1.5),
                markersize=3, alpha=0.6, color=:blue,
                legend=false, yticks=[], size=(1000, 400))
        vline!([90], linewidth=2, color=:red, linestyle=:dash, alpha=0.5)
    end
    
    mp4(anim, "evolucion_montecarlo.mp4", fps=15)
    println("Guardado: evolucion_montecarlo.mp4")
end


function main()
    println("\n")
    println("╔" * "═"^68 * "╗")
    println("║" * " "^10 * "SIMULACIÓN MONTE CARLO: GLÓBULOS ROJOS" * " "^20 * "║")
    println("║" * " "^15 * "Método de Metropolis - Ensamble Canónico" * " "^13 * "║")
    println("╚" * "═"^68 * "╝")
    
    # Configurar parámetros
    # NOTA: Para ver dinámica real con glóbulos rojos diamagnéticos,
    # necesitamos campo MUCHO más débil o temperatura MUCHO más alta
    params = SimulationParams(
        1000,          # N_globulos
        50_000,        # N_equilibrio
        200_000,       # N_produccion
        0.001,         # H_magnitud [Tesla] - Campo muy débil
        1000.0,        # T [K] - Temperatura elevada para ver fluctuaciones
        100            # save_interval
    )
    
    # Ejecutar simulación
    resultados = simular_ensemble(params)
    
    # Analizar y visualizar
    generar_resultados(resultados)
    
    println("\n¡Simulación finalizada con éxito!")
    println("Revisa los archivos generados:")
    println("  - resultados_completos.png")
    println("  - evolucion_montecarlo.mp4")
    println("  - energias.csv")
    println("  - configuracion_final.csv")
end

# Ejecutar
main()