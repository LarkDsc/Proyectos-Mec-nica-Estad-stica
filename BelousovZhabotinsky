using Statistics
using Plots
using Random

# Dimensiones de la imagen
nx, ny = 600, 450

# Parámetros de la reacción
alpha, beta, gamma = 1.0, 1.0, 1.0

function convolve2d_periodic(arr, kernel)
    """Convolución 2D con condiciones de frontera periódicas"""
    ny, nx = size(arr)
    ky, kx = size(kernel)
    result = zeros(ny, nx)
    
    half_ky, half_kx = ky ÷ 2, kx ÷ 2
    
    for i in 1:ny
        for j in 1:nx
            sum_val = 0.0
            for ki in 1:ky
                for kj in 1:kx
                    # Índices con condiciones periódicas
                    pi = mod1(i + ki - half_ky - 1, ny)
                    pj = mod1(j + kj - half_kx - 1, nx)
                    sum_val += arr[pi, pj] * kernel[ki, kj]
                end
            end
            result[i, j] = sum_val
        end
    end
    
    return result
end

function update!(arr_current, arr_next)
    """Actualiza las concentraciones evolucionando en el tiempo"""
    
    # Kernel de promediado 3x3
    m = ones(3, 3) / 9.0
    
    # Array para almacenar los promedios locales
    s = zeros(3, ny, nx)
    
    # Calcular el promedio de cada especie en las 9 celdas vecinas
    for k in 1:3
        s[k, :, :] = convolve2d_periodic(arr_current[k, :, :], m)
    end
    
    # Aplicar las ecuaciones de reacción
    arr_next[1, :, :] = s[1, :, :] .+ s[1, :, :] .* (alpha .* s[2, :, :] .- gamma .* s[3, :, :])
    arr_next[2, :, :] = s[2, :, :] .+ s[2, :, :] .* (beta .* s[3, :, :] .- alpha .* s[1, :, :])
    arr_next[3, :, :] = s[3, :, :] .+ s[3, :, :] .* (gamma .* s[1, :, :] .- beta .* s[2, :, :])
    
    # Asegurar que las concentraciones estén en [0,1]
    clamp!(arr_next, 0.0, 1.0)
    
    return nothing
end

# Inicializar los arrays con cantidades aleatorias de A, B y C
Random.seed!(42)
arr_current = rand(3, ny, nx)
arr_next = zeros(3, ny, nx)

# Número de iteraciones
n_frames = 200

# Para la Figura 3: Series temporales de concentraciones promedio
mean_A = Float64[]
mean_B = Float64[]
mean_C = Float64[]

# Para la Figura 4: Distribuciones en frames específicos
# Guardaremos las distribuciones en frame 1, 50, 100, 150, 200
frames_to_save = [1, 50, 100, 150, 200]
distributions = Dict{Int, Dict{String, Vector{Float64}}}()

println("Ejecutando simulación y recolectando datos...")

# Crear animación
anim = @animate for i in 1:n_frames
    # Alternar entre arr_current y arr_next
    if i % 2 == 1
        update!(arr_current, arr_next)
        current_arr = arr_next
    else
        update!(arr_next, arr_current)
        current_arr = arr_current
    end
    
    # Guardar concentraciones promedio (Figura 3)
    push!(mean_A, mean(current_arr[1, :, :]))
    push!(mean_B, mean(current_arr[2, :, :]))
    push!(mean_C, mean(current_arr[3, :, :]))
    
    # Guardar distribuciones completas en frames específicos (Figura 4)
    if i in frames_to_save
        distributions[i] = Dict(
            "A" => vec(current_arr[1, :, :]),
            "B" => vec(current_arr[2, :, :]),
            "C" => vec(current_arr[3, :, :])
        )
    end
    
    # Visualización
    heatmap(current_arr[1, :, :], 
            c=:winter, 
            axis=nothing, 
            border=:none,
            size=(800, 600),
            aspect_ratio=:equal,
            title="BZ Reaction - Frame $i")
end

# Guardar la animación como MP4
println("Guardando animación...")
mp4(anim, "bz_reaction.mp4", fps=30)

println("Generando Figura 3: Series temporales...")

p1 = plot(1:n_frames, mean_A, 
          label="[A]", 
          linewidth=2,
          xlabel="Iteración temporal",
          ylabel="Concentración promedio",
          title="Evolución temporal de concentraciones promedio",
          legend=:topright,
          size=(800, 500),
          color=:blue)

plot!(p1, 1:n_frames, mean_B, 
      label="[B]", 
      linewidth=2,
      color=:red)

plot!(p1, 1:n_frames, mean_C, 
      label="[C]", 
      linewidth=2,
      color=:green)

savefig(p1, "figura3_series_temporales.png")
println("Figura 3 guardada como 'figura3_series_temporales.png'")

println("Generando Figura 4: Histogramas...")

# Crear subplots para diferentes tiempos
plots_array = []

for (idx, frame) in enumerate(sort(collect(keys(distributions))))
    data_A = distributions[frame]["A"]
    
    p = histogram(data_A,
                  bins=50,
                  xlabel="Concentración [A]",
                  ylabel="Frecuencia",
                  title="Frame $frame",
                  legend=false,
                  color=:skyblue,
                  alpha=0.7,
                  xlims=(0, 1))
    
    push!(plots_array, p)
end

# Combinar en una figura con múltiples subplots
p_combined = plot(plots_array..., 
                  layout=(2, 3),
                  size=(1200, 800),
                  plot_title="Evolución de la distribución de concentraciones")

savefig(p_combined, "figura4_histogramas.png")
println("Figura 4 guardada como 'figura4_histogramas.png'")

println("\n========== ESTADÍSTICAS FINALES ==========")
println("Concentración promedio final:")
println("  [A] = $(round(mean_A[end], digits=4))")
println("  [B] = $(round(mean_B[end], digits=4))")
println("  [C] = $(round(mean_C[end], digits=4))")

println("\nDesviación estándar final:")
final_arr = n_frames % 2 == 0 ? arr_current : arr_next
println("  σ_A = $(round(std(final_arr[1, :, :]), digits=4))")
println("  σ_B = $(round(std(final_arr[2, :, :]), digits=4))")
println("  σ_C = $(round(std(final_arr[3, :, :]), digits=4))")

println("\n¡Análisis completo!")
println("Archivos generados:")
println("  - bz_reaction.mp4")
println("  - figura3_series_temporales.png")
println("  - figura4_histogramas.png")
