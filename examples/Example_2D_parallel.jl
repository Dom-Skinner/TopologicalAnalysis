using TopologicalAnalysis
using Random
using Plots

# Here we demonstrate how to do paralellism in 2D and test timings. 
# In Example_3D there is an example of how to do this for the 3D case.


# Reuse the example data from Example_2D
function generate_random_pts()
    # generates a MotifArray for a random distribution of points
    pos = rand(1000,2)
    delaunay_inf = find_delaunay_network(pos, periodic=false, alpha = 0)
    return compute_motifs(delaunay_inf)
end

function generate_perturbed_grid()
    # generates a MotifArray for a perturbed grid of points
    pos = [[x,y] for x in range(0,1,31), y in range(0,1,31)][:]
    pos .+= [0.01*randn(2) for i = 1:length(pos)]
    delaunay_inf = find_delaunay_network(pos, periodic=false, alpha = 0)
    return compute_motifs(delaunay_inf)
end

# Take 5 realizations of the random points and 5 of the perturbed grid in one total array
total_motif_array = vcat([generate_random_pts() for i = 1:5], [generate_perturbed_grid() for i = 1:5])

# Compute the flip graph for these in serial
@time flip_graph_combined = compute_flip(total_motif_array...; restrict = 0,thresh=0.5)

# Now in parallel. First add the processes and load TopologicalAnalysis on all of them
using Distributed
addprocs(4)
@everywhere using TopologicalAnalysis

# Redo the timing
@time flip_graph_combined = compute_flip(total_motif_array...; restrict = 0,thresh=0.5)

# Go back to one process
rmprocs(workers())

# Compute distance matrices for diffusion and optimal transport distances in serial
@time d_diff = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=false)
@time d_OT = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=true)

# Now in parallel.
# We don't expect much of a speed up from d_diff from this data set as it was already fast, but d_OT should see a significant speed up
addprocs(4)
@everywhere using TopologicalAnalysis

@time d_diff = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=false)
@time d_OT = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=true)

# Check timing as you scale to more processes to make sure it's working!