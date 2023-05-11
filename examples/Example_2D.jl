# Example script showing how to use the package LocalCellularStructure

# First load the required packages. For instructions on how to install a custom package see, https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/11-developing-julia-packages
using LocalCellularStructure
using Random


save_directory = "./"

# First let us create some data with some data, say random points in 2D.
Random.seed!(1234)
Positions = rand(1000,2) # 1000 random points.

# First we wish to compute the Delaunay. The function below returns a structure containing
# all information relevant to the delaunay triangulation.
# We tell it that the data is not periodic, and also, by setting alpha=0 we use the default
# value of alpha. We can change this as needed.
delaunay_info = find_delaunay_network(Positions, periodic=false, alpha = 0)

# delaunay_info is a structure of (custom) type TopologicalNetwork.
# delaunay_info.simplices is a N×3 list of vertices making up the N triangles of the delauany.
# delaunay_info.not_edge are the vertices identified as interior points (can be useful to plot
#   to make sure the alpha value has been correctly chosen).
# There are other fields, which we can ignore for now.

# one can save a TopologicalNetwork type as a HDF5 file by calling
save(save_directory*"test_delaunay_network.h5",delaunay_info)

# and one can load it again by
delaunay_info_reloaded = load(save_directory*"test_delaunay_network.h5")

# Next, we compute the motifs for a network, using the function compute_motifs
motif_array = compute_motifs(delaunay_info)

# motif_array is of type MotifArray
# motif_array.idx contains the original vertex number, and motif_array.tvec 
# contains the corresponding weinberg vector.

# we can save a MotifArray, much as we can save a TopologicalNetwork
save(save_directory*"motif_array.h5",motif_array)

# and load it again by
motif_array_reloaded = load(save_directory*"motif_array.h5")

# a MotifArray is associated with a particular Delauany triangulation.
# suppose we want to combine multiple MotifArray's to get a motif distribution
# In this case call avg_motif which takes any number of arguments 
# (pretend for now motif_array and motif_array_reloaded are different)
motif_distribution = avg_motif(motif_array,motif_array_reloaded)
# This returns a type of MotifDist, where MotifDist.map is a dictionary 
# taking a weinberg vector to a frequency. MotifDist's can be saved in much the same way as before


# Now on to computing the flip graph. Pass in a MotifDist, or any number of MotifArray's, and
# it will take all observed motifs, and calculate whether any pair of them are connected with 
# an edge, resulting in a graph contained in a custom structure FlipGraph.
flip_graph = compute_flip(motif_distribution; restrict = 0,thresh=0.0)
# setting restrict = n ignores motifs that you have only seen n times or fewer; only to be used
# for large systems where the flip graph gets too big/too expensive to calculate etc.
# thresh sets a threshold to remove vertices from the flip graph that are not important 
# according to the page rank algorithm. The lower the threshold the more vertices are kept.

# flip graphs can be saved/loaded as before.
save(save_directory*"flip_graph.h5",flip_graph)


# Now we are ready to compute a distance between two MotifDist's or two MotifArray's
# First stack them into a single vector
stacked_motif_array = vcat(motif_array,motif_array_reloaded)
# Then call the distance calculation
d = calculate_distance_matrix(flip_graph, stacked_motif_array, optimal_transport=true)
# as expected d is a matrix of zeros, because motif_array and motif_array_reloaded are the same thing.


# Now let us do a more realistic use case.
# We will take 5 realizations of a random distribution of points, and 5 realizations of a perturbed grid
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

# take 5 realizations of the random points and 5 of the perturbed grid in one total array
total_motif_array = vcat([generate_random_pts() for i = 1:5], [generate_perturbed_grid() for i = 1:5])

# compute the flip graph for these. N.B. writing the ... after total_motif_array is like calling
# compute_flip(total_motif_array[1],total_motif_array[2],...,total_motif_array[10];restrict = 0,thresh=0.5)
# See "splat" operation at https://docs.julialang.org/en/v1/base/base/
flip_graph_combined = compute_flip(total_motif_array...; restrict = 0,thresh=0.5)

# Now compute the combined distance matrix
d = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=true)

# To save, create a record of which inputs were random points and which were perturbed grids;
total_string_record = vcat(["Random_point_"*string(i) for i = 1:5],["Perturbed_grid_"*string(i) for i = 1:5])
using CSV, DataFrames
CSV.write(save_directory*"result.csv",DataFrame(d,total_string_record))

# We can further analyze this distance matrix using other Julia packages. 
# For instance, we might wish to compute the MDS embedding and plot this.
using Plots, MultivariateStats

# Pass the distnace matrix into a MDS algorithm (see MultivariateStats package)
MDS_coords = MultivariateStats.transform(MultivariateStats.fit(MDS,d, maxoutdim=2, distances=true))

p = scatter(MDS_coords[1,1:5],MDS_coords[2,1:5],label="Random points",
    xlabel="MDS PC 1", ylabel="MDS PC 2",grid=false,aspect_ratio=true)
scatter!(p,MDS_coords[1,6:10],MDS_coords[2,6:10],label="Perturbed grid")
