function get_files_dir(Data_dir)
      DIR = readdir(Data_dir)
      str_arr = DIR[[occursin("_avg.txt",d) for d in DIR]]
      str_arr = [String(SubString(d,1:(length(d)-10))) for d in str_arr].*"_"
      return str_arr
end

function initialize()

	## ---- For Drosophila -----
	Data_dir = "/home/gridsan/dskinner/Packing/Data/Drosophila/"
	network_save_file = Data_dir*"w_network_Droso"
	str_arr = ["WT_1frame_","WT_2frame_","WT_3frame_"]
	weight_arr_in = vcat([readin(Data_dir*str) for str in str_arr]...)
    N = Int(round(sum([length(w) for w in weight_arr_in])/length(weight_arr_in)))
	weight_start = Dict(amalg2(vcat([readin(Data_dir*str)[4:7] for str in str_arr]...)))
	weight_end = Dict(amalg2(vcat([readin(Data_dir*str)[end-3:end] for str in str_arr]...)))

	#=	
	## ---- For monodisperse ellipsoids -----
	network_save_file = Data_dir*"w_network"
	str_mat = ["Ells1/Spheres_","Ells1.5/Ells1.5_","Ells2/Ells2_","Ells2.5/Ells2.5_"]
	weight_arr_in = [subsample_dist(w,20) for w in weight_in]
	=#

	#=	
	## ----- For Active Brownian Particles ----
	Data_dir = "/home/gridsan/dskinner/Packing/Data/ABP/"
	str_arr = get_files_dir(Data_dir)
	network_save_file = Data_dir*"w_network_ABP"
	weight_arr_in = vcat([readin(Data_dir*str) for str in str_arr]...)
	=#


	#=
	## ----- For 3D Biofilm data----
	Data_dir = "/home/gridsan/dskinner/Packing/Data/BiofilmVib/"
	#Data_dir = "/home/gridsan/dskinner/Packing/Data/BiofilmSizes/"
	avg_files = filter(x->occursin("_avg.txt",x), readdir(Data_dir))
	network_save_file = Data_dir*"w_network"
	weight_arr_in = vcat([readin(Data_dir*str,0) for str in avg_files]...)
	=#

    	 #=
    	## ---- For swarming data
    	Data_dir = "/home/gridsan/dskinner/Packing/Data/Swarm/"
	avg_files = filter(x->occursin("_avg.txt",x), readdir(Data_dir))
	network_save_file = Data_dir*"w_network"
	weight_arr_in = vcat([readin(Data_dir*str,0) for str in avg_files]...)
	 =#
    
	#=
	## ----- For Polydisperse ellipsoids ----
	Data_dir = "/home/gridsan/dskinner/Packing/Data/"
	network_save_file = Data_dir*"Poly/w_network_Poly"
	str_arr = get_files_dir(Data_dir*"Poly/")
    	weight_poly = weight_in(Data_dir*"Poly/",str_arr)
	str_arr = get_files_dir(Data_dir*"Mono/")
    	weight_mono = weight_in(Data_dir*"Mono/",str_arr)
    	weight_arr_in = vcat(weight_poly,weight_mono)
	=#
    df = CSV.read("/home/gridsan/dskinner/Packing/TransportPlan.csv",delim=",";header=false)
    w_geo = [df[:,i] for i in 1:size(df,2)]
	g,weight = calculate_distance_matrix_parallel(network_save_file,weight_arr_in)
	w_naive = geodesic(weight_start,weight_end,1/9:1/9:8/9.,network_save_file,N)
	weight = vcat(weight,w_geo,w_naive)
	println("Length of array is: ", length(weight))
	return g,weight
end


function triangle_index(k)
	# Helps us to iterate through upper triangular part of matrix
	T = n ->Int(0.5*n*(n+1))
	kk = Int(ceil( sqrt(2*k + 0.25) - 0.5))
	return k - T(kk-1),kk + 1
end

function rem_self_edges!(g)
    # Removes self loops
    for e in collect(edges(g))
        if src(e) == dst(e)
            rem_edge!(g,e)
        end
    end
end
