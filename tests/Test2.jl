using Revise
using LocalCellularStructure
using Random

save_dir = homedir()*"/Documents/TopologicalAnalysis/data/"
Random.seed!(1234)
Positions = [randn(3) for k in 1:1000]

tg = find_delaunay_network(Positions, periodic=false, alpha = 0,
                                tol=0.6, edge_keep=false)

save(save_dir*"tmp_tn.h5",tg)
tg2 = load(save_dir*"tmp_tn.h5")

ma = compute_motifs(tg2)
save(save_dir*"tmp_ma.h5",ma)
ma2 = load(save_dir*"tmp_ma.h5")



fg = compute_flip(ma2; restrict = 0, edge_keep = false,thresh=1.5)
save(save_dir*"tmp_fg.h5",fg)
fg2 = load(save_dir*"tmp_fg.h5")
connected_flip_graph(fg)
