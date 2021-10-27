using CSV

include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

function write_raw_data_Ells(Data_dir)
    for n in 1:10
        dat_in = CSV.read(Data_dir*"PackLSD.2D"*string(n)*".dat";delim=' ',
        ignorerepeated=true,skipto=7,silencewarnings=true)
        Positions = [[dat_in[k,1], dat_in[k,2]] for k in 1:size(dat_in,1)]
        weinberg2D_wrap(Positions, Data_dir*"Ells_"*string(n),true)
    end
end

Data_dir = "/Users/Dominic/Documents/2d Cells/Data/"

Names = Data_dir.*["Ells1/","Ells1.5/","Ells2/","Ells2.5/","Ells3/","Ells3.5/","Ells4/"]

for name in Names
    write_raw_data_PV(name)
end
