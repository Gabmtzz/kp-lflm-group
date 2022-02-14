module KPpack
# ===============================================================================================================================
# libraries
using LinearAlgebra
using CSV
using DataFrames
using DelimitedFiles
using SparseArrays, BlockBandedMatrices
using Arpack

#export Materials, parMat, DOS, DiagM, Plotbands, PlotDOS,

#materials parameters
mutable struct Materials
    material #material or alloy name with composition: alloy_comp
    g1  #Luttinger parameter gamma_1
    g2  #Luttinger parameter gamma_
    g3  #Luttinger parameter gamma_1
    Eg  #material's band gap
    Ep  # Ep parameter
    F   # F parameter
    k     # kappa Luttinger parameter
    delta # spin-orbit valence band split
    VBO   #valence band offset
end

# =================================================================================
# Bowing parameters
mutable struct BowPar
    alloy
    cEg1
    cEg2
    cEp
    cF
    cDelta
    cVBO
end
# ====================================================================================
# using for the formation of quantum wells
struct mat
    material # name of alloy or material in the same format that structure "Materials"
    size     # size of the layer
end
# ===========================================================================================
#diagonalizes the Hamiltonian matrix and find the DOS
include("solver.jl")
#Creates the matrix of Kane with/without spin-orbit coupling
include("Hamiltonian.jl")
#adquires the material parameters
include("Params.jl")
#Shows the band structure and the DOS
include("plot.jl")
#finite diffrerences kp 8 band hamiltonian matrix
include("FinitediffQW.jl")

end # module
