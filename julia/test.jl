push!(LOAD_PATH, pwd());
import fmmv

n=10000

sources = rand(3, n)
charges = rand(n)

pot = fmmv.fmmv3d(sources, charges)       
