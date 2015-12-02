push!(LOAD_PATH, pwd());
import fmmv

n=10000

sources = rand(2, n)
charges = rand(n)

pot = fmmv.fmmv2d(sources, charges)       
