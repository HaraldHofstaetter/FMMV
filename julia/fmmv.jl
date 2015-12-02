module fmmv 

export fmmv2d #, fmmv2d_initialize, fmmv2d_evaluate, fmmv2d_finalize
export fmmv_complex2d #, fmmv_complex2d_initialize, fmmv_complex2d_evaluate, fmmv_complex2d_finalize
export fmmv3d #, fmmv3d_initialize, fmmv3d_evaluate, fmmv3d_finalize
export FmmvOptions, fmmvGetDefaultOptions2d, fmmvGetDefaultOptions3d

type FmmvOptions 
	beta::Cdouble # default: 0 

	pM::Cint # default:6 
	pL::Cint # default:6 
	s::Cint 

	ws::Cint #  default: 1 
	reducedScheme::Cint

	scale::Cdouble # default: 1 
	splitThreshold::Cint  
	splitTargetThreshold::Cint
	levels::Cint
	directEvalThreshold::Cint
	periodicBoundaryConditions::Cint # default: 0 
	extrinsicCorrection::Cint # default: 0 
	useHilbertOrder::Cint # default: 0, i.e. use Molton order 
	directEvalAccuracy::Cint # default: 2 (double) resp. 1 (single) 
	useFarfieldNearfieldThreads::Cint # default: 0 
	
	PAPIeventSet::Cint
end

function Base.show(io::IO, opt::FmmvOptions)
    for x in fieldnames(FmmvOptions)
        println(io, x, ": ", opt.(x))
    end
end    


function fmmvGetDefaultOptions2d(T::Type{Float32}) 
     ccall((:fmmvGetDefaultOptions, :libfmmv2df), FmmvOptions, () )
end        

function fmmvGetDefaultOptions2d(T::Type{Float64}) 
      ccall((:fmmvGetDefaultOptions, :libfmmv2d), FmmvOptions, () )
end   

function fmmvGetDefaultOptions3d(T::Type{Float32}) 
     ccall((:fmmvGetDefaultOptions, :libfmmv3df), FmmvOptions, () )
end        

function fmmvGetDefaultOptions3d(T::Type{Float64}) 
     ccall((:fmmvGetDefaultOptions, :libfmmv3d), FmmvOptions, () )
end    





function fmmv2d{T<:Union{Float32,Float64}}(sources::Array{T, 2},
                   charges::Array{T, 1};
                   dipoleMoments::Array{T, 2} = Array(T,0,0),
                   targets::Array{T, 2} = Array(T,0,0),
                   computeGradients::Bool = false,
                   options::FmmvOptions = fmmvGetDefaultOptions2d(T),
                   getStatistics::Bool = false)
    dim = 2               
    if size(sources)[1]!=dim
         throw(DimensionMismatch("expected size(sources) = ($dim,n), got $(size(sources))"))
    end
    nSources = size(sources)[2]
    if length(charges) != nSources
        throw(DimensionMismatch("expected length(charges) = $(nSources), got $(length(charges))"))
    end
    empty = Array(T,0,0)
    withDipoleMoments = ( dipoleMoments!=empty )
    if withDipoleMoments
        if size(dipoleMoments) != (dim, nSources)
             throw(DimensionMismatch("expected size(dipoleMoments) = ($dim,$(nSources)), got $(size(dipoleMoments))"))
        end
    end    
    withTargets = ( targets!=empty )
    nTargets = nSources
    if withTargets
        if size(targets)[1] != dim
             throw(DimensionMismatch("expected size(targets) = ($dim,n), got $(size(targets))"))
        end 
        nTargets = size(targets)[2]
    end

    pot = Array(T, nTargets)
    grad = empty
    if computeGradients
        grad = Array(T, dim, nTargets)
    end

    err = "for error message"
    NULL = convert(Ptr{Void}, 0)

    if T==Float32
        ccall((:fmmv2df, :libfmmv2df), 
            Void, (Cint, Ptr{T}, Ptr{T}, Ptr{T}, Cint, Ptr{T}, Ptr{T}, Ptr{T}, 
            Ptr{FmmvOptions}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL),
            &options, NULL, &err)
    else            
        ccall((:fmmv2d, :libfmmv2d), 
            Void, (Cint, Ptr{T}, Ptr{T}, Ptr{T}, Cint, Ptr{T}, Ptr{T}, Ptr{T},
            Ptr{FmmvOptions}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL), 
            &options, NULL, &err)
    end 

    if computeGradients
        return pot, grad
    else
        return pot
    end    
end    


function fmmv3d{T<:Union{Float32,Float64}}(sources::Array{T, 2},
                   charges::Array{T, 1};
                   dipoleMoments::Array{T, 2} = Array(T,0,0),
                   targets::Array{T, 2} = Array(T,0,0),
                   computeGradients::Bool = false,
                   options::FmmvOptions = fmmvGetDefaultOptions2d(T),
                   getStatistics::Bool = false)
    dim = 3               
    if size(sources)[1]!=dim
         throw(DimensionMismatch("expected size(sources) = ($dim,n), got $(size(sources))"))
    end
    nSources = size(sources)[2]
    if length(charges) != nSources
        throw(DimensionMismatch("expected length(charges) = $(nSources), got $(length(charges))"))
    end
    empty = Array(T,0,0)
    withDipoleMoments = ( dipoleMoments!=empty )
    if withDipoleMoments
        if size(dipoleMoments) != (dim, nSources)
             throw(DimensionMismatch("expected size(dipoleMoments) = ($dim,$(nSources)), got $(size(dipoleMoments))"))
        end
    end    
    withTargets = ( targets!=empty )
    nTargets = nSources
    if withTargets
        if size(targets)[1] != dim
             throw(DimensionMismatch("expected size(targets) = ($dim,n), got $(size(targets))"))
        end 
        nTargets = size(targets)[2]
    end

    pot = Array(T, nTargets)
    grad = empty
    if computeGradients
        grad = Array(T, dim, nTargets)
    end

    err = "for error message"
    NULL = convert(Ptr{Void}, 0)

    if T==Float32
        ccall((:fmmv3df, :libfmmv3df), 
            Void, (Cint, Ptr{T}, Ptr{T}, Ptr{T}, Cint, Ptr{T}, Ptr{T}, Ptr{T}, 
            Ptr{FmmvOptions}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL),
            &options, NULL, &err)
    else            
        ccall((:fmmv3d, :libfmmv3d), 
            Void, (Cint, Ptr{T}, Ptr{T}, Ptr{T}, Cint, Ptr{T}, Ptr{T}, Ptr{T},
            Ptr{FmmvOptions}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL), 
            &options, NULL, &err)
    end 

    if computeGradients
        return pot, grad
    else
        return pot
    end    
end    


function fmmv_complex2d{T<:Union{Float32,Float64}}(sources::Array{Complex{T}, 1};
                   charges::Array{Complex{T}, 1} = Array(Complex{T},0),
                   dipoleMoments::Array{Complex{T}, 1} = Array(Complex{T},0),
                   targets::Array{Complex{T}, 1} = Array(Complex{T},0),
                   computeGradients::Bool = false,
                   options::FmmvOptions = fmmvGetDefaultOptions2d(T),
                   getStatistics::Bool = false)
    nSources = length(sources)
    empty = Array(Complex{T},0)
    withCharges = ( charges!=empty )
    if withCharges
        if length(charges) != nSources
            throw(DimensionMismatch("expected length(charges) = $(nSources), got $(length(charges))"))
        end
    end    
    withDipoleMoments = ( dipoleMoments!=empty )
    if withDipoleMoments
        if length(dipoleMoments) != nSources
             throw(DimensionMismatch("expected length(dipoleMoments) = $(nSources), got $(length(dipoleMoments))"))
        end
    end    
    withTargets = ( targets!=empty )
    nTargets = nSources
    if withTargets
        nTargets = size(targets)[2]
    end

    pot = Array(Complex{T}, nTargets)
    grad = empty
    if computeGradients
        grad = Array(Complex{T}, nTargets)
    end

    err = "for error message"
    NULL = convert(Ptr{Void}, 0)

    if T==Float32
        ccall((:fmmv_complex2df, :libfmmv2df), 
            Void, (Cint, Ptr{Complex{T}}, Ptr{Complex{T}}, Ptr{Complex{T}}, Cint, Ptr{Complex{T}}, Ptr{Complex{T}}, Ptr{Complex{T}}, 
            Ptr{FmmvOptions}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, (withCharges?charges : NULL), (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL),
            &options, NULL, &err)
    else            
        ccall((:fmmv_complex2d, :libfmmv2d), 
            Void, (Cint, Ptr{Complex{T}}, Ptr{Complex{T}}, Ptr{Complex{T}}, Cint, Ptr{Complex{T}}, Ptr{Complex{T}}, Ptr{Complex{T}},
            Ptr{FmmvOptions}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, (withCharges?charges : NULL), (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL), 
            &options, NULL, &err)
    end 

    if computeGradients
        return pot, grad
    else
        return pot
    end    
end    



end
