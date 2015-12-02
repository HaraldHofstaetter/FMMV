module fmmv 

export fmmv2d #, fmmv2d_initialize, fmmv2d_evaluate, fmmv2d_finalize
#export fmmv_complex2d, fmmv_complex2d_initialize, fmmv_complex2d_evaluate, fmmv_complex2d_finalize
export fmmv3d #, fmmv3d_initialize, fmmv3d_evaluate, fmmv3d_finalize


function fmmv2d{T<:Union{Float32,Float64}}(sources::Array{T, 2},
                   charges::Array{T, 1};
                   dipoleMoments::Array{T, 2} = Array(T,0,0),
                   targets::Array{T, 2} = Array(T,0,0),
                   computeGradients::Bool = false,
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
            Ptr{Void}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL),
            NULL, NULL, &err)
    else            
        ccall((:fmmv2d, :libfmmv2d), 
            Void, (Cint, Ptr{T}, Ptr{T}, Ptr{T}, Cint, Ptr{T}, Ptr{T}, Ptr{T},
            Ptr{Void}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL), 
            NULL, NULL, &err)
    end 


    if computeGradients
        return pot, grad
    else
        return pot
    end    
end    


function fmmv3d{T<:Union{Float32,Float64}}(sources::Array{T, 3},
                   charges::Array{T, 1};
                   dipoleMoments::Array{T, 3} = Array(T,0,0),
                   targets::Array{T, 3} = Array(T,0,0),
                   computeGradients::Bool = false,
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
            Ptr{Void}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL),
            NULL, NULL, &err)
    else            
        ccall((:fmmv3d, :libfmmv3d), 
            Void, (Cint, Ptr{T}, Ptr{T}, Ptr{T}, Cint, Ptr{T}, Ptr{T}, Ptr{T},
            Ptr{Void}, Ptr{Void}, Ptr{Ptr{Cchar}} ),
            nSources, sources, charges, (withDipoleMoments?dipoleMoments : NULL),
            nTargets, (withTargets?targets : NULL), pot, (computeGradients?grad : NULL), 
            NULL, NULL, &err)
    end 


    if computeGradients
        return pot, grad
    else
        return pot
    end    
end    


end
