using Documenter
using Distances
using TDARipserer

if VERSION ≥ v"1.8-DEV" || VERSION < v"1.7-DEV" || !Sys.islinux()
    @warn "Doctests only run on Linux and Julia 1.7"
else
    DocMeta.setdocmeta!(
        TDARipserer, :DocTestSetup, :(using TDARipserer; using Distances); recursive=true
    )
    doctest(TDARipserer)
end
