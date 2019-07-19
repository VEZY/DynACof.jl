
function import_parameters(path,Names)
    if path == "package"
        paths= repeat(["package"],length(Names))
        paths_names= keys(Names)
        paths = NamedTuple{paths_names}(paths)       
    else
        paths= map(x -> normpath(string(path,"/",x)),Names)
    end
    if Names.tree==""
        Parameters= (constants= read_param_file(:constants,paths.constants), site= read_param_file(:site, paths.site),
                     soil= read_param_file(:soil, paths.soil), coffee= read_param_file(:coffee,paths.coffee), tree= ())
    else
        Parameters= (constants= read_param_file(:constants,paths.constants), site= read_param_file(:site, paths.site),
                     soil= read_param_file(:soil, paths.soil), coffee= read_param_file(:coffee,paths.coffee), 
                     tree= read_param_file(:tree,paths.tree))
    end
    Parameters
end


"""
    read_param_file(structure::Symbol,filepath::String="package")

Read DynACof parameter files and create the structure according to its structure.
If parameters are missing from the file, the structure is filled with the default values.

# Arguments
- `structure::Int64`: The structure type. Must be one of `constants`, `site`, `soil`, `coffee`, `tree`
- `filepath::Float64`: The path to the parameter file

# Return
The corresponding structure with the values read from the parameter file.

# Examples
```julia
julia> read_param_file(:constants)
constants(0.0010130000000000007, 0.622, 101.325, 0.5, 9.81, 287.0586, 8.314, 273.15, 0.41, 1.0000000000000006e-6, 1367.0, 5.670367e-8, 0.018, 4.57, 2.45, 0.4,2.15e-5)

julia> read_param_file(:site)
DynACof.site("Aquiares", "1979/01/01", 9.93833, -83.72861, 6, 1040.0, 25.0, 0.58, 0.144)
```
"""
function read_param_file(structure::Symbol,filepath::String="package")
    if basename(filepath)=="package"
        return eval(:($structure()))
    else
        params= evalfile(filepath)
        b= map((x,y) -> :($x = $y),collect(keys(params)),collect(values(params)))
        return eval(:($structure(;$(b...))))
    end
end
