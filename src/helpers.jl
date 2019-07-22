"""
    is_missing(MetData, "Date")
Find if a column is missing from a DataFrame.

# Arguments
- `data::DataFrame`: a DataFrame
- `column::String`: a column name

# Return
A boolean: `true` if the column is missing, `false` if it is present.

# Examples
```julia
df= DataFrame(A = 1:10)
is_missing(df,"A")
false
is_missing(df,"B")
true
```
"""
function is_missing(data::DataFrame,column::String)::Bool
  columns= names(data)
  for i in 1:length(columns)
    is_in_col= columns[i] == Symbol(column)
    if is_in_col
      return false
    end
  end
  return true
end


"""
    is_missing(Dict("test"=> 2), "test")
Find if a key is missing from a dictionary.

# Arguments
- `data::Dict`: a dictionary
- `column::String`: a key (parameter) name

# Return
A boolean: `true` if the key is missing, `false` if it is present.

# Examples
```julia
Parameters= Dict("Stocking_Coffee"=> 5580)
is_missing(Parameters,"Stocking_Coffee")
false
is_missing(Parameters,"B")
true
```
"""
function is_missing(data::Dict,column::String)::Bool
  try
    data[column]
  catch error
    if isa(error, KeyError)
      return true
    end
  end
  return false
end



"""
    warn_var("Date","Start_Date from Parameters","warn")
Warn or stop execution if mandatory meteorology input variables are not provided.
It helps the user to know which variable is missing and/or if there are replacements

# Arguments
- `Var::String`: Input variable name
- `replacement::String`: Replacement variable that is used to compute `"Var"`
- `type::String`: Type of error to return : either

# Note
* This function helps to debug the model when some mandatory meteorological variables
are missing from input: either an error (default), or a warning.
* If the `"replacement"` variable is not provided in the meteorology file either, this function
will return an error with a hint on which variables can be provided to compute `"Var"`

# Examples
```julia
warn_var("Date","Start_Date from Parameters","warn")
```
"""
function warn_var(Var::String,replacement::String,type::String="error")
  if type=="error"
    error(string("$Var missing from input Meteo. Cannot proceed unless provided.",
                 " Hint: $Var can be computed alternatively using $replacement if provided in Meteo file")
               )
  else
    println("$Var missing from input Meteo. Computed from $replacement")
  end
end

"""
    warn_var("Date")
Stop execution if mandatory meteorology input variable is not provided.

# Arguments
- `Var::String`: Input variable name

"""
function warn_var(Var::String)
  error("$Var missing from input Meteo. Cannot proceed unless provided.")
end

function cos°(x::Float64)::Float64
  cos(x*π/180.0)
end

function sin°(x::Float64)::Float64
  sin(x*π/180.0)
end

function tan°(x::Float64)::Float64
  tan(x*π/180.0)
end

function acos°(x::Float64)::Float64
  acos(x)*180.0/π
end

function asin°(x::Float64)::Float64
  asin(x)*180.0/π
end


function atan°(x::Float64)::Float64
  atan(x)*180.0/π
end


"""
Trigonometric Functions (degree)

These functions give the obvious trigonometric functions. They respectively compute the cosine, sine, tangent,
arc-cosine, arc-sine, arc-tangent with input and output in degree.

# Returns
The output in degree

# Details
The conversions between radian to degree is: ``x\\cdot \\frac{\\pi,180}``


# Examples
```julia
# cosinus of an angle of 120 degree:
cos°(120)
# should yield -0.5, as in the base version:
cos(120*π/180)
```
"""
cos°,sin°,tan°,acos°,asin°,atan°





function struct_to_tuple(structure,instance)
  structure_names= fieldnames(structure)
  eval_names= map(x -> :($(instance).$x), structure_names)
  # (; zip(structure_names, structure_values)...)
  NamedTuple{structure_names}(map(eval,eval_names))
end