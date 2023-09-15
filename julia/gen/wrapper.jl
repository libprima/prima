# Script to parse PRIMA headers and generate Julia wrappers.
using PRIMA_jll
using Clang
using Clang.Generators
using JuliaFormatter

function main()

  cd(@__DIR__)
  include = joinpath(PRIMA_jll.artifact_dir, "include", "prima")
  headers = [joinpath(include, "prima.h")]

  options = load_options(joinpath(@__DIR__, "prima.toml"))
  options["general"]["output_file_path"] = joinpath("..", "interface", "wrappers.jl")

  args = get_default_args()
  push!(args, "-I$include")
  
  ctx = create_context(headers, args, options)
  build!(ctx)

  path = options["general"]["output_file_path"]
  format_file(path, YASStyle())
  return nothing
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
