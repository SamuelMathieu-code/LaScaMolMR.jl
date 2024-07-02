push!(LOAD_PATH, "../src/")
# import Pkg; Pkg.add("Documenter")
using Documenter, LaScaMolMR

makedocs(sitename="LaScaMolMR.jl", format = Documenter.HTML(prettyurls = false),  authors = "Samuel Mathieu")
