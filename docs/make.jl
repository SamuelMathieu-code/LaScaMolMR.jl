push!(LOAD_PATH, "../src/")
using Documenter, LaScaMolMR


makedocs(sitename="LaScaMolMR.jl", format = Documenter.HTML(prettyurls = false))