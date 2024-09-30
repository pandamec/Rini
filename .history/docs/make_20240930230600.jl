using Documenter
using Rini

makedocs(
    sitename = "Rini",
    format = Documenter.HTML(),
    modules = [Rini]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
