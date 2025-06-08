
module Build

using JuliaTools

# These are the dependencies for Callisto
# which are not in the General registry nor the standard library.
function deps_for_callisto()
    b(x,y) = github_repo("bittide/$(x).jl", y)
    b(x) = github_repo("bittide/$(x).jl")
    local_package_repos = [
        b("JuliaTools"),
        b("LogData"),
        b("Piecewise"),
        b("Topology"),
        b("StateSystems"),
        b("BittideTools"),
    ]
    return local_package_repos
end

function main()
    working_dir = joinpath(pwd(), "working")
    package_dir = pwd()
    local_package_repos = deps_for_callisto()
    build_and_test_package(working_dir, package_dir, local_package_repos)
end


end


using .Build
Build.main()


