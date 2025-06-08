
module Build

using Pkg
using Pkg.TOML

mutable struct Repo
    url::String
    branch::String
end

github_repo(reponame, branch) = Repo("https://github.com/$(reponame).git", branch)
github_repo(reponame) = github_repo(reponame, "main")
match_repo_name(s) = match(r"(/)([a-zA-Z0-9_-]+)(\.jl)", s)
get_package_name(r::Repo) = extract_package_name_from_path(r.url)

"""
    extract_package_name_from_path(package_path)

Extracts the package name from a package path string.
Assumes the path ends with "/PackageName.jl".
Returns the "PackageName" or an empty string if the pattern is not matched.
"""
function extract_package_name_from_path(package_path)
    m = match_repo_name(package_path)
    if isnothing(m)
        return ""
    end
    return m[2]
end


"""
    make_package_list(project_toml_file::String) -> Vector{String}

Parses a `Project.toml` file and returns a list of package names
listed under its `[deps]` section.
"""
function get_deps_from_toml(project_toml_file)
    parsed_toml = TOML.parsefile(project_toml_file)
    deps = get(parsed_toml, "deps", Dict{String, Any}())
    return collect(keys(deps))
end


#
# Pkg.add will clone the repo, using the correct branch
# to somewhere in ~/.julia/packages, and then add the package
# to the active environment
#
add_package(r::Repo) = Pkg.add(PackageSpec(url = r.url, rev = r.branch))
dev_package(package_path) = Pkg.develop(PackageSpec(path = package_path))

"""
    build_and_test_package(working_dir, package_dir, local_package_repos)

Sets up an environment, installs dependencies, and tests a given Julia package.

Example:
  package_dir = "/home/somewhere/Callisto.jl")
  working_dir = "/home/emptydir"

The package_dir is not changed, but is used to find the Project.toml
and install the required dependencies.  Then the package is dev'ed at
that directory and tested. Dependencies and a Julia environment is stored
in the working_dir, which should be empty. The working_dir will be created
if it does not exist. The local_package_repos list must be in order, so
that they can be installed in the given order.

Steps:
1. Saves the current active Pkg project.
2. Creates `environment` and `repos` subdirectories within `working_dir`.
3. Parses the `Project.toml` of the `package_dir` to get required dependencies.
4. Separates dependencies into local and non-local.
5. Activates the new environment in `working_dir/environment`.
6. Installs local dependencies. Each one is cloned into `working_dir/repos, and installed using Pkg.dev
7. Installs non-local dependencies using Pkg.add
8. Installs the main package from `package_dir` in development mode.
9. Runs `Pkg.test()` for the main package.
10. Restores the original active Pkg project.
"""
function build_and_test_package(working_dir, package_dir, local_package_repos)
    save_project = dirname(Base.active_project())
    try
        environment_dir = joinpath(working_dir, "environment")
        repos_dir = joinpath(working_dir, "repos")
        mkpath(repos_dir)

        required_packages = get_deps_from_toml(joinpath(package_dir, "Project.toml"))
        local_package_names = get_package_name.(local_package_repos)
        nonlocal_package_names = [a for a in required_packages if !(a in local_package_names)]

        Pkg.activate(environment_dir)
        Pkg.instantiate()

        # You cannot just look at Project.toml to find the local
        # dependencies, because if one of those local packages has another local dependency,
        # that will not be found.
        for r in local_package_repos
            add_package(r)
        end
        Pkg.add(nonlocal_package_names)

        # now install the main package and test it
        dev_package(package_dir) # this is the main package
        package_name = extract_package_name_from_path(package_dir)
        Pkg.test(package_name)

    finally
        # revert back to usual dir
        Pkg.activate(save_project)
    end
    return
end

function main_test()
    home = rstrip(ENV["HOME"], ['/', ' '])
    working_dir = joinpath(home, "transfer")
    package_dir = joinpath(home, "opensource/github/bittide/Callisto.jl")
    local_package_repos = deps_for_callisto()
    build_and_test_package(working_dir, package_dir, local_package_repos)
end

function deps_for_callisto()
    b(x,y) = github_repo("bittide/$(x).jl", y)
    b(x) = github_repo("bittide/$(x).jl")
    local_package_repos = [
        b("JuliaTools"),
        b("PlotKitCairo"),
        b("PlotKitAxes", "dev"),
        b("PlotKitDiagrams", "dev"),
        b("PlotKitCharts", "dev"),
        b("PlotKit", "dev"),
        b("LogData"),
        b("Piecewise"),
        b("Topology"),
        b("StateSystems"),
        b("BittideTools"),
    ]
    return local_package_repos
end

function main()
    println("current working directory = ", pwd())
    println("current active project = ", Base.active_project())
    working_dir = joinpath(pwd(), "working")
    package_dir = pwd()
    local_package_repos = deps_for_callisto()
    build_and_test_package(working_dir, package_dir, local_package_repos)
end


end


using .Build
Build.main()
