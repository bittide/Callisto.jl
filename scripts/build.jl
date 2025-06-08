
module Build

using Pkg
using Pkg.TOML

mutable struct Repo
    url::String
    branch::String
end


github_repo(reponame, branch) = Repo("https://github.com/$(reponame).git", branch)
github_repo(reponame) = github_repo(reponame, "main")

#

# PUBLIC repos. Callisto or PlotKit have these as dependencies
    # BittideTools
    # JuliaTools
    # LogData
    # Piecewise
    # PlotKitAxes
    # PlotKitCairo
    # PlotKitCharts
    # PlotKitDiagrams
    # StateSystems
    # Topology
    # PlotKit
    # Callisto

# PRIVATE repos
  # These will be needed for testing because
  # they are deps of Elara. So either the CI for Elara
  # uses a token or we make them public.
    # ReadElf
    # RiscVSim

  # To be released to users, but private for now:
  # These have no dependents
    # Callisto2
    # CallistoLinear
    # CallistoVisualization
    # Elara

  # Internal tools which do not need to be public
  # These have no dependents
    # StandardPlots
    # PlotKitGL
    # PlotKitPresentations
    # GoogleColors

# use: PkgDependency.tree("JuliaTools"; reverse = true)
# to see the dependents.





all_available_local_packages() = [
    github_repo("bittide/JuliaTools.jl"),
    github_repo("bittide/PlotKitCairo.jl"),
#    github_repo("bittide/PlotKitGL.jl"),
    github_repo("bittide/PlotKitAxes.jl", "dev"),
    github_repo("bittide/PlotKitDiagrams.jl", "dev"),
    github_repo("bittide/PlotKitCharts.jl", "dev"),
    github_repo("bittide/PlotKit.jl", "dev"),
#    github_repo("bittide/PlotKitPresentations.jl"),
#    github_repo("bittide/GoogleColors.jl"),
    github_repo("bittide/LogData.jl"),
    github_repo("bittide/Piecewise.jl"),
#    github_repo("bittide/CallistoLinear.jl"),
    github_repo("bittide/Topology.jl"),
    github_repo("bittide/StateSystems.jl"),
    github_repo("bittide/BittideTools.jl"),
    github_repo("bittide/Callisto.jl", "dev"),
#    github_repo("bittide/Callisto2.jl"),
#    github_repo("bittide/CallistoVisualization.jl", "timingplot"),
#    github_repo("bittide/StandardPlots.jl"),
#    github_repo("bittide/ReadElf.jl"),
#    github_repo("bittide/RiscVSim.jl"),
#    github_repo("bittide/Elara.jl"),
]



function git_clone(repo::Repo, dirname = "")
    git_exe = Sys.which("git")
    if isnothing(git_exe)
        error("git command not found in PATH.")
    end
    cmd = `$(git_exe) clone -b $(repo.branch) $(repo.url) $(dirname)`
    run(cmd)
end

match_repo_name(s) = match(r"(/)([a-zA-Z0-9_-]+)(\.jl)", s)
match_repo(r::Repo) = match_repo_name(r.url)
get_package_name(r::Repo) = extract_package_name_from_path(r.url)

function get_local_directory(r::Repo, repodir)
    m = match_repo(r)
    if !isnothing(m)
        return joinpath(repodir, m[2] * ".jl")
    end
    return nothing
end

"""
    get_repo(package_repo::Repo, repodir)

Changes the current directory to `dirname` and checks if a directory corresponding
to the package name derived from `local_package_repo` exists. If the directory
does not exist, it attempts to clone the repository specified by `local_package_repo`.

"""
function get_repo(package_repo::Repo, repodir)
    local_directory = get_local_directory(package_repo, repodir)
    if !isnothing(local_directory) && !isdir(local_directory)
        git_clone(package_repo, local_directory)
    end
    return local_directory
end


"""
    make_package_list(project_toml_file::String) -> Vector{String}

Parses a `Project.toml` file and returns a list of package names
listed under its `[deps]` section.
"""
function make_package_list(project_toml_file)
    parsed_toml = TOML.parsefile(project_toml_file)
    deps = get(parsed_toml, "deps", Dict{String, Any}())
    return collect(keys(deps))
end

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
    separate_local_and_nonlocal_packages(required_packages)

Given a list of `required_packages` names, this function categorizes them into
local and non-local packages based on `all_available_local_packages()`.

Returns:
- A list of repository paths for local packages (e.g., "bittide/PackageName.jl").
- A list of package names for non-local packages (assumed to be in the general registry).

The returned local package repos are in the same order as they are in
all_available_local_packages(). The idea is they are ordered by dependency,
so can be installed in that order.
"""
function separate_local_and_nonlocal_packages(required_packages)
    local_package_repos = []
    local_package_names = []
    nonlocal_package_names = []

    for package_repo in all_available_local_packages()
        package_name = get_package_name(package_repo)
        if package_name in required_packages
            push!(local_package_repos, package_repo)
            push!(local_package_names, package_name)
        end
    end
    # now find remaining packages and they must be
    # in the general registry
    for package in required_packages
        if !(package in local_package_names)
            push!(nonlocal_package_names, package)
        end
    end
    return local_package_repos, nonlocal_package_names
end



"""
    install_local_packages(repos_dir, local_package_repos)

For each package repository in `local_package_repos`:
1. Clone the repo into `repos_dir`
2. Adds the package to the current environment in development mode (`Pkg.develop`).
"""
function install_local_packages(repos_dir, local_package_repos)
    for package_repo in local_package_repos
        package_path = get_repo(package_repo, repos_dir)
        dev_package(package_path)
    end
end



function install_all_local_packages(repos_dir)
    local_package_repos = all_available_local_packages()
    install_local_packages(repos_dir, local_package_repos)
end


function install_nonlocal_packages(nonlocal_packages)
    for package in nonlocal_packages
        install_package(package)
    end
end


function dev_package(package_path)
    try
        Pkg.develop(PackageSpec(path = package_path))
    catch e
        println("Error installing (dev) $package_path: $e")
    end
end

function install_package(package_name)
    try
        Pkg.add(package_name)
    catch e
        println("Error installing $package_name: $e")
    end
end

function mkpath_if_not_exists(path)
    if !isdir(path)
        mkpath(path)
    end
end


"""
    build_and_test_package(working_dir, package_dir)

Sets up an environment, installs dependencies, and tests a given Julia package.

Example:
  package_dir = "/home/somewhere/Callisto.jl")
  working_dir = "/home/emptydir"

The package_dir is not changed, but is used to find the Project.toml
and install the required dependencies.  Then the package is dev'ed at
that directory and tested. Dependencies and a Julia environment is stored
in the working_dir, which should be empty. The working_dir will be created
if it does not exist.

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
function build_and_test_package(working_dir, package_dir)
    save_project = dirname(Base.active_project())
    try
        environment_dir = joinpath(working_dir, "environment")
        repos_dir = joinpath(working_dir, "repos")
        mkpath_if_not_exists(repos_dir)

        required_packages = make_package_list(joinpath(package_dir, "Project.toml"))
        local_package_repos, nonlocal_package_names = separate_local_and_nonlocal_packages(required_packages)

        Pkg.activate(environment_dir)
        Pkg.instantiate()

        # You cannot just do install_local_packages(repos_dir, local_package_repos)
        # because if one of those local packages has another local dependency,
        # that will not be found
        install_all_local_packages(repos_dir)
        install_nonlocal_packages(nonlocal_package_names)

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
    build_and_test_package(working_dir, package_dir)
end


function main()
    println("current working directory = ", pwd())
    println("current active project = ", Base.active_project())
    working_dir = joinpath(pwd(), "working")
    package_dir = pwd()
    build_and_test_package(working_dir, package_dir)
end


end


using .Build
Build.main()
