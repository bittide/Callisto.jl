
module Build

using Pkg
using Pkg.TOML
local_packages() = ["bittide/JuliaTools.jl",
                    "bittide/PlotKitCairo.jl",
                    "bittide/PlotKitGL.jl",
                    "bittide/PlotKitAxes.jl",
                    "bittide/PlotKitDiagrams.jl",
                    "bittide/PlotKitCharts.jl",
                    "bittide/PlotKit.jl",
                    "bittide/PlotKitPresentations.jl",
                    "bittide/GoogleColors.jl",
                    "bittide/LogData.jl",
                    "bittide/Piecewise.jl",
                    "bittide/CallistoLinear.jl",
                    "bittide/Topology.jl",
                    "bittide/StateSystems.jl",
                    "bittide/BittideTools.jl",
                    "bittide/Callisto.jl",
                    "bittide/Callisto2.jl",
                    "bittide/CallistoVisualization.jl",
                    "bittide/StandardPlots.jl",
                    "bittide/ReadElf.jl",
                    "bittide/RiscVSim.jl",
                    "bittide/Elara.jl",
                    ]

is_package_installed(pname) = pname in keys(Pkg.project().dependencies)

function dev_package(package_name::String)
    if !is_package_installed(package_name)
        println("Installing: $package_name")
        Pkg.add(package_name)
    end
end

function dev_package(dir_name::String)
    #home = rstrip(ENV["HOME"], ['/', ' '])
    home = rstrip(pwd(), ['/', ' '])
    m = match(r"(/)([a-zA-Z0-9]*)(\.jl)", dir_name)
    if isnothing(m)
        println("Not installing", dir_name)
        return
    end
    package_name = m[2]
    #    if !is_package_installed(package_name)
    println("Installing (dev): $package_name")
    Pkg.develop(PackageSpec(path = joinpath(home, dir_name)))
    #end
end

function git_clone(reponame)
    cmd = `/usr/bin/git clone https://github.com/$(reponame).git`
    run(cmd)
end

function get_repos(dirname)
    cd(dirname)
    for package in local_packages()
        m = match(r"(/)([a-zA-Z0-9]*)(\.jl)", package)
        if !isnothing(m) && !isdir(m[2] * ".jl")
            git_clone(package)
        end
    end
end

function mainx()
    println("current working directory = ", pwd())
    println("current active project = ", Base.active_project())
    for a in readdir()
        println(a)
    end
    Pkg.Registry.add()
    mkdir("others")
    cd("others")
    git_clone("bittide/JuliaTools")
    git_clone("bittide/StateSystems")
    dev_package("others/StateSystems.jl")
    dev_package("others/StateSystems.jl")
end

function make_package_list(project_toml_file)
    parsed_toml = TOML.parsefile(project_toml_file)
    package_list = keys(parsed_toml["deps"])
    return collect(package_list)
end

function extract_package_name_from_path(package_path)
    m = match(r"(/)([a-zA-Z0-9]*)(\.jl)", package_path)
    if isnothing(m)
        return ""
    end
    return m[2]
end

function install_all_packages(environment_dir, repos_dir, package_list)
    Pkg.activate(environment_dir)
    Pkg.instantiate()
    # The assumption here is that local_packages are in
    # dependency order, so that we can install them in that order.
    installed_packages = Dict()
    for package_repo in local_packages()
        package_name = extract_package_name_from_path(package_repo)
        if package_name in package_list
            package_path = joinpath(repos_dir, package_name * ".jl")
            dev_package(package_path)
            installed_packages[package_name] = true
        end
    end
    # now find remaining packages and they must be
    # in the general registry
    for package in package_list
        if !(package in keys(installed_packages))
            install_package(package)
        end
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

function main_test()
    home = rstrip(ENV["HOME"], ['/', ' '])
    repos_dir = joinpath(home, "transfer")
    package_dir = joinpath(home, "opensource/github/bittide/Callisto.jl")
    environment_dir = joinpath(home, "transfer/environment")

    get_repos(repos_dir)
    package_list = make_package_list(joinpath(package_dir, "Project.toml"))
    install_all_packages(environment_dir, repos_dir, package_list)
    dev_package(package_dir) # this is the main package
    package_name = extract_package_name_from_path(package_dir)
    Pkg.test(package_name)

    # revert back to usual dir
    Pkg.activate(joinpath(home, "opensource/github/bittide/experimental_scripts/src/sandbox/Sandbox"))
end


end

using .Build
Build.main_test()
