
module Build

using Pkg

local_packages() = ["opensource/github/bittide/JuliaTools.jl",
                    "opensource/github/bittide/PlotKitCairo.jl",
                    "opensource/github/bittide/PlotKitGL.jl",
                    "opensource/github/bittide/PlotKitAxes.jl",
                    "opensource/github/bittide/PlotKitDiagrams.jl",
                    "opensource/github/bittide/PlotKitCharts.jl",
                    "opensource/github/bittide/PlotKit.jl",
                    "opensource/github/bittide/PlotKitPresentations.jl",
                    "opensource/github/bittide/GoogleColors.jl",
                    "opensource/github/bittide/LogData.jl",
                    "opensource/github/bittide/Piecewise.jl",
                    "opensource/github/bittide/CallistoLinear.jl",
                    "opensource/github/bittide/Topology.jl",
                    "opensource/github/bittide/StateSystems.jl",
                    "opensource/github/bittide/BittideTools.jl",
                    "opensource/github/bittide/Callisto.jl",
                    "opensource/github/bittide/Callisto2.jl",
                    "opensource/github/bittide/CallistoVisualization.jl",
                    "opensource/github/bittide/StandardPlots.jl",
                    "opensource/github/bittide/ReadElf.jl",
                    "opensource/github/bittide/RiscVSim.jl",
                    "opensource/github/bittide/Elara.jl",
                    ]

is_package_installed(pname) = pname in keys(Pkg.project().dependencies)

function install_package(package_name::String)
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

function main()
    println("current working directory = ", pwd())
    println("current active project = ", Base.active_project())
    for a in readdir()
        println(a)
    end
    Pkg.Registry.add()
    #
    # todo: remove Project.toml before this starts.
    #
    dev_package("others/StateSystems.jl")
  
end



end

using .Build
Build.main()
