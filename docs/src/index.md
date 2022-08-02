# Callisto.jl

## Introduction

`Callisto.jl` is a simulator for [bittide](https://www.bittide.io/), a system
architecture for synchronous distributed computing.  A detailed description of
the bittide mechanism is available in the papers:

  * S. Lall, C. Cascaval, M. Izzard and T. Spalink. Modeling and Control of bittide Synchronization.  Proceedings of the American Control Conference, 2022. [arXiv:2109.14111](https://arxiv.org/abs/2109.14111)
    
  * S. Lall, C. Cascaval, M. Izzard and T. Spalink. Resistance Distance and	Control Performance for bittide Synchronization. Proceedings of the European Control Conference, 2022. [arXiv:2111.05296](https://arxiv.org/abs/2111.05296)

The model and algorithm used by Callisto follows closely that described in the
first of these papers. 

## Resources

* The [source code for Callisto](https://github.com/bittide/Callisto.jl)
* The [documentation](https://www.bittide.io/Callisto.jl/)
* The [CallistoVisualization package](https://github.com/bittide/CallistoVisualization.jl)
* The [PlotKit package](https://github.com/bittide/PlotKit.jl)



## Installation 

Start Julia, and at the REPL prompt install Callisto as follows.

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/bittide/Callisto.jl")
```

## Quickstart

Now you can run a simple simulation:

```julia
julia> using Callisto
julia> c = CalOpts()
julia> x = callisto(c)
```

You can view the output:

```julia
julia> xc = parse_callisto_logx(c, x...)
julia> using PlotKit
julia> pzip(a) = Point.(zip(a.x, a.y))
julia> save(drawplot(pzip.(xc.freq)), "frequency.png")
```

This will save a plot of the frequency against time in the file `frequency.png` in 
the Julia working directory.

## Common options


```julia
julia> c = CalOpts(; kwargs...)
```

`CalOpts` takes a number of optional keyword arguments. These include

- `topology`: a tuple, whose first entry is the string name, and remaining
   entries are numerical parameters. For example,
     + `("triangle")`
     + `("mesh", 3, 4)`
     + `("line", 6)`
- `kp`: proportional gain
- `ki`: integral gain
- `latency`: either a number, in which case it is the latency of every link,
   or a list, whose `i`th entry is the latency of the `i`th link
- `tmax`: simulation end time
- `poll_period`: the delay, in localticks, between successive control updates
- `control_delay`: the delay, in localticks, at a node between making a measurement
   and the frequency update
   
   
   






