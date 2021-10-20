# NuclearReactionNetwork.lj Documentation

```@contents
```

## LinearInterpolations

```@meta
CurrentModule = LinearInterpolations
```

```@docs
LinearInterpolations
NcapLerp
TrajectoryLerp
get_rate(ncap_lerp::NcapLerp, temp::Float64)
get_pfunc(ncap_lerp::NcapLerp, temp::Float64)
get_temp_density(traj_lerp::TrajectoryLerp, time::Float64)
```

## LinearSolvers

```@meta
CurrentModule = LinearSolvers
```

```@docs
LinearSolvers
LinearSolver
LS_UMFPACK
LS_MKLPardisoSolver
solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LinearSolver)
solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LS_UMFPACK)
solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LS_MKLPardisoSolver)
```

## Astro

```@meta
CurrentModule = Astro
```

```@docs
Astro
CurrentTrajectory
get_current_trajectory(trajectory::TrajectoryLerp, time::Float64)
```

## Network

```@meta
CurrentModule = Network
```

```@docs
Network
NetworkBoundary
NetworkIndex
NetworkIndex(networkboundary::NetworkBoundary)
step_time!(time::Time)
get_networksize(net_idx::NetworkIndex)
get_networksize(networkboundary::NetworkBoundary)
zn_to_index(z::Int, n::Int, net_idx::NetworkIndex)
zn_in_network(z::Int, n::Int, net_idx::NetworkIndex)
Time
step_time!
```

## ReactionTypes

```@meta
CurrentModule = ReactionTypes
```

```@docs
ReactionTypes
ProbDecay
NeutronCapture
AlphaDecay
Photodissociation
ReactionData
initialize_reactions()
check_eq_reaction(lhs::AbstractReaction, rhs::AbstractReaction)
```

## NetworkDatas

```@meta
CurrentModule = NetworkDatas
```

```@docs
NetworkDatas
NetworkData
OutputInfo
IncludedReactions
fill_jacobian!(nd::NetworkData; use_yproposed::Bool=false)
update_ydot!(nd::NetworkData; use_yproposed::Bool=false)
```

## InOut

```@meta
CurrentModule = InOut
```

```@docs
InOut
Result
dump_y(nd::NetworkData)
dump_ya(nd::NetworkData)
dump_result(nd::NetworkData)
dump_iteration(nd::NetworkData, itation::Int)
initialize_network_data(path::String)
```

## NetworkSolve

```@meta
CurrentModule = NetworkSolve
```

```@docs
newton_raphson_iteration!(nd::NetworkData, Δy::Vector{Float64})
clip_abundance!(nd::NetworkData, min_val::Float64=0.0)
SolveNetwork!(nd::NetworkData)
```

## Index

```@index
```