### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 6d6b6472-4d5e-11eb-0e08-b3c26ad37f08
using Perfusion

# ╔═╡ c5497f8e-4d4e-11eb-3087-29485fb0bee1
begin
	using MAT, Plots, PlutoUI, Statistics
	gr(size=(600,400), html_output_format=:png)
end

# ╔═╡ 55f56d58-4d64-11eb-3c16-f56d0e192f6c
md"""
### Residual Sum of Squares
"""

# ╔═╡ 5fb9eb5c-4d64-11eb-2d4d-99d72a2140c5
md"""
# END

Everything below are the cogs behind the machine.
"""

# ╔═╡ 97f67b12-4d5e-11eb-332d-d9521cadce9d
figopts = (framestyle = :grid, gridalpha=0.5, gridstyle=:dot, linewidth = 3.5, 
        tickfontsize = 10, fg_legend = :transparent, legendfontsize = 9, 
		legend = :bottomright)

# ╔═╡ 583e98ec-4d5e-11eb-373f-adcf5b0f16be
modelfunc = Dict(
    "extendedtofts" => model_tofts,
    "tofts" => model_tofts,
    "uptake" => model_uptake,
    "exchange" => model_exchange
)

# ╔═╡ ba0a75ec-4d4e-11eb-383f-0bd33c2f5c5e
mat = matread("dro.mat");

# ╔═╡ 98c1115e-4d51-11eb-2880-65657a732856
allmodels = collect(keys(mat["fits"]))

# ╔═╡ 38393f20-4d54-11eb-0248-e3937edde1f7
md"""
# Fits on Digital Reference Object (DRO)

This document contains three figures.

1. Parameter maps for a model chosen from the list below.
  + If multiple models are chosen, then the first selected model will be shown.
  + The parameter map can be changed using the dropdown list below the figure

2. The fitted curve for each chosen model
  + The chosen voxel is shown by the crosshairs in the parameter map
  + The sliders will move the crosshairs

3. The residual sum of squares (static figure)

### Choose which models to show

The models represent:

- `exchange`: Two compartment exchange model (2CXM)
  + Parameters: Fp, PS, ve, vp, T, Te, Tp
- `extendedtofts`: Extended Tofts model
  + Parameters: Kt, ve, vp, kep
- `uptake`: Compartmental tissue uptake model (CTUM)
  + Parameters: Fp, PS, vp
- `tofts`: Standed Tofts Model
  + Parameters: Kt, ve, vp

Multiple models can be selected by holding ctrl while clicking.

$(@bind models MultiSelect(allmodels, default = allmodels))
"""

# ╔═╡ e395e576-4d63-11eb-1ca4-556bafe55491
let
	p = []
	cmax = 0.0005
	for model in allmodels
		map = mat["rss"][model]
		push!(p, heatmap(map; clim=(0, cmax), title = model))
	end
	plot(p...)
end

# ╔═╡ ae8f734c-4d54-11eb-2a2b-f33ec0accbe0
t = mat["t"];

# ╔═╡ cbc0af94-4d5e-11eb-3e9d-61955bc84d32
ca = mat["ca"];

# ╔═╡ c08e627a-4d53-11eb-0dce-639de5501eb4
ct = mat["ct"];

# ╔═╡ a894973e-4d53-11eb-121d-194fe32116fe
(nx, ny, nt) = size(ct)

# ╔═╡ caf170f4-4d53-11eb-1a6a-bb3fd836acdc
(midx, midy, midt) = round.(Int, [nx, ny, nt] ./ 2)

# ╔═╡ 1572b524-4d52-11eb-22a7-09978ae54280
fits = mat["fits"];

# ╔═╡ bb48c55a-4d51-11eb-2e44-db95a10f4e7a
allparams = collect(keys(fits[first(models)]))

# ╔═╡ 8a07502c-4d53-11eb-0acc-7bba93ebb256
md"""
$(@bind param Select(allparams, default = last(allparams)))
$(@bind y Slider(3:ny-3; default = midy))
$(@bind x Slider(3:nx-3; default = midx))
"""

# ╔═╡ 514a525c-4d50-11eb-07c6-97593cc5aec5
let
	model = first(models)
	map = copy(fits[model][param])
	cmax = quantile(filter(!isnan, map), 0.9)
	
	map[x-3:x-2, y] .= 0
    map[x+2:x+3, y] .= 0
    map[x, y-3:y-2] .= 0
    map[x, y+2:y+3] .= 0
	
	heatmap(map; yflip = true, clims = (0, cmax), title = "$model: $param")
end

# ╔═╡ b180e12a-4d5f-11eb-1636-5939375ec6dd
md"""
Fitted values at voxel (x = $x, y = $y)
"""

# ╔═╡ 9c775c6a-4d54-11eb-2f23-d3b183c95a4f
let
	conc = ct[x,y,:]
	p = scatter(t, conc; ylabel = "Time [min]", xlabel = "Concentration [mM]", label = nothing, palette = :seaborn_colorblind, figopts...)
	
	ylim = 1.1 .* extrema(conc)
	
	for model in models
        names = collect(keys(fits[model]))
        values = [fits[model][name][x,y] for name in names]
        params = NamedTuple{Tuple(Symbol.(names))}(values)    
        if !isnan(values[1])
            fittedcurve = modelfunc[model](; t, ca, params)
            plot!(t, fittedcurve; ylim, label = model, figopts...)
        end
    end
	p
end

# ╔═╡ Cell order:
# ╟─38393f20-4d54-11eb-0248-e3937edde1f7
# ╟─514a525c-4d50-11eb-07c6-97593cc5aec5
# ╟─8a07502c-4d53-11eb-0acc-7bba93ebb256
# ╟─b180e12a-4d5f-11eb-1636-5939375ec6dd
# ╟─9c775c6a-4d54-11eb-2f23-d3b183c95a4f
# ╟─55f56d58-4d64-11eb-3c16-f56d0e192f6c
# ╟─e395e576-4d63-11eb-1ca4-556bafe55491
# ╟─5fb9eb5c-4d64-11eb-2d4d-99d72a2140c5
# ╟─97f67b12-4d5e-11eb-332d-d9521cadce9d
# ╟─bb48c55a-4d51-11eb-2e44-db95a10f4e7a
# ╟─98c1115e-4d51-11eb-2880-65657a732856
# ╟─583e98ec-4d5e-11eb-373f-adcf5b0f16be
# ╠═a894973e-4d53-11eb-121d-194fe32116fe
# ╠═caf170f4-4d53-11eb-1a6a-bb3fd836acdc
# ╠═ae8f734c-4d54-11eb-2a2b-f33ec0accbe0
# ╠═cbc0af94-4d5e-11eb-3e9d-61955bc84d32
# ╠═c08e627a-4d53-11eb-0dce-639de5501eb4
# ╠═1572b524-4d52-11eb-22a7-09978ae54280
# ╠═ba0a75ec-4d4e-11eb-383f-0bd33c2f5c5e
# ╠═6d6b6472-4d5e-11eb-0e08-b3c26ad37f08
# ╠═c5497f8e-4d4e-11eb-3087-29485fb0bee1
