using CSV
using DataFrames
using GLMakie

const ROOT = @__DIR__
const EXTRACTED = joinpath(ROOT, "extracted")
const XLSX_CACHE = joinpath(ROOT, "xlsx_cache")
const OUT_PNG = joinpath(ROOT, "akira_all_data_overview.png")

const TXT_FILES = [
    "130618-03 3 traces.txt",
    "130618-03 Si1.txt",
    "130618-03 Si2.txt",
    "130618-03 Si3.txt",
    "Fig5A.txt",
    "JN2016_Fig5A_150413_01_4494.txt",
    "JN2016_Fig5B_130917_01_3351.txt",
    "JN2016_Fig5C_130917_01_3310.txt",
    "JN2016_Fig5E_150413_01_3856.txt",
    "JN2016_Fig5F_150413_01_3714.txt",
    "JN2016_Fig5G_150413_01_3460.txt",
]

const TSV_FILES = [
    "Fig_2Aii.tsv",
    "Fig_2Aiii.tsv",
    "Fig_2_Aiv.tsv",
    "Fig_2B.tsv",
    "Fig_2C.tsv",
    "Fig_2D.tsv",
]

const TRACE_COLORS = [
    RGBf(0.12, 0.47, 0.71),
    RGBf(0.17, 0.63, 0.17),
    RGBf(0.84, 0.15, 0.16),
    RGBf(0.58, 0.40, 0.74),
    RGBf(0.55, 0.34, 0.29),
    RGBf(0.89, 0.47, 0.76),
]

clean_label(x) = replace(strip(String(x), ['"', '\'']), "_" => " ")

function crop_png_whitespace!(path::AbstractString)
    ps = """
    Add-Type -AssemblyName System.Drawing
    \$path = '$path'
    \$img = [System.Drawing.Bitmap]::FromFile(\$path)
    [int]\$minX = \$img.Width
    [int]\$minY = \$img.Height
    [int]\$maxX = -1
    [int]\$maxY = -1
    for (\$y = 0; \$y -lt \$img.Height; \$y++) {
        for (\$x = 0; \$x -lt \$img.Width; \$x++) {
            \$c = \$img.GetPixel(\$x, \$y)
            if ((\$c.R -ne 255) -or (\$c.G -ne 255) -or (\$c.B -ne 255)) {
                if (\$x -lt \$minX) { \$minX = \$x }
                if (\$y -lt \$minY) { \$minY = \$y }
                if (\$x -gt \$maxX) { \$maxX = \$x }
                if (\$y -gt \$maxY) { \$maxY = \$y }
            }
        }
    }
    if (\$maxX -ge \$minX -and \$maxY -ge \$minY) {
        [int]\$width = \$maxX - \$minX + 1
        [int]\$height = \$maxY - \$minY + 1
        \$rect = New-Object System.Drawing.Rectangle([int]\$minX, [int]\$minY, [int]\$width, [int]\$height)
        \$cropped = \$img.Clone(\$rect, \$img.PixelFormat)
        \$img.Dispose()
        \$cropped.Save(\$path, [System.Drawing.Imaging.ImageFormat]::Png)
        \$cropped.Dispose()
    } else {
        \$img.Dispose()
    }
    """
    run(`powershell -NoProfile -Command $ps`)
end

function first_line_has_header(path::AbstractString)
    open(path, "r") do io
        line = readline(io)
        any(isletter, line)
    end
end

function load_txt_table(path::AbstractString)
    if first_line_has_header(path)
        df = DataFrame(CSV.File(path; delim='\t', normalizenames=false))
        names_clean = [clean_label(name) for name in names(df)]
        rename!(df, Symbol.(names_clean))
        return names_clean, Matrix{Float64}(df)
    end

    df = DataFrame(CSV.File(path; delim='\t', header=false, normalizenames=false))
    names_clean = ["col$(i)" for i in 1:ncol(df)]
    rename!(df, Symbol.(names_clean))
    return names_clean, Matrix{Float64}(df)
end

function downsample(x::AbstractVector, ys::Vector{<:AbstractVector}; max_points::Int=8000)
    n = length(x)
    stride = max(1, ceil(Int, n / max_points))
    idx = 1:stride:n
    return x[idx], [y[idx] for y in ys]
end

function finite_pair(x::AbstractVector, y::AbstractVector)
    mask = isfinite.(x) .& isfinite.(y)
    return x[mask], y[mask]
end

function plot_txt_panel!(ax::Axis, filename::String)
    path = joinpath(EXTRACTED, filename)
    headers, data = load_txt_table(path)
    ncols = size(data, 2)

    if iseven(ncols) && all(occursin("Time", headers[i]) for i in 1:2:ncols)
        for pair_idx in 1:(ncols ÷ 2)
            x = data[:, 2pair_idx - 1]
            y = data[:, 2pair_idx]
            x, y = finite_pair(x, y)
            isempty(x) && continue
            x = x .- first(x)
            x_ds, y_ds_vec = downsample(x, [y])
            lines!(ax, x_ds, y_ds_vec[1], color=TRACE_COLORS[mod1(pair_idx, length(TRACE_COLORS))], linewidth=1.0)
        end
    else
        for y_idx in 2:ncols
            x = data[:, 1]
            y = data[:, y_idx]
            x, y = finite_pair(x, y)
            isempty(x) && continue
            x = x .- first(x)
            x_ds, y_ds_vec = downsample(x, [y])
            lines!(ax, x_ds, y_ds_vec[1], color=TRACE_COLORS[mod1(y_idx - 1, length(TRACE_COLORS))], linewidth=1.0)
        end
    end

    ax.title = replace(splitext(filename)[1], "_" => " ")
    ax.xlabel = "Time (s)"
    ax.ylabel = "Value"
    ax.titlesize = 16
    ax.xlabelsize = 12
    ax.ylabelsize = 12
    ax.xticklabelsize = 10
    ax.yticklabelsize = 10
end

function read_tsv(path::AbstractString)
    df = DataFrame(CSV.File(path; delim='\t', normalizenames=false))
    return [String(n) for n in names(df)], df
end

function numeric_series(df::DataFrame, xcol::Int, ycol::Int)
    xs = Float64[]
    ys = Float64[]
    for row in eachrow(df)
        x = row[xcol]
        y = row[ycol]
        if ismissing(x) || ismissing(y) || x == "" || y == ""
            continue
        end
        try
            push!(xs, x isa AbstractString ? parse(Float64, x) : Float64(x))
            push!(ys, y isa AbstractString ? parse(Float64, y) : Float64(y))
        catch
        end
    end
    return xs, ys
end

function plot_tsv_panel!(ax::Axis, filename::String)
    path = joinpath(XLSX_CACHE, filename)
    headers, df = read_tsv(path)
    pairs = if filename == "Fig_2Aiii.tsv"
        [(1, 2), (4, 5)]
    elseif filename in ("Fig_2_Aiv.tsv", "Fig_2D.tsv")
        [(1, 2), (3, 4)]
    else
        [(1, 2)]
    end

    for (i, (xidx, yidx)) in enumerate(pairs)
        x, y = numeric_series(df, xidx, yidx)
        isempty(x) && continue
        scatter!(ax, x, y; color=(TRACE_COLORS[mod1(i, length(TRACE_COLORS))], 0.75), markersize=7)
    end

    ax.title = replace(splitext(filename)[1], "_" => " ")
    ax.xlabel = clean_label(headers[first(first(pairs))])
    ax.ylabel = "Value"
    ax.titlesize = 16
    ax.xlabelsize = 12
    ax.ylabelsize = 12
    ax.xticklabelsize = 10
    ax.yticklabelsize = 10
end

function build_figure()
    n_panels = length(TXT_FILES) + length(TSV_FILES)
    ncols = 3
    nrows = ceil(Int, n_panels / ncols)
    fig = Figure(size=(2200, 3300), backgroundcolor=:white)
    Label(fig[1, :], "Akira bundle overview: extracted numeric datasets", fontsize=26, font=:bold, padding=(0, 0, 10, 10))

    panel_idx = 1
    for filename in TXT_FILES
        row = ceil(Int, panel_idx / ncols) + 1
        col = mod1(panel_idx, ncols)
        ax = Axis(fig[row, col], xgridvisible=true, ygridvisible=true, backgroundcolor=RGBf(0.99, 0.99, 0.99))
        plot_txt_panel!(ax, filename)
        panel_idx += 1
    end

    for filename in TSV_FILES
        row = ceil(Int, panel_idx / ncols) + 1
        col = mod1(panel_idx, ncols)
        ax = Axis(fig[row, col], xgridvisible=true, ygridvisible=true, backgroundcolor=RGBf(0.99, 0.99, 0.99))
        plot_tsv_panel!(ax, filename)
        panel_idx += 1
    end

    while panel_idx <= nrows * ncols
        row = ceil(Int, panel_idx / ncols) + 1
        col = mod1(panel_idx, ncols)
        Box(fig[row, col], color=:white, strokewidth=0)
        panel_idx += 1
    end

    resize_to_layout!(fig)
    save(OUT_PNG, fig)
    crop_png_whitespace!(OUT_PNG)
    return fig
end

const AKIRA_OVERVIEW_NOAUTORUN = get(ENV, "AKIRA_OVERVIEW_NOAUTORUN", "0") == "1"

if !AKIRA_OVERVIEW_NOAUTORUN
    fig = build_figure()
    display(fig)
end
