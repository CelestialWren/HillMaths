using Plots

struct Hill
    x_coefficent::Float64
    y_coefficent::Float64
    max_height::Float64
    centre_offset::Tuple{Float64,Float64}
    x_range::StepRange
    y_range::StepRange

    Hill(
        x_coefficent::Real,
        y_coefficent::Real,
        max_height::Real,
        centre_offset::Tuple{<:Real,<:Real},
        x_domain::Tuple{<:Integer,<:Integer},
        y_domain::Tuple{<:Integer,<:Integer},
    ) =
        x_domain[1] > x_domain[2] ?
        error("The lower bound of the domain for x must be less than the upper bound") :
        y_domain[1] > y_domain[2] ?
        error("The lower bound of the domain for y must be less than the upper bound") :
        new(
            x_coefficent,
            y_coefficent,
            max_height,
            centre_offset,
            range(x_domain[1], x_domain[2], step = 1),
            range(y_domain[1], y_domain[2], step = 1),
        )
end

function hill_height(hill::Hill, x::Real, y::Real)::Float64
    hill.max_height - (
        hill.x_coefficent * (x - hill.centre_offset[1])^2 +
        hill.y_coefficent * (y - hill.centre_offset[2])^2
    )
end

function hill_height(hills::Tuple{Hill, Hill}, x::Real, y::Real)::Float64
    max(hills[1].max_height - (
        hills[1].x_coefficent * (x - hills[1].centre_offset[1])^2 +
        hills[1].y_coefficent * (y - hills[1].centre_offset[2])^2), (hills[2].max_height - (
            hills[2].x_coefficent * (x - hills[2].centre_offset[1])^2 +
            hills[2].y_coefficent * (y - hills[2].centre_offset[2])^2)))
end

function partial_x(hill::Hill, point::Tuple{<:Real,<:Real})::Float64
    if point[1] < hill.x_range[1] ||
       point[1] > last(hill.x_range) ||
       point[2] < hill.y_range[1] ||
       point[2] > last(hill.y_range)
        error("The point passed in must be in the range of X & Y values, it was $point")
    else
        2 * hill.centre_offset[1] * hill.x_coefficent - 2 * point[1] * hill.x_coefficent
    end
end

function partial_y(hill::Hill, point::Tuple{<:Real,<:Real})::Float64
    if point[1] < hill.x_range[1] ||
       point[1] > last(hill.x_range) ||
       point[2] < hill.y_range[1] ||
       point[2] > last(hill.y_range)
        error("The point passed in must be in the range of X & Y values, it was $point")
    else
        2 * hill.centre_offset[2] * hill.y_coefficent - 2 * point[2] * hill.y_coefficent
    end
end

function directional_derivative(
    hill::Hill,
    point::Tuple{<:Real,<:Real},
    angle_degrees::Real,
)::Float64
    partial_x(hill, point) * cosd(angle_degrees) +
    partial_y(hill, point) * sind(angle_degrees)
end

function gradient(hill::Hill, point::Tuple{<:Real,<:Real})::Tuple{Float64,Float64}
    (partial_x(hill, point), partial_y(hill, point))
end

function vec_to_polar(vector::Tuple{<:Real,<:Real})::Tuple{Float64,Float64}
    (sqrt(vector[1]^2 + vector[2]^2), acotd(vector[1] / vector[2]))
end

function estimate_closeness(
    hill::Hill,
    point::Tuple{<:Fractional,<:Fractional},
    grade_degrees::Float64,
)
    [
        (
            (cosd(angle), sind(angle)),
            abs(grade_degrees - atand(directional_derivative(hill, point, angle), sqrt(cosd(angle)^2 + sind(angle)^2))),
        ) for angle = 0:359]
end
function cell_to_move_to(
    hill::Hill,
    point::Tuple{<:Fractional,<:Fractional},
    grade_degrees::Float64,
)
    closeness_estimate::Vector{Tuple{Tuple{<:Fractional,<:Fractional},Float64}} = estimate_closeness(hill, point, grade_degrees)
    minima = 1
    for i in eachindex(closeness_estimate)
        if closeness_estimate[minima][2] > closeness_estimate[i][2]
            minima = i
        end
    end
    (
        closeness_estimate[minima][1][1] + point[1],
        closeness_estimate[minima][1][2] + point[2],
    )
end

function level_curve(hill::Hill, z::Real)
    local point = missing
    for val in hill.centre_offset[2]:1.0:convert(Float64,last(hill.y_range))
        if abs(z-hill_height(hill, hill.centre_offset[1], val)) < 1
            point = (hill.centre_offset[1], 2*hill.centre_offset[2]-val)
        end
    end
    if ismissing(point)
        error("could not find z-val")
    end
    is_loop = false
    curve::Vector{Tuple{<:Fractional,<:Fractional}} = [point]
    while (last(curve) != curve[1] || length(curve) == 1) &&
              let cell = cell_to_move_to(hill, last(curve), 0.0)
                cell[1] > hill.x_range[1] && cell[1] < last(hill.x_range) && cell[2] > hill.y_range[1] && cell[2] < last(hill.y_range)
              end &&
              length(curve) < 1000 && (abs(1+hill.centre_offset[1]-last(curve)[1]) > 1 || last(curve)[2] <= hill.centre_offset[2])
        push!(curve, cell_to_move_to(hill, last(curve), 0.0))
        is_loop = last(curve) == curve[1]
    end
    vcat(curve, reverse([(hill.centre_offset[1] + (hill.centre_offset[1]-points[1]), points[2]) for points = curve]))
end

function plot_level_curve(hill::Hill, z::Real)
    plot3d!(
        [
            (pair[1], pair[2], hill_height(hill, pair...)) for
            pair in level_curve(hill, z)
        ],
        linecolor = :black
    )
end

function path_of_ascent(hill::Hill, point::Tuple{<:Fractional, <:Fractional}) 
    is_loop = false
    curve::Vector{Tuple{<:Fractional,<:Fractional}} = [point]
    while (last(curve) != curve[1] || length(curve) == 1) &&
              let cell = cell_to_move_to(hill, last(curve), 90.0)
                cell[1] > hill.x_range[1] && cell[1] < last(hill.x_range) && cell[2] > hill.y_range[1] && cell[2] < last(hill.y_range)
              end &&
              length(curve) < 1000 && hill.max_height - hill_height(hill,last(curve)...) > 1
        push!(curve, cell_to_move_to(hill, last(curve), 90.0))
        is_loop = last(curve) == curve[1]
    end
  curve
end

function plot_path_of_ascent(hill::Hill, point::Tuple{<:Fractional, <:Fractional})
    plot3d!(
        [
            (pair[1], pair[2], hill_height(hill, pair...)) for
            pair in path_of_ascent(hill, point)
        ],
        linecolor = :black, legend = :none
    )
end

function walking_path(hill::Hill, point::Tuple{<:Fractional, <:Fractional}) 
    is_loop = false
    curve::Vector{Tuple{<:Fractional,<:Fractional}} = [point]
    while (last(curve) != curve[1] || length(curve) == 1) &&
              let cell = cell_to_move_to(hill, last(curve), 10.0)
                cell[1] > hill.x_range[1] && cell[1] < last(hill.x_range) && cell[2] > hill.y_range[1] && cell[2] < last(hill.y_range)
              end &&
              length(curve) < 100
        push!(curve, cell_to_move_to(hill, last(curve), 10.0))
        is_loop = last(curve) == curve[1]
    end
  curve
end

function plot_walking_path(hill::Hill, point::Tuple{<:Fractional, <:Fractional})
    plot3d!(
        [
            (pair[1], pair[2], hill_height(hill, pair...)) for
            pair in walking_path(hill, point)
        ],
        linecolor = :black, legend = :none
    )
end

function plot_hill(hill::Hill)
    surface(
        hill.x_range,
        hill.y_range,
        (x = hill.x_range, y = hill.y_range) -> hill_height(hill, x, y), zlims=(0,1.2*hill.max_height))
end
function plot_hill!(hill::Hill)
    surface!(
        hill.x_range,
        hill.y_range,
        (x = hill.x_range, y = hill.y_range) -> hill_height(hill, x, y), zlims=(0,1.2*hill.max_height))
end
plotlyjs()
hillA = Hill(.003, .005, 40, (0, 0), (-120, 120), (-200, 200))
hillB = Hill(.003, .005, 100, (350, 30), (100, 600), (-200, 200))
plot = plot_hill(hillA)
plot_hill!(hillB)
plot_level_curve(hillA, 30)
plot_path_of_ascent(hillA, (-43.0,-56.0))
plot_walking_path(hillA, (39.0,-42.0))
plot_level_curve(hillB, 60)
plot_path_of_ascent(hillB, (244.0,-50.0))
plot_walking_path(hillB, (355.0,-28.0))
display(plot)