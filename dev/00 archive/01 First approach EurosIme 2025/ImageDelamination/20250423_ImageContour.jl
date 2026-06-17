using Pkg
Pkg.add(["Images", "ImageIO", "ImageEdgeDetection", "ImageView"])
Pkg.add("ImageFiltering")
using Images, ImageIO, ImageFiltering, ImageView

# Sobel kernels
kx = [-1 0 1; -2 0 2; -1 0 1] ./ 4
ky = [-1 -2 -1; 0 0 0; 1 2 1] ./ 4

# Compute Sobel edge magnitude
function sobel_edges(img::AbstractMatrix{<:Gray}; threshold::Float64 = 0.5)
    gx = imfilter(img, kx)
    gy = imfilter(img, ky)
    edge_mag = sqrt.(gx.^2 .+ gy.^2)
    edge_mag ./= maximum(edge_mag)  # Normalize to [0, 1]
    contour = edge_mag .> threshold
    return contour
end

# Load image and apply
function extract_contour(image_path::String; threshold::Float64 = 0.5)
    img = Gray.(load(image_path))
    contour = sobel_edges(img; threshold=threshold)
    return contour
end

# Example usage
contour_image = extract_contour("Programs/20250304-B_x5_5000.png", threshold=0.2)

# View result
imshow(contour_image)

Pkg.add("ImageMorphology")

using Images, ImageIO, ImageMorphology, ImageView

"""
Extract the contour of a region based on intensity range.
- intensity_range: Tuple like (0.3, 0.7) to select pixels in that range
"""
function extract_region_contour(image_path::String; intensity_range=(0.3, 0.7))
    # Load and convert to grayscale
    img = Gray.(load(image_path))
    
    # Select region of interest by thresholding
    region_mask = (img .>= intensity_range[1]) .& (img .<= intensity_range[2])

    # Optional: fill holes and clean small components
    region_mask = imfill(region_mask, 4)
    region_mask = component_area_filter(region_mask, (area -> area > 50))

    # Extract the outer contour by subtracting eroded mask from original mask
    eroded = imerode(region_mask, ones(3,3))
    contour = region_mask .& .!eroded  # border = region - eroded interior

    return contour
end

# Example usage
contour_image = extract_region_contour("Programs/20250304-B_x5_5000.png", intensity_range=(0.3, 0.7))

# Show result
imshow(contour_image)

####
using Pkg
Pkg.add("StatsBase")
using Images, ImageIO, ImageMorphology, ImageView, StatsBase
using ImageMorphology: erode
"""
Filter out connected components smaller than min_area.
"""
function filter_large_components(mask::BitMatrix, min_area::Int)
    labeled = label_components(mask)
    counts = countmap(labeled)

    # Keep only component labels with area > min_area (excluding background 0)
    allowed_labels = Set(k for (k, v) in counts if k != 0 && v > min_area)
    return map(i -> i in allowed_labels, labeled)
end

"""
Extract the contour of a region based on intensity range.
"""
function extract_region_contour(image_path::String; intensity_range=(0.3, 7), min_area=50)
    img = Gray.(load(image_path))

    # Threshold by intensity
    region_mask = (img .>= intensity_range[1]) .& (img .<= intensity_range[2])

    # Fill small holes
    region_mask = imfill(region_mask, 4)

    # Filter out small regions
    region_mask = filter_large_components(region_mask, min_area)

    # Extract contour: region minus eroded version
   
   eroded = erode(region_mask, trues(3, 3))
    contour = region_mask .& .!eroded

    return contour
end

# Example usage:
contour = extract_region_contour("Programs/20250304-B_x5_5000.png", intensity_range=(0.2, 0.6), min_area=200)
imshow(contour)