using Images
using FileIO

function process_image(input_path::String, output_path::String)
    # Load the image
    img = load(input_path)
    
    # Convert to RGB array
    img_array = channelview(img)
    
    # Create output array with same size
    height, width = size(img)[1:2]
    output        = similar(img)
    soft_red      = RGB(0.85, 0.45, 0.45)    # Soft pastel red
    soft_green    = RGB(0.45, 0.85, 0.45)  # Soft pastel green
    soft_blue     = RGB(0.5, 0.7, 0.9)     # Soft pastel blue for high contrast

    # Process each pixel
    for i in 1:height
        for j in 1:width
            # Get RGB values (normalized between 0 and 1)
            r = Float64(img_array[1,i,j])
            g = Float64(img_array[2,i,j])
            b = Float64(img_array[3,i,j])
            
            # Calculate intensity (simple average method)
            intensity = (r + g + b) / 3
            
            # Threshold for high/low contrast (you can adjust this value)
            threshold = 0.2
        
            
            if intensity > threshold
                # High intensity -> Green
                output[i,j] = soft_blue
            else
                # Low intensity -> Red
                output[i,j] = soft_red
            end
        end
    end
    
    # Save the processed image
    save(output_path, output)
    return output
end

function process_image_with_split(input_path::String, output_path::String, split_x::Int,threshold::Float64)
    # Load the image
    img = load(input_path)
    
    # Convert to RGB array
    img_array = channelview(img)
    
    # Create output array with same size
    height, width = size(img)[1:2]
    output = similar(img)
    
    # Define soft colors
    soft_blue = RGB(0.5, 0.7, 0.9)     # For one zone (e.g., left side)
    soft_green = RGB(0.45, 0.85, 0.45) # For other zone (e.g., right side)
    soft_red      = RGB(0.85, 0.45, 0.45) 
    # Process each pixel
    for i in 1:height
        for j in 1:width
            # Get RGB values (normalized between 0 and 1)
            r = Float64(img_array[1,i,j])
            g = Float64(img_array[2,i,j])
            b = Float64(img_array[3,i,j])
            
            # Calculate intensity (simple average method)
            intensity = (r + g + b) / 3
            
            # Threshold for high/low contrast
            
            
            if intensity > threshold
                # High intensity areas split by x-coordinate
                if j < split_x
                    output[i,j] = soft_red   # Left zone
                else
                    output[i,j] = soft_green  # Right zone
                end
            else
                # Low intensity -> soft green (or could be another color)
                output[i,j] = soft_blue
            end
        end
    end
    
    # Save the processed image
    save(output_path, output)
    return output
end

# Example usage

split_position= 550
threshold = 0.3
input_file = "/Users/gianpaulrincon/Documents/GitHub/Rini/entwicklung/ImageDelamination/20250304-B_x5_1000.png"   # Replace with your input file path
output_file = "2050304-B_1000_5x.png"  # Output file path
process_image_with_split(input_file, output_file, split_position,threshold)

split_position= 310
threshold=0.3
input_file = "/Users/gianpaulrincon/Documents/GitHub/Rini/entwicklung/ImageDelamination/20250304-B_x5_3000.png"   # Replace with your input file path
output_file = "2050304-B_3000_5x.png"  # Output file path
process_image_with_split(input_file, output_file, split_position,threshold)

split_position= 360
threshold=0.3
input_file = "/Users/gianpaulrincon/Documents/GitHub/Rini/entwicklung/ImageDelamination/20250304-B_x5_5000.png"   # Replace with your input file path
output_file = "2050304-B_5000_5x.png"  # Output file path
process_image_with_split(input_file, output_file, split_position,threshold)

split_position= 360
threshold=0.3
nput_file = "/Users/gianpaulrincon/Documents/GitHub/Rini/entwicklung/ImageDelamination/20250304-B_x5_55000.png"   # Replace with your input file path
output_file = "2050304-B_55000_5x.png"  # Output file path
process_image_with_split(input_file, output_file, split_position,threshold)




