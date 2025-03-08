using Images
using FileIO

function process_image(input_path::String, output_path::String)
    # Load the image
    img = load(input_path)
    
    # Convert to RGB array
    img_array = channelview(img)
    
    # Create output array with same size
    height, width = size(img)[1:2]
    output = similar(img)
    soft_red = RGB(0.85, 0.45, 0.45)    # Soft pastel red
    soft_green = RGB(0.45, 0.85, 0.45)  # Soft pastel green
    soft_blue = RGB(0.5, 0.7, 0.9)     # Soft pastel blue for high contrast

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

# Example usage
input_file = "/Users/gianpaulrincon/Documents/GitHub/Rini/entwicklung/ImageDelamination/20250307-D_Start_10x.png"   # Replace with your input file path
output_file = "2050307_D_Start.png"  # Output file path
process_image(input_file, output_file)


input_file = "/Users/gianpaulrincon/Documents/GitHub/Rini/entwicklung/ImageDelamination/20250307-D_10000_10x.png"   # Replace with your input file path
output_file = "2050307_D_10000.png"  # Output file path
process_image(input_file, output_file)
