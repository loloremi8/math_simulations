import os
from PIL import Image

# Directory containing the frames
frames_dir = "frames"
output_filename = "wave_propagation_3D.gif"

# Get a sorted list of frame filenames
frame_files = sorted([f for f in os.listdir(frames_dir) if f.endswith('.png')])

# Load the frames into a list
frames = [Image.open(os.path.join(frames_dir, frame)) for frame in frame_files]

# Save the frames as a GIF
frames[0].save(output_filename, save_all=True, append_images=frames[1:], duration=50, loop=0)
print(f"GIF saved as {output_filename}")