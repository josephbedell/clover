#!/usr/bin/env python3
from PIL import Image

# Open the image file
with Image.open('four_leaf_clover.jpeg') as img:
    # Save the image in a different format (e.g., PNG)
    img.save('four_leaf_clover.png', 'PNG')