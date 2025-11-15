import ephem
from PIL import Image, ImageDraw
import math
from math import floor, ceil, radians
import numpy as np
np.math = math
import tetra3
import astropy
import time

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import fit_wcs_from_points

t3 = tetra3.Tetra3()

start_time = time.perf_counter()
img = Image.open("stars-0000.jpg")
w, h = img.size
print(f"Size: {w}x{h}")
end_time = time.perf_counter()

elapsed = end_time - start_time
print(f"Elapsed time: {elapsed:.6f} seconds")

start_time = time.perf_counter()
solution = t3.solve_from_image(img, return_matches=True)
end_time = time.perf_counter()
elapsed = end_time - start_time
print(f"Solution time: {elapsed:.6f} seconds")

print(f'RA = {solution["RA"]:.4f} {ephem.hours(radians(solution["RA"]))}')
print(f'Dec = {solution["Dec"]:.4f} {ephem.degrees(radians(solution["Dec"]))}')
print(f'Roll = {solution["Roll"]:.4f}')
cra = solution["RA"]
cdec = solution["Dec"]

matched_stars = np.array(solution["matched_stars"])
matched_centroids = np.array(solution["matched_centroids"])

star_x = matched_centroids[:,1]
star_y = matched_centroids[:,0]
star_xy = (star_x, star_y)

star_ra = np.array(matched_stars[:,0]) * u.deg
star_dec = np.array(matched_stars[:,1]) * u.deg

center = SkyCoord(solution["RA"] * u.deg, solution["Dec"] * u.deg)
world_coords = SkyCoord(ra=star_ra, dec=star_dec, frame='icrs')

start_time = time.perf_counter()
wcs = fit_wcs_from_points(
        star_xy,
        world_coords,
        proj_point=center,
        projection='TAN',
        sip_degree=2)
end_time = time.perf_counter()
elapsed = end_time - start_time
print(f"Fit time: {elapsed:.6f} seconds")

ra, dec = wcs.all_pix2world(w / 2., h / 2., 0)

print(f"Center RA  = {ra:.4f}")
print(f"Center Dec = {dec:.4f}")

n = len(matched_stars)
print(f"We matched {n} stars.")

draw = ImageDraw.Draw(img)
for i in range(n):
    ra, dec, _ = matched_stars[i]
    y, x = matched_centroids[i]
    nx, ny = wcs.all_world2pix(ra, dec, 0, ra_dec_order=True)
    ex = abs(x-nx)
    ey = abs(y-ny)
    print(ra, dec, x, y, ex, ey)

for i in range(floor(cra-6), ceil(cra+6), 2):
    for j in range(floor(cdec-6), ceil(cdec+6), 2):
        try:
            ra = i * u.deg
            dec0 = j  * u.deg
            dec1 = (j+2) * u.deg
            x0, y0 = wcs.all_world2pix(ra, dec0, 0, ra_dec_order=True)
            x1, y1 = wcs.all_world2pix(ra, dec1, 0, ra_dec_order=True)
            draw.line([x0, y0, x1, y1], fill=(255, 255, 255))
        except:
            pass

for j in range(floor(cdec-6), ceil(cdec+6), 2):
    for i in range(floor(cra-6), ceil(cra+6), 2):
        try:
            ra0 = i * u.deg
            ra1 = (i+2) * u.deg
            dec = j  * u.deg
            x0, y0 = wcs.all_world2pix(ra0, dec, 0, ra_dec_order=True)
            x1, y1 = wcs.all_world2pix(ra1, dec, 0, ra_dec_order=True)
            draw.line([x0, y0, x1, y1], fill=(255, 255, 255))
        except:
            pass

center_x = w / 2
center_y = h / 2
arrow_length = 40
angle = 0 - solution["Roll"]
angle_rad = math.radians(angle)
end_x = center_x + arrow_length * math.sin(angle_rad)
end_y = center_y - arrow_length * math.cos(angle_rad)
draw.line([center_x, center_y, end_x, end_y], fill="red", width=2)

# Draw the arrowhead (pointed V)
arrowhead_length = 10
arrowhead_angle = math.radians(30)

# First barb
barb1_x = end_x - arrowhead_length * math.sin(angle_rad + arrowhead_angle)
barb1_y = end_y + arrowhead_length * math.cos(angle_rad + arrowhead_angle)
draw.line([end_x, end_y, barb1_x, barb1_y], fill="red", width=2)

# Second barb
barb2_x = end_x - arrowhead_length * math.sin(angle_rad - arrowhead_angle)
barb2_y = end_y + arrowhead_length * math.cos(angle_rad - arrowhead_angle)
draw.line([end_x, end_y, barb2_x, barb2_y], fill="red", width=2)

# Add the label "N" at the tip
font_size = 15
try:
    from PIL import ImageFont
    font = ImageFont.truetype("arial.ttf", font_size)
except IOError:
    font = ImageFont.load_default()

# Calculate text position with a small offset from the arrow tip
text_offset = 10 # Increased offset to avoid overlap with arrowhead
text_x = end_x + text_offset * math.sin(angle_rad)
text_y = end_y - text_offset * math.cos(angle_rad)

# Adjust text position to center it roughly
# Note: draw.textsize is deprecated, using font.getbbox for modern Pillow
bbox = font.getbbox("N")
text_width = bbox[2] - bbox[0]
text_height = bbox[3] - bbox[1]

text_x -= text_width / 2
text_y -= text_height / 2

draw.text((text_x, text_y), "N", fill="white", font=font)

# Second arrow (E)
angle_e = angle - 90
angle_rad_e = math.radians(angle_e)
end_x_e = center_x + arrow_length * math.sin(angle_rad_e)
end_y_e = center_y - arrow_length * math.cos(angle_rad_e)
draw.line([center_x, center_y, end_x_e, end_y_e], fill="red", width=2)

# Arrowhead for the second arrow
barb1_x_e = end_x_e - arrowhead_length * math.sin(angle_rad_e + arrowhead_angle)
barb1_y_e = end_y_e + arrowhead_length * math.cos(angle_rad_e + arrowhead_angle)
draw.line([end_x_e, end_y_e, barb1_x_e, barb1_y_e], fill="red", width=2)

barb2_x_e = end_x_e - arrowhead_length * math.sin(angle_rad_e - arrowhead_angle)
barb2_y_e = end_y_e + arrowhead_length * math.cos(angle_rad_e - arrowhead_angle)
draw.line([end_x_e, end_y_e, barb2_x_e, barb2_y_e], fill="red", width=2)

# Label for the second arrow
text_x_e = end_x_e + text_offset * math.sin(angle_rad_e)
text_y_e = end_y_e - text_offset * math.cos(angle_rad_e)

bbox_e = font.getbbox("E")
text_width_e = bbox_e[2] - bbox_e[0]
text_height_e = bbox_e[3] - bbox_e[1]

text_x_e -= text_width_e / 2
text_y_e -= text_height_e / 2

draw.text((text_x_e, text_y_e), "E", fill="white", font=font)

img.save("overlay.jpg")
