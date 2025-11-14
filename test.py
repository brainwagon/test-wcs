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
arrow_length = 20
angle = 0 - solution["Roll"]
angle_rad = math.radians(angle)
end_x = center_x + arrow_length * math.sin(angle_rad)
end_y = center_y - arrow_length * math.cos(angle_rad)
draw.line([center_x, center_y, end_x, end_y], fill="red", width=3)

img.save("overlay.jpg")
