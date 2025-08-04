import math 
from pyproj import CRS, Transformer
import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline
import matplotlib.pyplot as plt

#lon for longitude, lat for latitude, r for radius in kilometer
""" pinner() function calculates the tangent points of two circles on a transverse mercator projection.
    It takes the coordinates and radii of two circles and returns the tangent points in longitude and latitude.
    The function uses the pyproj library to handle coordinate transformations and the math library for calculations.
    The tangent points are calculated in a custom transverse mercator projection centered between the two circles.
    The function returns a list containing the coordinates of the tangent points in longitude and latitude.
    The tangent points are calculated using the properties of circles and the geometry of the transverse mercator projection.
    The function is designed to work with geographic coordinates and is useful for applications that require precise calculations
    of tangent points between circles on a map."""
def pinner(lon1, lat1, r1, lon2, lat2, r2):

    #Defining the central meridian and latitude for custom transverse merator projection
    #Method: Finding the middle point between 2 points

    center_lon = (lon1 + lon2) / 2
    center_lat = (lat1 + lat2) / 2
    l1 = r1 * 1000
    l2 = r2 * 1000

    #Defining the custom transverse merator(TMERC) projection around the center_lon and center_lat
    #Method:
    #+proj=tmerc is for choosing the TMERC projection
    #+lat_0 and +lon_0 is for choosing the center of the projection
    #+k_0 is the scaling factor of the projection, 1.0 means no scaling
    #+x_0 and +y_0 sets false easting/northing, for preventing negative values, but not required for the situation, so both 0
    #+ellps = WGS84 specifies the earth shape model used for projection, which is an ellipsoid
    #+units = m sets the units to be default aka meters
    #+no_defs means no default setting is used in sense of the normal projection -- all i know it is a relic of old versions of pyproj

    custom_tmerc_proj_string = (
    f"+proj=tmerc +lat_0={center_lat} +lon_0={center_lon} +k_0=1.0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    )

    #creating a TMERC projection based on the parameters above

    custom_tmerc_crs = CRS(custom_tmerc_proj_string)

    #Standard WGS84 Geographic CRS -- the standard latitude and longitude system for the input of the function

    wgs84_crs = CRS("EPSG:4326")

    #Building the transformers between the standard WGS84 and the custom TMERC projection
    #The first two inputs indicates the different projections
    #always_xy makes sure all systems work with x being horizontal(West-East) and y being vertial(North-South)

    transformer_to_tmerc = Transformer.from_crs(wgs84_crs, custom_tmerc_crs, always_xy=True)
    transformer_from_tmerc = Transformer.from_crs(custom_tmerc_crs, wgs84_crs, always_xy=True)

    ###Checker for the projection conditions
    #print(f"Custom T-Merc centered at ({center_lon:.4f}, {center_lat:.4f}) with k_0=1.0")

    #Convert Lon/Lat points to custom projected (x, y) on TMERC projection
    x1, y1 = transformer_to_tmerc.transform(lon1, lat1)
    x2, y2 = transformer_to_tmerc.transform(lon2, lat2)

    # * Value checker of inputs:
    #print(f"\nCircle 1 center in T-Merc: ({x1:.2f}, {y1:.2f}) meters, Radius: {l1} m")
    #print(f"Circle 2 center in T-Merc: ({x2:.2f}, {y2:.2f}) meters, Radius: {l2} m")

    #The calculations of the new coordinates

    #Find common external center coordinates

    xf = (x1 * l2 - x2 * l1)/(l2 - l1)
    yf = (y1 * l2 - y2 * l1)/(l2 - l1)

    #Find the gradients -- (mP and mQ) of the two common tangents, while A, B, C and d are variables for calculations
    A = (x2 - xf) ** 2 - l2 ** 2
    B = -2 * (x2 - xf) * (y2 - yf)
    C = (y2 - yf) ** 2 - l2 ** 2
    d = B ** 2 - 4 * A * C
    mP = ( -B + math.sqrt(d) )/ (2 * A)
    mQ = ( -B - math.sqrt(d) )/ (2 * A)

    # * Value Checker for the 2D gradients
    #print(f"{mP:.3f}")
    #print(f"{mQ:.3f}")

    #Computing the coordinates -- xP1, yP1, xP2, yP2, xQ1, yQ1, xQ2, yQ2
    #S1, S2, T1, T2, T3, T4 are variables for calculations

    S1 = 1 + mP ** 2
    S2 = 1 + mQ ** 2
    T1 = -2 * x1 + 2 * mP * (yf - mP * xf - y1)
    T2 = -2 * x2 + 2 * mP * (yf - mP * xf - y2)
    T3 = -2 * x1 + 2 * mQ * (yf - mQ * xf - y1)
    T4 = -2 * x2 + 2 * mQ * (yf - mQ * xf - y2)

    xP1 = -T1 / (2 * S1)
    yP1 = mP * xP1 + (yf - mP * xf)
    xP2 = -T2 / (2 * S1)
    yP2 = mP * xP2 + (yf - mP * xf)
    xQ1 = -T3 / (2 * S2)
    yQ1 = mQ * xQ1 + (yf - mQ * xf)
    xQ2 = -T4 / (2 * S2)
    yQ2 = mQ * xQ2 + (yf - mQ * xf)

    tangent_points_in_xy = [
                            (xP1, yP1), (xP2, yP2), (xQ1, yQ1), (xQ2, yQ2)
                            ]
    
    # * Value Checker of  (x, y) coordinates in 2D plane
    #for i,(x, y) in enumerate(tangent_points_in_xy):
        #print(f"Coord Point{i + 1}: ({x:.2f},{y:.2f})")

    #Act as the return of the function for later smoothing steps

    coord_returner_P = [[], []]
    coord_returner_Q = [[], []]

    #Convert back to lon/lat

    for i, (tx, ty) in enumerate(tangent_points_in_xy):
        lon_t, lat_t = transformer_from_tmerc.transform(tx, ty)
        #Separate lon and lats into 2 different lists for later interpolation uses

        if i < 2:

            # * Value-checker: print(f"Point P{i+1}: ({lon_t:.2f}, {lat_t:.2f})")

            coord_returner_P[0].append(lon_t)
            coord_returner_P[1].append(lat_t)
        else:

            # * Value_checker: print(f"Point Q{i-1}: ({lon_t:.2f}, {lat_t:.2f})")

            coord_returner_Q[0].append(lon_t)
            coord_returner_Q[1].append(lat_t)

    #Returns an array that contains 4 arrays: [[P-line lon], [P-line lat], [Q-line lon], [Q-line lat]]

    return [coord_returner_P[0], coord_returner_P[1], coord_returner_Q[0], coord_returner_Q[1]]

###Testing for the function:pinner()
#pinner(108.8, 19, 100, 110.4, 22.1, 150)




#LoP represents the list of points and their corresponding radii
#generator() function takes a list of points and their radii, and generates the tangent points for each pair of consecutive points.
#It uses the pinner() function to calculate the tangent points and then averages them to create smooth curves for the P and Q paths.
#The function returns a list containing the x and y coordinates of the P and Q paths.
#It also plots the P and Q paths using matplotlib for visualization.
#The function is designed to work with geographic coordinates and is useful for applications that require precise calculations
#of tangent points between circles on a map.
#It is useful for applications that require precise calculations of tangent points between circles on a map. 

def generator(LoP: list):
    #Put different Plon, Plat, Qlon, Qlat together for a array to be separated eventually to interpolate a P-curve and Q-curve

    logger = [[], [], [], []]
    for i in range(len(LoP) - 1):
        holder = pinner(LoP[i][0], LoP[i][1], LoP[i][2], LoP[i+1][0], LoP[i+1][1], LoP[i+1][2])
        for j in range(4):
            logger[j] = logger[j] + holder[j]

    # * Value Test: print(logger)

    #Averaging out 2 tangential points on same circle for better curve making
    averager = [[], [], [], []]
    for i in range(len(averager)):
        averager[i].append(logger[i][0])
        for j in range(1, len(logger[0]) - 1, 2):
            average_point = (logger[i][j] + logger[i][j+1]) / 2
            averager[i].append(average_point)
        averager[i].append(logger[i][-1])       

    # * Value Test: print(averager)

    P_and_Q = []

    #P_and_Q contains 2 arrays, every index is the array of size 4
    #P_and_Q = [
    #           [[x-path of P], [x-curve of P], [y-path of P] [y-curve of P]], 
    #           [[x-path of Q], [x-curve of Q], [y-path of P],[y-curve of Q]]
    #           ] 
    for i in range(0, len(averager), 2):
        x_path = np.array(averager[i])
        y_path = np.array(averager[i+1])

        #Parametrization
        t_dots = np.arange(len(x_path))

        #Creating t to x/y spline, k=3 for the spline function to be cubic
        spl_x = make_interp_spline(t_dots, x_path, k=3)
        spl_y = make_interp_spline(t_dots, y_path, k=3)

        #linspace maker
        t_smooth = np.linspace(t_dots.min(), t_dots.max(), 200)

        #Making a x_curve and y_curve based on linspace of t
        x_curve = spl_x(t_smooth)
        y_curve = spl_y(t_smooth)

        P_and_Q.append([x_path, x_curve, y_path ,y_curve])
    
    #Pyplot parameters
    plt.figure(figsize=(10, 6))
    plt.plot(P_and_Q[0][1], P_and_Q[0][3], '-', color='red', linewidth=2, label='P-path')
    plt.plot(P_and_Q[1][1], P_and_Q[1][3], '-', color='blue', linewidth=2, label='Q-path')
    plt.title('Parametric Spline for Non-Monotonic X-Values')
    plt.xlabel('X-coordinate')
    plt.ylabel('Y-coordinate')
    plt.axis('equal') # Important for paths to maintain aspect ratio
    plt.legend()
    plt.grid(True)
    plt.show()
       

test_list = [[108.8, 19, 100], [110.4, 22.1, 150], [115.7, 26.1, 270], [120.6, 28.8, 400]]
generator(test_list)

