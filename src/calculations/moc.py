import calculations.init as init
import calculations.nozzle as nozzle
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

# ------------ Constants ------------------
gamma = init.GAMMA
n_lines = init.N_CHARACTERISTICS
mach_exit = nozzle.mach_design_exit(init.AE_AT)
r_throat = init.R_THROAT
pm_a = init.PRANDTL_MEYER_A
theta_a = init.THETA_A

# ------------ Functions ------------------
def mach_angle(mach):
    """
    Calculate the Mach angle (in degrees) for a given Mach number.
    """
    if mach <= 1:
        raise ValueError("Mach number must be greater than 1.")
    return np.degrees(np.asin(1 / mach))


def prandtl_meyer_angle(mach):
    """
    Calculate the Prandtl-Meyer angle (in degrees) for a given Mach number.
    """
    if mach < 1:
        raise ValueError("Mach number must be greater than or equal to 1.")
    term1 = np.sqrt((gamma + 1) / (gamma - 1))
    term2 = np.atan(np.sqrt((gamma - 1) * (mach**2 - 1) / (gamma + 1)))
    term3 = np.atan(np.sqrt(mach**2 - 1))
    return float(np.degrees(term1 * term2 - term3))


def mach_from_prandtl_meyer_angle(angle):
    """
    Calculate Mach number from the Prandtl-Meyer angle (in degrees).
    """
    function = lambda mach: prandtl_meyer_angle(mach) - angle
    mach = fsolve(function, 2.0)[0]
    return mach

def characteristic_points(n_lines):
    """
    Calculate the characteristic points for the Method of Characteristics.
    """
    # global points
    # points = {
    #     "theta"     # flow angle (theta) in degrees,
    #     "pm"        # prandtl_meyer_angle in degrees,
    #     "mach"      # mach number
    #     "mu"        # mach angle
    #     "x"         # x location    
    #     "y"         # y location
    #     "isw"       #  is at the wall
    #     "iscl"      # is at the centerline
    # }

    charact_points = {}

    global n_points
    n_points = n_lines + n_lines * (n_lines + 1) / 2
    j = 1 + n_lines
    k = 0

    for i in range(1, int(n_points) + 1):
        point = {
            "index": i,        # index of the point
            "theta": None,      # flow angle (theta) in degrees
            "pm": None,         # prandtl-meyer angle in degrees
            "mach": None,       # mach number
            "mu": None,         # mach angle
            "x": None,          # x location
            "y": None,          # y location
            "isw": False,       # is at the wall
            "iscl": False       # is at the centerline
        }

        # Create characteristic points
        if i == j + k:
            point["isw"] = True
            k += 1
            j += n_lines - k
            
        charact_points[i] = point

    # Set centerline points accordingly (theta = 0 and y = 0)
    for i in range(1, int(n_points) + 1):
        if charact_points[i]["index"] == 1:
            charact_points[i]["iscl"] = True
            charact_points[i]["theta"] = 0
            charact_points[i]["y"] = 0
        elif i > 1:
            if charact_points[i - 1]["isw"] == True:
                charact_points[i]["iscl"] = True
                charact_points[i]["theta"] = 0
                charact_points[i]["y"] = 0
    return charact_points

         
def interpolate_angle(theta, n_lines):
    """
    Interpolate the angle theta for the Method of Characteristics.    
    """
    result = [0]
    temp = theta - (1 + (theta % 1))
    delta = temp / (n_lines - 2)
    for i in range(n_lines - 2):
        result.append(result[i] + delta)
    result.append(theta)
    return result

# def vector_to_string(vector):
#     """
#     Convert a vector to a string with 6 decimal places.
#     """
#     return " ".join([f"{x:.6f}" for x in vector])


# def vector_to_string_int(vector):
#     """
#     Convert a vector to a string with no decimal places.
#     """
#     return " ".join(map(str, vector))


def return_xy_intersection_point(xt, yt, tht_to, xb, yb, tht_b):
    """
    Calculate the intersection point of two lines defined by their angles and points.
    """
    tht_to_rad = np.radians(tht_to)
    tht_b_rad = np.radians(tht_b)

    tan_to = np.tan(tht_to_rad)
    tan_b = np.tan(tht_b_rad)

    x = (xt * tan_to - xb * tan_b + yb - yt) / (tan_to - tan_b)
    y = (tan_to * tan_b * (xt - xb) + tan_to * yb - tan_b * yt) / (tan_to - tan_b)

    return (x, y)





# ---------- Execute ------------------
def execute():

    print(f"\nMach Exit Design: {mach_exit:.5f}")

    pe_p0 = nozzle.static_pressure_ratio(mach_exit)
    print(f"Static Pressure Ratio at Exit (Pe_P0): {pe_p0:.5f}\n")

    theta_max = prandtl_meyer_angle(mach_exit) / 2
    print(f"Maximum wall angle (\u03B8_max): {theta_max:.5f} degrees\n")

    theta_array = interpolate_angle(theta_max, n_lines)
    print(f"Interpolated angles (\u03B8_array): {theta_array} degrees")
    charact_points = characteristic_points(n_lines)
    # print(f"Characteristic points: {charact_points}")
    wall_points = [i for i in charact_points if charact_points[i]["isw"]==True]
    print(f"\nWall points: {wall_points}")
    center_points = [i for i in charact_points if charact_points[i]["iscl"]==True]
    print(f"Center points: {center_points}")

    mach_a = mach_from_prandtl_meyer_angle(pm_a)
    print(f"\nMach number at point a: {mach_a:.5f}")
    mu_a = mach_angle(mach_a)
    print(f"Mach angle at point a: {mu_a:.5f} degrees")

    ### ==================
    ### BEGIN MOC
    ### ==================

    # Point a (throat)
    xa = 0                      # x location of point a (cm)
    ya = r_throat               # y location of point a (throat radius) (cm)
    theta_a                     # kick-off angle (deg)
    pm_a                        # same as kick-off angle (deg)
    mach_a                      # Mach number at point a

    # Point 1
    theta_1 = charact_points[1]["theta"]
    pm_1 = charact_points[1]["pm"] = theta_a + pm_a 
    mach_1 = charact_points[1]["mach"] = mach_from_prandtl_meyer_angle(pm_1)
    mu_1 = charact_points[1]["mu"] = mach_angle(mach_1)           # Mach angle at point 1 (deg)
    print(f"\nMach number at point 1: {mach_1:.5f}")
    print(f"Mach angle at point 1: {mu_1:.5f} degrees")
    
    # Slope at point a
    slope_a = ((theta_a - mu_a) + (theta_1 - mu_1)) / 2           # slope at point a (deg)

    # Point 1
    charact_points[1]["x"] = xa + r_throat * np.tan(np.radians(90 + slope_a))  # x location of point 1 (cm)
    pm_1
    mach_1 
    mu_1
    
    # Points 2 to n_lines + 1
    for i in range(2, n_lines + 2):
        # assign flow and pm angle (equal for first reflection)
        if charact_points[i]["isw"] == False:
            charact_points[i]["theta"] = theta_array[i-1]
            charact_points[i]["pm"] = theta_array[i-1]
            charact_points[i]["mach"] = mach_from_prandtl_meyer_angle(charact_points[i]["pm"])
            charact_points[i]["mu"] = mach_angle(charact_points[i]["mach"])

            # find theta_ax nu_ax and the Mach angle
            tv_ax = (2 * (theta_array[i-1]) + (charact_points[i-1]["theta"] - charact_points[i-1]["pm"])) / 2
            mach_ax = mach_from_prandtl_meyer_angle(tv_ax)
            mu_ax = mach_angle(mach_ax)

            # find left and right running characteristic slopes
            theta_b = (charact_points[i-1]["theta"] + charact_points[i-1]["mu"] + charact_points[i]["theta"] + charact_points[i]["mu"]) / 2
            theta_t = (tv_ax - mu_ax + charact_points[i]["theta"] - charact_points[i]["mu"]) / 2

            # find xy location of point i
            x_i, y_i = return_xy_intersection_point(xa, ya, theta_t, charact_points[i-1]["x"], charact_points[i-1]["y"], theta_b)
            charact_points[i]["x"] = x_i
            charact_points[i]["y"] = y_i

        # extend for the wall location
        else:
            charact_points[i]["theta"] = charact_points[i-1]["theta"]
            charact_points[i]["pm"] = charact_points[i-1]["pm"]
            charact_points[i]["mach"] = charact_points[i-1]["mach"]
            charact_points[i]["mu"] = charact_points[i-1]["mu"]

            # find left and second left running (both point up), the first left is the interpolated wall angle (max for first reflection)
            theta_t = theta_max
            theta_b = (charact_points[i-1]["theta"] + charact_points[i-1]["mu"] + charact_points[i]["theta"] + charact_points[i]["mu"]) / 2
            
            # find xy location of the wall point
            x_i, y_i = return_xy_intersection_point(xa, ya, theta_t, charact_points[i-1]["x"], charact_points[i-1]["y"], theta_b)
            charact_points[i]["x"] = x_i
            charact_points[i]["y"] = y_i
    
    # remaining points  
    j = 0
    k = 1
    for i in range(n_lines + 2, int(n_points) + 1):
        if charact_points[i]["iscl"] == True:
            charact_points[i]["pm"] = charact_points[i - (n_lines - j)]["theta"] + charact_points[i - (n_lines - j)]["pm"]
            charact_points[i]["mach"] = mach_from_prandtl_meyer_angle(charact_points[i]["pm"])
            charact_points[i]["mu"] = mach_angle(charact_points[i]["mach"])

            # find left and right running characteristic slopes
            theta_t = (charact_points[i - (n_lines - j)]["theta"] - charact_points[i - (n_lines - j)]["mu"] + charact_points[i]["theta"] - charact_points[i]["mu"]) / 2
            theta_b = 0

            # find xy location of centerline point
            x_i, y_i = return_xy_intersection_point(charact_points[i - (n_lines - j)]["x"], charact_points[i - (n_lines - j)]["y"], theta_t, charact_points[i- (n_lines - j + 1)]["x"], charact_points[i - (n_lines - j + 1)]["y"], theta_b)
            charact_points[i]["x"] = x_i
            charact_points[i]["y"] = y_i

        elif charact_points[i]["iscl"] == False and charact_points[i]["isw"] == False:
            # assign flow angle and pm angle
            charact_points[i]["theta"] = theta_array[k]
            charact_points[i]["pm"] = (charact_points[i - (n_lines - j)]["theta"] + charact_points[i - (n_lines - j)]["pm"] - (charact_points[i-1]["theta"] - charact_points[i-1]["pm"])) / 2
            
            # find the Mach number and Mach angle
            charact_points[i]["mach"] = mach_from_prandtl_meyer_angle(charact_points[i]["pm"])
            charact_points[i]["mu"] = mach_angle(charact_points[i]["mach"])

            # find left and right running characteristic slopes
            theta_t = (charact_points[i - (n_lines - j)]["theta"] - charact_points[i - (n_lines - j)]["mu"] + charact_points[i]["theta"] - charact_points[i]["mu"]) / 2
            theta_b = (charact_points[i-1]["theta"] + charact_points[i-1]["mu"] + charact_points[i]["theta"] + charact_points[i]["mu"]) / 2

            #find x and y location of point i
            x_i, y_i = return_xy_intersection_point(charact_points[i - (n_lines - j)]["x"], charact_points[i - (n_lines - j)]["y"], theta_t, charact_points[i-1]["x"], charact_points[i-1]["y"], theta_b)
            charact_points[i]["x"] = x_i
            charact_points[i]["y"] = y_i
            k = k + 1

        elif charact_points[i]["isw"] == True:
            charact_points[i]["theta"] = charact_points[i-1]["theta"]
            charact_points[i]["pm"] = charact_points[i-1]["pm"]
            charact_points[i]["mach"] = charact_points[i-1]["mach"]
            charact_points[i]["mu"] = charact_points[i-1]["mu"]
            
            # find left and right running characteristic slopes
            theta_t = charact_points[i-1]["theta"]
            theta_b = (charact_points[i-1]["theta"] + charact_points[i-1]["mu"] + charact_points[i]["theta"] + charact_points[i]["mu"]) / 2

            #find x and y location of point i
            x_i, y_i = return_xy_intersection_point(charact_points[i - (n_lines - j)]["x"], charact_points[i - (n_lines - j)]["y"], theta_t, charact_points[i-1]["x"], charact_points[i-1]["y"], theta_b)
            charact_points[i]["x"] = x_i
            charact_points[i]["y"] = y_i
            k = 1
            j = j + 1


    # Print the results
    print("\nCharacteristic Points:")
    for i in range(1, int(n_points) + 1):
        point = charact_points[i]
        print(f"Point {point['index']}: x = {point['x']:.5f}, y = {point['y']:.5f}, Mach = {point['mach']:.5f}, PM = {point['pm']:.5f}, Mu = {point['mu']:.5f}, Theta = {point['theta']:.5f} degrees")

    plt.figure(figsize=(10, 6))
    for i in range(1, int(n_points) + 1):
        point = charact_points[i]
        plt.plot(point['x'], point['y'], 'o')
        # plt.text(point['x'], point['y'], f"{point['index']}", fontsize=9, ha='right')
    plt.xlabel("X (cm)")
    plt.ylabel("Y (cm)")
    plt.title("Characteristic Points")
    plt.grid()

    
    x_plot = [charact_points[i]["x"] for i in range(1, int(n_points) + 1)]
    x_plot.insert(0, xa)  # Include the throat point
    y_plot = [charact_points[i]["y"] for i in range(1, int(n_points) + 1)]
    y_plot.insert(0, ya)  # Include the throat point

    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, y_plot, 'o')
    plt.xlabel("X (cm)")
    plt.ylabel("Y (cm)")
    plt.grid()



    # Generate grid for contour plot
    # Extract wall points for plotting the nozzle contour
    wall_x = [charact_points[i]["x"] for i in wall_points]
    wall_x.insert(0, xa)  # Include the throat point
    wall_y = [charact_points[i]["y"] for i in wall_points]
    wall_y.insert(0, ya)  # Include the throat point

    # Plot the nozzle contour
    plt.figure(figsize=(10, 6))
    plt.plot(wall_x, wall_y, label="Nozzle Wall Contour", color="blue")
    plt.xlabel("X Location (cm)")
    plt.ylabel("Y Location (cm)")
    plt.title("Nozzle Wall Contour")
    plt.grid()
    plt.legend()
    

    x_values = [charact_points[i]["x"] for i in range(1, int(n_points) + 1)]
    x_values.insert(0, xa)  # Include the throat point
    y_values = [charact_points[i]["y"] for i in range(1, int(n_points) + 1)]
    y_values.insert(0, ya)  # Include the throat point
    mach_values = [charact_points[i]["mach"] for i in range(1, int(n_points) + 1)]
    mach_values.insert(0, mach_a)  # Include the Mach number at the throat

    # Create a grid of x and y values
    x_grid = np.linspace(min(x_values), max(x_values), 100)
    y_grid = np.linspace(min(y_values), max(y_values), 100)
    X, Y = np.meshgrid(x_grid, y_grid)

    # Interpolate Mach values onto the grid
    from scipy.interpolate import griddata
    MACH = griddata((x_values, y_values), mach_values, (X, Y), method='linear')

    # Plot the contour
    plt.figure(figsize=(10, 6))
    contour = plt.contourf(X, Y, MACH, levels=50, cmap='viridis')
    plt.colorbar(contour, label="Mach Number")
    plt.xlabel("X Location (cm)")
    plt.ylabel("Y Location (cm)")
    plt.title("Mach Number Contour")
    plt.grid()

    
    # Output the results to a file
    with open("./data/nozzle_contour.txt", "w") as f:
        for i in range(0, len(wall_x)):
            x = wall_x[i]
            y = np.sqrt(wall_y[i])
            z = 0 
            f.write(f"{x:.5f} {y:.5f} {z}\n")

    plt.show()