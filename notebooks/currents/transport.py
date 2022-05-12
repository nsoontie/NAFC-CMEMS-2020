# A module to support transport calculations
# Nancy Soontiens May 2022

from geopy.distance import distance
import numpy as np
import xarray as xr


def depth_integrate(var, mask, e3t, depth_axis=0):
    """Return depth integrated variable
    """
    return np.sum(var*mask*e3t, axis=depth_axis)


def interpolate_transect(var, lons, lats, method='linear'):
    """Interpolate a variable along a transect defined by lons, lats"""
    lons = xr.DataArray(lons, dims='transect')
    lats = xr.DataArray(lats, dims='transect')
    return var.interp(longitude=lons, latitude=lats, method=method)


def dot(a,b):
    """dot product between vectors a, b where a and b are complex"""
    return np.real(a)*np.real(b) + np.imag(a)*np.imag(b)


def proj_vector(vector, ref_vector):
    """project vector onto ref_vector, where vector and ref_vector are 
    complex numbers"""
    return (dot(vector, ref_vector)/dot(ref_vector, ref_vector)*ref_vector)


def proj_vector_perp(vector, ref_vector):
    """project vector onto the vector rotated 90 deg CCW from ref_vector"""
    return proj_vector(vector, 1j*ref_vector)


def change_axis(vector, ref_vector):
    """Rotate vector to along and across components of ref_vector"""
    along = np.sign(dot(vector, ref_vector)) * \
            np.abs(proj_vector(vector, ref_vector))
    across = np.sign(dot(vector, 1j*ref_vector)) * \
             np.abs(proj_vector_perp(vector, ref_vector))
    return along + 1j*across


def vectorize_transect(lons, lats, dir_east_perp, dir_north_perp):
    """Translate lon and lat coordinates of a transect into dx and dy 
    distance in metres.  The orientation of the vector is assigned based
    on dir_east_perp and dir_north_perp which take on values of -1 or 1 and 
    represent the east and north direction of the vector that is perpendicular
    to the transect. For example, dir_east_perp=1 and dir_north_perp=-1 means the 
    vector perpendicular to this transect is in the 4th quadrant"""
    across_vector = dir_east_perp + 1j*dir_north_perp
    along_vector  = -1j*across_vector # rotate 90 CW to get along vector
    dir_east_along = np.real(along_vector)
    dir_north_along = np.imag(along_vector)
    n = len(lons) - 1
    vector = np.empty((n,),dtype=complex)
    for i in range(n):
        dx = distance((lons[i], lats[i]), (lons[i+1], lats[i])).m
        dy = distance((lons[i], lats[i]), (lons[i], lats[i+1])).m
        vector[i] = dir_east_along*dx + dir_north_along*1j*dy
    return vector


def calculate_transport(u_transect, v_transect,
                        transect,
                        barotropic=True,
                        mask_transect=None, e3t_transect=None,
                        depth_axis=0):
    """Caclulate transport along a transect.
    u_transect, v_transect are the zonal and meridional velocities interpolated to 
    the transect coordinates.
    transect is an array of complex numbers dx +j*dy where dx and dy represent the
    x displacement and y displacment for each transect segment.
    By default, u_transect and v_transect should be depth-averaged.
    For depth-dependent velocities, pass baratopric=False and the land mask and 
    vertical grid spacing interpolated to the transect.
   """
    n = len(transect)
    transport = 0
    for i in range(n):
        vel = u_transect.isel(transect=i) +1j*v_transect.isel(transect=i)
        ref_vector = transect[i]
        rotated_vel = change_axis(vel, ref_vector)
        if not barotropic:
            rotated_vel = depth_integrate(rotated_vel, 
                                          mask_transect.isel(transect=i), 
                                          e3t_transect.isel(transect=i),
                                          depth_axis=depth_axis) 
        transport += np.imag(rotated_vel)*np.abs(ref_vector)
    return transport


def get_transect_transport(u, v, mask, e3t, transect,
                           depth_axis=0,
                           barotropic=True,
                           num_points=100):
    lons = np.linspace(transect.lon1.values[0],
                       transect.lon2.values[0],
                       num=num_points)
    lats = np.linspace(transect.lat1.values[0],
                       transect.lat2.values[0],
                       num=num_points)
    dir_east_perp = transect['Dir-E'].values[0]
    dir_north_perp = transect['Dir-N'].values[0]
    transect_vector = vectorize_transect(lons, lats,
                                         dir_east_perp,
                                         dir_north_perp)

    if barotropic:
        u = depth_integrate(u, mask, e3t, depth_axis=depth_axis)
        v = depth_integrate(v, mask, e3t, depth_axis=depth_axis)
        mask_transect = None
        e3t_transect = None
    else:
        mask_transect = interpolate_transect(mask, lons, lats, method='nearest')
        e3t_transect = interpolate_transect(e3t, lons, lats)

    u_transect = interpolate_transect(u, lons, lats)
    v_transect = interpolate_transect(v, lons, lats)

    transport = calculate_transport(u_transect, v_transect,
                                    transect_vector,
                                    barotropic=barotropic,
                                    mask_transect=mask_transect,
                                    e3t_transect=e3t_transect,
                                    depth_axis=depth_axis)

    return transport
