# A module to support transport calculations
# Nancy Soontiens May 2022

from geopy.distance import distance
import numpy as np
import xarray as xr


def depth_integrate(var, mask, e3t):
    """Return depth integrated variable
    """
    return (var*mask*e3t).sum(dim='depth', skipna=True)


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
    return dot(vector, ref_vector)/dot(ref_vector, ref_vector)*ref_vector


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


def get_rotated_velocities(u_transect, v_transect, transect):
    """Rotate zonal and meridonal velocites u_transect and v_transect to
    across and along transect components"""
    all_dx = np.sum(np.real(transect))
    all_dy = np.sum(np.imag(transect))
    ref_vector = all_dx + 1j*all_dy
    
    vels = u_transect + 1j*v_transect
    rotated_vels = change_axis(vels, ref_vector)
    # Isolate along and across components
    across_vel = np.imag(rotated_vels)
    across_vel.attrs = {'units': 'm/s',
                        'short_name': 'v_across',
                        'long_name': 'Across-transect Velocity'}
    along_vel = np.real(rotated_vels)
    along_vel.attrs = {'units': 'm/s',
                        'short_name': 'v_along',
                        'long_name': 'Along-transect Velocity'}
    return along_vel, across_vel


def calculate_transport(across_vel,
                        transect,
                        mask_transect,
                        e3t_transect,
                        depth_transect,
                        barotropic=True
                        ):
    """Caclulate transport along a transect.
    across_vel is the across tranport velocity interpolated to the 
    the transect coordinates.
    transect is an array of complex numbers dx +j*dy where dx and dy represent the
    x displacement and y displacment for each transect segment.
    By default, across_vel should be depth-averaged.
    For depth-dependent velocities, pass baratopric=False.
   """
    n = len(transect)
    transport = 0
    for i in range(n):
        ref_vector = transect[i]
        across = across_vel.isel(transect=i)
        if not barotropic:
            depth_integrated_vel = depth_integrate(across, 
                                                   mask_transect.isel(transect=i), 
                                                   e3t_transect.isel(transect=i))
        else:
            H = depth_transect.isel(transect=i) 
            depth_integrated_vel = across*H
            
        tadd = depth_integrated_vel*np.abs(ref_vector)
        tadd = tadd.where(~np.isnan(tadd), 0)
        transport += tadd
    transport.attrs = {'units': 'm^3/s',
                       'short_name': 'transport',
                       'long_name': 'Volume Transport'}
    return transport


def get_transect_transport(u, v, mask, e3t, depth, transect,
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
        H = depth.where(depth!=0)
        u = depth_integrate(u, mask, e3t)/H
        v = depth_integrate(v, mask, e3t)/H
 
    mask_transect = interpolate_transect(mask, lons, lats, method='nearest')
    e3t_transect = interpolate_transect(e3t, lons, lats)
    depth_transect = interpolate_transect(depth, lons, lats)

    u_transect = interpolate_transect(u, lons, lats)
    v_transect = interpolate_transect(v, lons, lats)

    along_vel, across_vel = get_rotated_velocities(u_transect,
                                                   v_transect,
                                                   transect_vector)
    transport = calculate_transport(across_vel,
                                    transect_vector,
                                    mask_transect,
                                    e3t_transect,
                                    depth_transect,
                                    barotropic=barotropic
                                    )

    return transport, along_vel, across_vel
