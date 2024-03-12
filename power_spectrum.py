import numpy as np


def fft_radial_distance(shape, is_3d=False):
    print("shape:", shape, "is_3d:", is_3d)
    if is_3d:
        kz, ky, kx = np.meshgrid(np.fft.fftfreq(shape), np.fft.fftfreq(shape), np.fft.rfftfreq(shape)) #fftn uses last dimension as axis that is reduced
        k = np.sqrt(kx**2 + ky**2 + kz**2)*shape  # Calculate radial distance in Fourier space
    else:
        kx, ky = np.meshgrid(np.fft.rfftfreq(shape), np.fft.fftfreq(shape)) #fft2 uses first dimension as axis that is reduced
        k = np.sqrt(kx**2 + ky**2)*shape  # Calculate radial distance in Fourier space
    return k

 
def radial_bin_average(image_fft, radial_dist):
    # Calculate the magnitude of the FFT (power spectrum)
    magnitude = np.abs(image_fft)**2

    # Bin the values based on radial distance
    radial_dist_int = np.rint(radial_dist).astype(int)
    
    #print(radial_dist_int.shape, magnitude.shape, image_fft.shape)
    tbin = np.bincount(radial_dist_int.ravel(), magnitude.ravel())
    nr = np.bincount(radial_dist_int.ravel())


    radial_profile = tbin / nr
    return radial_profile, nr

def radial_bin_average_cross(image_fft1, image_fft2, radial_dist):
    # Calculate the magnitude of the FFT (power spectrum)
    magnitude = np.abs(image_fft1*np.conj(image_fft2))

    # Bin the values based on radial distance
    radial_dist_int = np.rint(radial_dist).astype(int)
    
    tbin = np.bincount(radial_dist_int.ravel(), magnitude.ravel())
    nr = np.bincount(radial_dist_int.ravel())
    radial_profile = tbin / nr
    return radial_profile

#fftimages is a list of 2D or 3D images.  The shape of the images must be the same
def average_fft_power_spectrum(fftimages, shape, is_3d=False):
    num_images = len(fftimages)
    radial_dist = fft_radial_distance(shape, is_3d)
    
    # Calculate the maximum bin index, rounding up to the nearest integer
    #max_bin_index = int(np.ceil(np.max(radial_dist)))
    
    # Initialize the sum of radial power spectra
    #sum_radial_power_spectrum = np.zeros(max_bin_index)
    
    image_count = 0
    for fftimage in fftimages:
        radial_power_spectrum, nummodes = radial_bin_average(fftimage, radial_dist)
        if image_count == 0:
            sum_radial_power_spectrum = np.zeros(len(radial_power_spectrum))
        sum_radial_power_spectrum += radial_power_spectrum
        image_count += 1

    
    # Average across all images
    avg_radial_power_spectrum = sum_radial_power_spectrum / num_images
    return avg_radial_power_spectrum, nummodes

def average_fft_cross_power_spectrum(fftimages1, fftimages2, shape, is_3d=False):
    num_images1 = len(fftimages1); num_images2 = len(fftimages2)

    if(num_images1 != num_images2):
        print("number of images in each set must be the same:", num_images1, num_images2, "Aborting")
        return
    
    radial_dist = fft_radial_distance(shape, is_3d)
    
    # Calculate the maximum bin index, rounding up to the nearest integer
    max_bin_index = int(np.ceil(np.max(radial_dist)))
    
    # Initialize the sum of radial power spectra
    sum_radial_power_spectrum = np.zeros(max_bin_index)

    for fftimage1, fftimage2 in zip(fftimages1, fftimages2):
        radial_power_spectrum = radial_bin_average_cross(fftimage1, fftimage2, radial_dist)
        sum_radial_power_spectrum += radial_power_spectrum
    
    # Average across all images
    avg_radial_power_spectrum = sum_radial_power_spectrum / num_images1
    return avg_radial_power_spectrum



#removes the CIC window function for DM (3d only) 
def removeCICWindow(Rf, BoxSize, Na):
    dk = 2. * np.pi / BoxSize
    dx = BoxSize / Na
    kz, ky, kx = np.meshgrid(np.fft.fftfreq(Na), np.fft.fftfreq(Na), np.fft.rfftfreq(Na))
    kz, ky, kx = kz * Na * dk, ky * Na * dk, kx * Na * dk

    denom = (0.5 * dx * kx) * (0.5 * dx * ky) * (0.5 * dx * kz)
    W =  np.sin(0.5 * dx * kx) * np.sin(0.5 * dx * ky) * np.sin(0.5 * dx * kz)
    cond = (np.abs(denom) > 1e-10)
    W[cond] /= denom[cond]
    W[~cond] = 1

    return Rf/W**2
