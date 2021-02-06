import numpy as np




def FraunhoferDiffraction(nearField):

    temp = np.fft.fft2(nearField)

    temp /= GenShift(temp.shape[0])

    return np.fft.fftshift(temp)



def GenShift(npts):

    if npts % 2 == 1 :
        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(npts)*(npts-1)/npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T

    else:
        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(npts)*(npts)/npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T
