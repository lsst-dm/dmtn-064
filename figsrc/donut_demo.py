import numpy as np
import galsim
import matplotlib.pyplot as plt
from astropy.utils.console import ProgressBar

def donut_demo():
    rng = galsim.BaseDeviate(123)
    u = galsim.UniformDeviate(rng)

    Ellerbroek_alts = [0.0, 2.58, 5.16, 7.73, 12.89, 15.46]  # km
    Ellerbroek_weights = [0.652, 0.172, 0.055, 0.025, 0.074, 0.022]
    Ellerbroek_interp = galsim.LookupTable(Ellerbroek_alts, Ellerbroek_weights,
                                           interpolant='linear')
    spd = []  # Wind speed in m/s
    dirn = [] # Wind direction in radians
    r0_500 = [] # Fried parameter in m at a wavelength of 500 nm.
    for i in range(len(Ellerbroek_alts)):
        spd.append(u()*20)  # Use a random speed between 0 and 20
        dirn.append(u()*360*galsim.degrees)  # And an isotropically distributed wind direction.
        r0_500.append(0.2*Ellerbroek_weights[i]**(-3./5))
        print("Adding layer at altitude {:5.2f} km with velocity ({:5.2f}, {:5.2f}) m/s, "
              "and r0_500 {:5.3f} m."
              .format(Ellerbroek_alts[i], spd[i]*dirn[i].cos(), spd[i]*dirn[i].sin(), r0_500[i]))
    atm = galsim.Atmosphere(r0_500=r0_500, speed=spd, direction=dirn, altitude=Ellerbroek_alts,
                            rng=rng, screen_size=102.4, screen_scale=0.1)

    ab = [0,0]+[0.2]*12
    opt = galsim.OpticalScreen(diam=8.1, aberrations=ab)

    aper = galsim.Aperture(diam=8.1, obscuration=0.2, lam=500.0)
    donut_aper = galsim.Aperture(diam=8.1, obscuration=0.2, lam=500.0,
                                 oversampling=0.5, pad_factor=8.0)

    print("Making pure optics images")

    focal_opt_img = galsim.PhaseScreenList(opt).makePSF(
        lam=500.0, aper=aper
    ).drawImage(nx=128, ny=128, scale=0.01).array

    # Make d(aberrations)/d(z).  This is just made up, not associated with any real data.
    dabdz = [0,0]+[0.1]*12
    dabdz[4] = 27.5

    gsp = galsim.GSParams(maximum_fft_size=8192)

    intra_opt = galsim.OpticalScreen(diam=8.1, aberrations=[a-b for a,b in zip(ab, dabdz)])
    intra_opt_img = galsim.PhaseScreenList(intra_opt).makePSF(
        lam=500.0, aper=donut_aper, gsparams=gsp
    ).drawImage(nx=64, ny=64, scale=0.2).array

    extra_opt = galsim.OpticalScreen(diam=8.1, aberrations=[a+b for a,b in zip(ab, dabdz)])
    extra_opt_img = galsim.PhaseScreenList(extra_opt).makePSF(
        lam=500.0, aper=donut_aper, gsparams=gsp
    ).drawImage(nx=64, ny=64, scale=0.2).array

    # Next make instantaneous opt+atm psfs
    print("Making instantaneous optics+atm images")

    focal_instant_img = galsim.PhaseScreenList([opt, atm]).makePSF(
        lam=500.0, aper=aper, exptime=0.0
    ).drawImage(nx=128, ny=128, scale=0.01).array

    intra_instant_img = galsim.PhaseScreenList([intra_opt, atm]).makePSF(
        lam=500.0, aper=donut_aper, exptime=0.0, gsparams=gsp
    ).drawImage(nx=64, ny=64, scale=0.2).array

    extra_instant_img = galsim.PhaseScreenList([extra_opt, atm]).makePSF(
        lam=500.0, aper=donut_aper, exptime=0.0, gsparams=gsp
    ).drawImage(nx=64, ny=64, scale=0.2).array

    # And finally, the integrated opt+atm psfs
    print("Making integrated optics+atm focal image")
    exptime = 200.0
    with ProgressBar(exptime/0.025) as bar:
        focal_psf = galsim.PhaseScreenList([opt, atm]).makePSF(
            lam=500.0, aper=aper, exptime=exptime, _bar=bar)
        focal_img = focal_psf.drawImage(nx=128, ny=128, scale=0.01).array
    print("Making integrated optics+atm donut images")
    with ProgressBar(2*exptime/0.025) as bar:
        intra_psf = galsim.PhaseScreenList([intra_opt, atm]).makePSF(
            lam=500.0, aper=donut_aper, exptime=exptime, gsparams=gsp, _bar=bar)
        extra_psf = galsim.PhaseScreenList([extra_opt, atm]).makePSF(
            lam=500.0, aper=donut_aper, exptime=exptime, gsparams=gsp, _bar=bar)
        intra_img = intra_psf.drawImage(nx=64, ny=64, scale=0.2).array
        extra_img = extra_psf.drawImage(nx=64, ny=64, scale=0.2).array

    # Make figure

    intra_min, intra_max = np.min(intra_opt_img), np.max(intra_opt_img)
    extra_min, extra_max = np.min(extra_opt_img), np.max(extra_opt_img)

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 8))
    axes[0,0].imshow(focal_opt_img, extent=128*0.5*0.01/0.2*np.r_[-1,1,-1,1])
    axes[1,0].imshow(intra_opt_img, extent=64*0.5*np.r_[-1,1,-1,1], vmin=intra_min, vmax=intra_max)
    axes[2,0].imshow(extra_opt_img, extent=64*0.5*np.r_[-1,1,-1,1], vmin=extra_min, vmax=extra_max)

    axes[0,1].imshow(focal_instant_img, extent=128*0.5*0.01/0.2*np.r_[-1,1,-1,1])
    axes[1,1].imshow(intra_instant_img, extent=64*0.5*np.r_[-1,1,-1,1])
    axes[2,1].imshow(extra_instant_img, extent=64*0.5*np.r_[-1,1,-1,1])

    axes[0,2].imshow(focal_img, extent=128*0.5*0.01/0.2*np.r_[-1,1,-1,1])
    axes[1,2].imshow(intra_img, extent=64*0.5*np.r_[-1,1,-1,1], vmin=intra_min, vmax=intra_max)
    axes[2,2].imshow(extra_img, extent=64*0.5*np.r_[-1,1,-1,1], vmin=extra_min, vmax=extra_max)

    axes[0,0].set_ylabel("Focal")
    axes[1,0].set_ylabel("Intra-focal")
    axes[2,0].set_ylabel("Extra-focal")
    axes[0,0].set_title("Optics only")
    axes[0,1].set_title("Instantaneous")
    axes[0,2].set_title("Integrated")

    fig.tight_layout()
    fig.savefig("donut_demo.png")


if __name__ == '__main__':
    donut_demo()
