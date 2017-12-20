..
  Technote content.

  See https://developer.lsst.io/docs/rst_styleguide.html
  for a guide to reStructuredText writing.

  Do not put the title, authors or other metadata in this document;
  those are automatically added.

  Use the following syntax for sections:

  Sections
  ========

  and

  Subsections
  -----------

  and

  Subsubsections
  ^^^^^^^^^^^^^^

  To add images, add the image file (png, svg or jpeg preferred) to the
  _static/ directory. The reST syntax for adding the image is

  .. figure:: /_static/filename.ext
     :name: fig-label

     Caption text.

   Run: ``make html`` and ``open _build/html/index.html`` to preview your work.
   See the README at https://github.com/lsst-sqre/lsst-technote-bootstrap or
   this repo's README for more info.

   Feel free to delete this instructional comment.

:tocdepth: 1

.. Please do not modify tocdepth; will be fixed when a new Sphinx theme is shipped.

.. sectnum::

.. note::

   **This technote is not yet published.**

   We develop a model for the HSC optical PSF through analysis of out-of-focus "donut" images.

.. _intro:

Introduction
============

We aim to create a point spread function (PSF) model for Hyper Suprime-Cam (HSC) that separates the
contribution of optics from the contributions of sensors and the atmosphere.  "Optics" in this
context includes all of: aberrations built into the HSC design, deviations from that design in the
realized surfaces, misalignments of surfaces, and non-flatness of and discontinuities between the
detectors populating the HSC focal surface.  Modeling the sensor discontinuities is of particular
interest, as these introduce discontinuities into the final delivered PSF.  These discontinuities
restrict current PSF interpolation algorithms to operating on single CCDs at a time.  By modeling
the discontinuities, we hope to enable interpolation of the non-optical PSF components across the
entire focal plane.

We assume that the optics design, deviations in the realized surfaces, and sensor discontinuities
can all be characterized by a single static reference wavefront model.  Misalignments will in
general vary from exposure to exposure, as the telescope and camera experience different
gravitational loads, temperatures, and potentially hysteresis.  The added effects due to
misalignments on the wavefront should be slowly varying in degree of misalignment, and pupil and
focal coordinates.  We believe we can model these effects with a small number of parameters for
every individual exposure.

In this technote, our goal is to demonstrate that we can measure the HSC wavefront from out-of-focus
"donut" images of stars, and that we can use this wavefront to predict the optical component of
in-focus PSFs.  We make this validation using both simulations and observations on donut and
in-focus images.  Ultimately, we intend to also use the out-of-focus exposures to learn the degrees
of freedom that capture the effects of misalignments, which would give us a model that can be
applied to any in-focus image, but that work is left to a future technote.

.. _model:

Forward Model
=============

To model the optical part of the PSF, we use (monochromatic) Fourier optics, described in this
section.

Both a geometric optics model and a Fourier optics model employ three distinct but related
coordinate systems.  Under the geometric model, a point source at a given angular position
:math:`\theta` emits parallel rays of light that, when not vignetted, fill the pupil (essentially
the primary mirror) at points :math:`u`, and then get focused by the optics into (ideally) a single
point :math:`x` on the focal plane.  (My favorite description of a telescope is exactly this: a
device that transforms angles into positions.  Also, note that :math:`\theta`, :math:`u`, and
:math:`x` are all 2-dimensional.)

In practice, aberrations are always present which means that the light from a single incoming angle
gets distributed over a non-zero area of the focal plane.  The point spread function :math:`I(x;
\theta)` describes the surface brightness distribution over focal plane coordinates :math:`x` for
light originating from a single point on the sky :math:`\theta`.

Fourier optics extends the geometric model by incorporating the effects of the wave nature of light.
Under Fourier optics we picture not parallel rays of incoming light, but incoming plane waves of
light with the planes perpendicular to the rays.  The job of the telescope in this picture is to
transform an incoming plane wave into a spherical wave that converges to a single point on the focal
plane.  Fourier optics essentially decomposes this spherical wave into plane waves with the correct
relative phases, which are then allowed to interfere along the focal plane, producing the
diffractive PSF.  The role of aberrations under this picture is to alter those phases from what they
would be for a perfect converging spherical wave.  In particular, we define the wavefront
:math:`W(u; \theta)` as the lag or lead (in dimensions of length) of the delivered wave from a
perfect converging spherical wave (to convert this length into a phase, multiply by :math:`\frac{2
\pi}{\lambda}`, where :math:`\lambda` is the wavelength of the light).

The final ingredient for our Fourier optics model is the pupil obscuration function :math:`P(u;
\theta)`, a binary-valued function that indicates, for a given incoming angle :math:`\theta`, which
points :math:`u` on the pupil are unobscured all the way to the focal plane.

The equation for the Fourier optics PSF is:

.. math::
    I(x; \theta) \propto \left| \mathfrak{F}\left[P(u; \theta) \exp\left(\frac{-2\pi i}{\lambda} W(u; \theta)\right)\right] \right|^2

where the Fourier transform operator :math:`\mathfrak{F}` is an integral over :math:`u` yielding
a function of :math:`x`.

The goal of the forward model is to infer the wavefront from a given PSF image and knowledge of the
obscuration function.  For this purpose, it is convenient to decompose the wavefront into a Zernike
polynomial series:

.. math::
    W(u; \theta) = \sum_{j=4}^{j_\mathrm{max}} a_j(\theta) Z_j(u).

where :math:`a_j(\theta)` are a set of :math:`\theta`-dependent coefficients, and :math:`Z_j(u)` is
the :math:`j`-th Zernike polynomial.  We start at :math:`Z_4` because this is the first Zernike term
that actually affects the profile of the PSF (:math:`Z_1` is a phase term, the effect of which
completely disappears after the :math:`|\cdot|^2` operation, and the :math:`Z_2, Z_3` terms are
perfectly degenerate with coordinate origin shifts).  For the moment, we will leave the functional
form for the field angle :math:`\theta` dependence of the Zernike coefficients unspecified, though
this will be important when we develop the wavefront model separating static effects and dynamic
effects eluded to above.

Donut PSFs
----------

The final delivered PSF also includes contributions from sensors (primarily in the form of charge
diffusion), and the atmosphere.  The atmospheric component nominally follows from the same Fourier
optics formalism as the optics component, with the main differences being that the atmospheric
wavefront aberrations are typically much larger than optical aberrations and they evolve on scales
of milliseconds instead of being essentially constant during an exposure.  An important consequence
of this, combined with the nonlinearity of the equation for :math:`I(x)`, is that the time averaged
wavefront is not simply related to the time averaged PSF.

In order to develop our wavefront-based optics PSF model, we need a way to measure the wavefront.
Measuring the wavefront ab initio from in-focus data is nearly impossible due to the confounding
atmospheric PSF.  However, by defocusing the camera, we can amplify the effects of the optics
wavefront while keeping the effects of the atmospheric PSF roughly constant.  The caveat to this
approach is that optical aberrations change when defocusing the camera.  Fortunately, this change is
linear with defocus (verified in Zemax), so an effective strategy is to measure both the intra-focal
(sensors placed in front of the focal plane) and extra-focal wavefront (sensors placed behind the
focal plane) and average them together to get the in-focus wavefront prediction.

.. figure:: /_static/donut_demo.png
    :name: donut_demo
    :target: ../../_static/donut_demo.png

    GalSim simulation illustrating how donut images enable optical wavefront reconstruction.  Note
    the difference in scale between the rows.  The first column shows simulated images including
    only optics effects.  The non-circularity of the images (both in and out of focus) is related to
    non-circularity in the  simulated wavefront aberrations.  The second column shows the same
    simulation, but including atmospheric effects (but with an exposure time of 0 seconds).  The
    change to the in-focus PSF is dramatic compared to the change in the donuts.  The right column
    integrates the atmospheric effects over 200 seconds.  While it would be very difficult to infer
    anything about the optics from the in-focus integrated image, the integrated donuts still retain
    most of the optics information.  The source code can be found in this technote's repository.

Donut Fits
----------

All of the wavefront inferences presented in this technote are derived from fitting a forward image
model :math:`M_p(\Theta)` to donut images :math:`D_p` (and, when available, pixel variances
:math:`V_p`) in order to minimize

.. math::

    \chi^2 = \sum_p \frac{\left(M_p(\Theta) - D_p\right)^2}{V_p}

where the subscript :math:`p` indexes individual pixels.  The parameters :math:`\Theta` of the
forward image model are the centroid :math:`x` and :math:`y`, the flux, the wavefront Zernike
coefficients :math:`a_j` up to some fixed order :math:`j_\mathrm{max}` (starting at :math:`j=4`),
and finally a single parameter to describe the additional blurring of the image due to the
atmospheric PSF (although using a single parameter to characterizes the atmospheric PSF contribution
is unlikely to be a good model for in-focus data, for fitting donuts it's a good approximation).  In
our case, we choose to model this additional PSF component by convolving the optical PSF model with
a Kolmogorov profile with Fried parameter :math:`r_0`.

Note that to perform the fit, a wavelength :math:`\lambda` and pupil obscuration function
:math:`P(u)` also need to be specified.  In the case of simulated data, the wavelength is known
exactly, and the pupil obscuration can be obtained from the same software used to create the
simulated data.  For real data, we use the observed filter effective wavelength and infer pupil
obscurations from pinhole images described below.

We perform the fit using the Levenberg-Marquardt algorithm implemented in the python package lmfit
(which wraps the implementation in scipy).

.. _zemax:

Model validation in Zemax
=========================

To validate the wavefront-based PSF modeling, we used the commercial raytracing package Zemax and
the S402C description of HSC and the Subaru telescope to simulate sets of intra-focal, extra-focal,
and in-focus images.  We used the Zemax HuygensPSF tool to simulate in-focus optical PSF images at a
resolution of 0.25 :math:`\mu m`, which roughly corresponds to 0.0028 arcseconds at the mean HSC
pixel-scale of 0.168 arcseconds per 15 :math:`\mu m` pixel.  To simulate donut images, we first
displaced (in Zemax) the camera from the primary mirror by :math:`\pm 0.9` mm along the optic axis.
We again used the HuygensPSF tool to generate an image, this time at the native pixel scale
resolution of 15 :math:`\mu m`.  To make these images more realistic, we used GalSim to convolve
them by a Kolmogorov profile and add a small amount of uncorrelated stationary Gaussian noise (this
extra convolution and noise addition also seems to help our fitter converge).  We also used Zemax to
determine the pupil obscuration function exactly.

After finding the wavefront coefficients by fitting the intra-focal and extra-focal donuts, we
averaged them together to produce an inferred wavefront for the in-focus optical PSF.  We then
compared the inferred in-focus wavefront and PSF to the simulation truth obtained from Zemax.

We validated the forward model approach using 2 configurations of the simulated camera/telescope:
one in which the optics are perfectly aligned, and one in which the camera is slightly shifted and
tilted with respect to the primary mirror optic axis.  We investigated 4 points in the field of view
spanning the complete range in incoming field angle of 0 to 0.75 degrees.

The following multi-panel figures show the results for one of these images sets.  In the particular
set shown, the camera has been misaligned and the field angle is near maximal at 0.75 degrees.

.. figure:: /_static/pf6_intra_vs_jmax.png
    :name: donut_vs_jmax
    :target: ../../_static/pf6_intra_vs_jmax.png

    Intra-focal donuts fit validation.  The top row shows forward model fits to the Zemax-simulated
    images using progressively more Zernike polynomials in the wavefront description.  The middle
    row shows the particular Zemax donut being fit (each column is identical in this row).  The
    bottom row shows the residuals (in the same arbitrary units as the data and model).  The
    analogous figure for extra-focal donuts is qualitatively similar, but omitted here for
    concision.

.. figure:: /_static/pf6_WF_vs_jmax.png
    :name: WF_vs_jmax
    :target: ../../_static/pf6_WF_vs_jmax.png

    Wavefront inference validation.  The top row shows the wavefront inferred from the donut image
    pairs using progressively more Zernike polynomials in the wavefront description (in units of
    waves).  The middle row shows the true wavefront from Zemax.  The bottom row shows the
    residuals.

.. figure:: /_static/pf6_PSF_vs_jmax.png
    :name: PSF_vs_jmax
    :target: ../../_static/pf6_PSF_vs_jmax.png

    Optical PSF inference validation.  The top row shows the optical PSF inferred from the donut
    image pairs using progressively more Zernike polynomials in the wavefront description.  The
    middle row shows the true optical PSF from Zemax.  The bottom row shows the residuals.

One important point here is that at no point did we do anything special to account for distortion in
the optics.  That is, the simulated images created by Zemax use focal plane coordinates with
physical units of mm, and in general have a nonlinear relationship with sky-coordinates or even a
tangent plane projection thereof.  Investigating the potential impact of distortion is on our list
of open questions.

.. _data:

Model application to real data
==============================

To test the optical PSF model on real data, we rely on a set of HSC engineering images: visits 69008
through 69072.  These images were taken in sets of three (which we will refer to as a "triplet"),
alternating through an in-focus exposure, an extra-focal exposure (with camera displaced by +0.9
mm), and an intra-focal exposure (with camera displaced by -0.9 mm).  Each exposure in a triplet was
taken at the same location on the sky, which allows us to directly compare intra- and extra-focal
donut images with in-focus PSFs for a fixed set of focal-plane locations, without needing to worry
about interpolating across the focal plane.  The telescope elevation and rotation angle were varied
from one triplet to the next.  All triplet images in this range were taken through the I2 filter.

Estimating the pupil obscuration function
-----------------------------------------

Recall that we assume the pupil obscuration function :math:`P(u; \theta)` is known a priori during a
given wavefront fit.  This function varies considerably across the HSC field-of-view due to the
significant vignetting of the camera, and also has significant contributions from the telescope
spiders and camera shadow.

Pinhole images
^^^^^^^^^^^^^^

To estimate the pupil obscuration and its variation, we use a set of HSC images taken through a
pinhole filter and illuminated by a flat screen.  Each pinhole forms an image of the pupil on the
focal plane.  We choose to describe this image using a combination of three circles and four
rectangles.  The first circle is used to indicate the boundary of the primary mirror (shown in blue
below).  The second circle indicates the shadow formed by HSC itself (shown in green below).  The
third circle shows where rays are clipped by the first HSC lens (shown in red below).  Using ds9, we
match the edges of these circles to the edges formed by the pinhole images and record the positions
of each circle center.

.. figure:: /_static/pinholes.png
    :name: Pinholes
    :target: ../../_static/pinholes.png

    Screenshot showing images through pinhole filter and circles used to characterize the pupil.

.. figure:: /_static/pinholes_zoom.png
    :name: Pinholes Zoom
    :target: ../../_static/pinholes_zoom.png

    Zoom-in on one pinhole image and circles used to model pupil.  Note the non-uniform illumination
    from the (suboptimal) flat field screen.

We next investigate how the circles centers relate to one another as a function of focal plane
position.  The plots below indicate that this variation is quite close to linear in the field
radius.

.. figure:: /_static/camera_shadow_displacement.png
    :name: camera_shadow_displacement
    :target: ../../_static/camera_shadow_displacement.png

    Camera shadow displacement with respect to primary mirror across the focal plane.

.. figure:: /_static/camera_shadow_fit.png
    :name: camera_shadow_fit
    :target: ../../_static/camera_shadow_fit.png

    Camera shadow displacement radial fit.

.. figure:: /_static/l1_displacement.png
    :name: l1_displacement
    :target: ../../_static/l1_displacement.png

    HSC lens clipping displacement with respect to primary mirror across the focal plane.

.. figure:: /_static/l1_fit.png
    :name: l1_fit
    :target: ../../_static/l1_fit.png

    HSC lens clipping displacement radial fit.

Pinhole != pupil
^^^^^^^^^^^^^^^^

While the pinhole images are a valuable source of data from the real HSC instrument, we note that
the images formed this way are not strictly the same as the pupil.  Each individual pinhole image is
formed by light encountering the optics (specifically, the primary mirror) from a variety of angles,
and then passing through the filter plane within a narrow transverse aperture (i.e., within the
pinhole).  By contrast, true pupil rays for a given field angle are initially parallel, and are
unconstrained in the filter plane (although in practice, because the filter plane is near the focal
plane, only a relatively small region of the filter will intersect the incoming beam for a given
field angle).  In the following figure, we use the raytracing software batoid [#]_ to confirm that
this distinction does indeed produce different images, though we have not, as yet, propagated these
differences to see their potential impact on wavefront or PSF inference.

.. figure:: /_static/Pupil_vs_pinhole.png
    :name: pu_vs_ph
    :target: ../../_static/Pupil_vs_pinhole.png

    Comparison of pinhole and pupil images for HSC.  Simulations produced by batoid.  For the pupil
    images, parallel rays are traced through optics and obscurations until the surviving rays impact
    the detector plane.  For the pinhole images ray fans over a small solid angle are traced in
    opposite directions starting from locations in the filter plane (the locations of pinholes).
    Rays impacting the focal plane that are exactly opposite surviving rays exiting the telescope
    entrance pupil are retained to produce the pinhole images.  Four different field angles are
    shown.

It should be possible to empirically measure the true pupil by using images of stars taken very far
from focus (much farther from focus than the donut images analyzed below).  As the defocus is
increased, interference effects become insignificant compared to geometric effects, allowing the
pupil to be cleanly observed.

.. [#] https://github.com/jmeyers314/batoid

Individual donut fitting
------------------------

For measuring wavefronts, it's important to select reasonably bright isolated objects that originate
from point sources and not extended sources.  To this end, we use three criteria to select donuts
for further analysis:

1. High signal-to-noise ratio.
2. The donut "hole" is significant.  This feature would be washed out for extended sources.
3. The object is not too large, which may indicate blending of neighboring sources.

For the above, we use base_CircularApertureFlux_25_0_flux / base_CircularApertureFlux_25_0_fluxSigma
as the signal-to-noise ratio statistic, i.e., the relative signal-to-noise ratio of a circular
aperture flux measurement with a radius of 25 pixels.  For the donut hole statistic, we use the
ratio base_CircularApertureFlux_3_0_flux / base_CircularApertureFlux_25_0_flux, and for the
size statistic, we use base_CircularApertureFlux_35_0_flux / base_CircularApertureFlux_25_0_flux.

.. figure:: /_static/donutSelection-0069030-046.png
    :name: selection
    :target: ../../_static/donutSelection-0069030-046.png

    Example of donut selection.  Blue donuts are selected for fitting; red donuts are rejected.
    These particular donuts are from HSC CCD 46, which is located towards the edge of the HSC field
    of view.  As such, there is significant vignetting visible.  Notice the rejected donut in the
    second row, which has a visible blend.  The other three rejected donuts clearly originate from
    extended sources.  Also notice the unflagged artifact in the selected donut in the middle column
    of the third row.  The fit to this donut may be unreliable.

We fit each selected donut independently using the model specified above.  To improve the
convergence of the model, we fit iteratively, increasing the value of :math:`j_\mathrm{max}` in each
iteration (from 4 to 11 to 15 to 21), and using the results of the previous iteration to initialize
the parameter values for each subsequent iteration.  Sample results for :math:`j_\mathrm{max} = 21`
are shown in the figure below.

.. figure:: /_static/donutGoodnessOfFit-0069030-046.png
    :name: fit
    :target: ../../_static/donutGoodnessOfFit-0069030-046.png

    Example donut fits.  The left column shows the data, the middle column shows the best fitting
    model, and the last column shows the residual.  Note that these donuts are all from the same
    exposure and the same CCD, which is why they all look similar.

While the fits are plausible, there is clearly structure in the data not being captured by the
model.  It may be possible to improve the fits by increasing :math:`j_\mathrm{max}`, at the expense
of increased computational time to perform the fits, and potentially increased degeneracy between
fit parameters.

Wavefront variation across the focal plane
------------------------------------------

Donut fits
^^^^^^^^^^

The following figures show the variation of donuts and fits over the focal plane.

.. figure:: /_static/donutStampData-0069030.png
    :name: dataFOV
    :target: ../../_static/donutStampData-0069030.png

    Donut data over an entire HSC field of view.  The patterns have a smoothly
    varying structure across the field.

.. figure:: /_static/donutStampModel-0069030.png
    :name: modelFOV
    :target: ../../_static/donutStampModel-0069030.png

    Best fitting models with :math:`j_\mathrm{max}=21`.

.. figure:: /_static/donutStampResid-0069030.png
    :name: residFOV
    :target: ../../_static/donutStampResid-0069030.png

    Residuals.

The residuals appear to vary smoothly over the focal plane.  Features are coherent over scales of
many CCDs (i.e., over 10s of arcminutes).  Some features can even be picked out over most of the
focal plane.

The longer the coherence scale of a feature, the closer to the pupil plane it must originate.  That
is, if features in the wavefront or wavefront residuals varied rapidly across the focal plane, they
could not originate on the primary mirror or in the lower atmosphere, as these effect incoming beams
roughly equally.  Conversely, a given wavefront effect located near the focal plane or high up in
the atmosphere will only affect a narrow range of focal plane positions (or equivalently, a narrow
range of incoming angles).

Following this logic, it appears that a significant proportion of the wavefront residuals may
originate near the Subaru primary mirror.

Wavefront
^^^^^^^^^

As in the Zemax tests, we predict the wavefront for in-focus images as the average of the inferred
intra-focal and extra-focal wavefronts.  (This assumes that the intra-focal and extra-focal camera
displacements are precisely equal).  The in-focus wavefront field-of-view variation, along with the
pupil obscuration function are shown in the following figure.

.. figure:: /_static/donutPairStampWavefront-0069028.png
    :name: wavefrontFOV
    :target: ../../_static/donutPairStampWavefront-0069028.png

    Model wavefront.

PSF
^^^

The purely optical PSF implied by the wavefront plotted above is shown below.

.. figure:: /_static/donutPairStampPsf-0069028.png
    :name: psfFOV
    :target: ../../_static/donutPairStampPsf-0069028.png

    Model optical PSF.  The size of each postage stamp is 0.64 arcseconds on a side.

Wavefront coefficients
^^^^^^^^^^^^^^^^^^^^^^

Another way to look at these results is to plot the pupil plane wavefront coefficients as functions
of focal plane location:

.. figure:: /_static/donutZernikePyramid-0069028.png
    :name: donutZernikePyramid
    :target: ../../_static/donutZernikePyramid-0069028.png

    Wavefront coefficients (i.e., :math:`a_j(\theta)` ) across the focal plane.

The coefficients of Zernike polynomial terms that vary like :math:`\cos(n \theta_\mathrm{pupil})` in
the pupil show variation roughly proportional to :math:`\cos(n \theta_\mathrm{FOV})` in the field of
view.  This is a simple consequence of the nearly circular symmetry of the HSC optical system.  The
amplitudes of the coefficients are also diminishing as the index increases, which presumably means
that despite the presence of significant residuals in the donut fits, we are capturing the most
important wavefront features.

Model prediction comparison to measured PSFs
--------------------------------------------

The triplets of extrafocal, intrafocal, and in-focus images enable a particular check on the
accuracy of the wavefront-based optical PSF reconstruction.  Since the in-focus images were taken
nearly contemporaneously with the donut images on the same field, the state of the optics (other
than focus) should be nearly identical in all three images of a triplet.  The in-focus PSF contains
a significant contribution from the atmosphere, of course, making a direct comparison of data to
model difficult. However, if we approximate the atmospheric contribution as a convolution by a
constant isotropic surface brightness profile, then we can convolve the model optics PSF by this
fiducial atmospheric PSF to produce a profile more directly comparable to the in-focus data.

The plots below show predicted and observed moments of the PSF across the HSC field of view.  The
moments of each prediction and observation are summarized in "whiskers", where the orientation of
each whisker indicates the orientation of the PSF defined by:

.. math::
    \beta = \arctan(e_2, e_1)/2

and the length of each whisker indicates the ellipticity defined by:

.. math::
    e = \sqrt{e_1^2+e_2^2}

and :math:`e1` and :math:`e2` are related to the second moments of the PSF by

.. math::
    e_1 = \frac{M_{xx} - M_{yy}}{M_{xx}+M_{yy}},

    e_2 = \frac{2 M_{xy}}{M_{xx}+M_{yy}}.

The whisker comparison plot for the triplet corresponding to the earlier figures in this section is
immediately below, followed by whisker plots for the other available triplets.

.. figure:: /_static/donutTripletWhisker-0069026.png
    :name: whisker26
    :target: ../../_static/donutTripletWhisker-0069026.png

    Whisker plot comparison for model derived from visits intra/extra focal visits 69028 and 69030
    (left) and infocus data taken in visit 69026.

.. figure:: /_static/donutTripletWhisker-0069008.png
    :name: whisker08
    :target: ../../_static/donutTripletWhisker-0069008.png

    Whisker plot for triplet (69008/69010/69012).

.. figure:: /_static/donutTripletWhisker-0069014.png
    :name: whisker14
    :target: ../../_static/donutTripletWhisker-0069014.png

    Whisker plot for triplet (69014/69016/69018).

.. figure:: /_static/donutTripletWhisker-0069032.png
    :name: whisker32
    :target: ../../_static/donutTripletWhisker-0069032.png

    Whisker plot for triplet (69032/69034/69036).

.. figure:: /_static/donutTripletWhisker-0069038.png
    :name: whisker38
    :target: ../../_static/donutTripletWhisker-0069038.png

    Whisker plot for triplet (69038/69040/69042).

.. figure:: /_static/donutTripletWhisker-0069050.png
    :name: whisker50
    :target: ../../_static/donutTripletWhisker-0069050.png

    Whisker plot for triplet (69050/69052/69054).

.. figure:: /_static/donutTripletWhisker-0069056.png
    :name: whisker56
    :target: ../../_static/donutTripletWhisker-0069056.png

    Whisker plot for triplet (69056/69058/69060).

.. figure:: /_static/donutTripletWhisker-0069062.png
    :name: whisker62
    :target: ../../_static/donutTripletWhisker-0069062.png

    Whisker plot for triplet (69062/69064/69066).

The whisker plot comparisons show clear correlations between the predicted PSF moments and the
observed PSF moments, indicating that we are on the right path towards optical PSF modeling.  While
the differences between predicted and observed whiskers are still under investigation, we would like
to point out that the current model has no freedom for atmospheric variation across the field of
view, or other contributions to the delivered PSF such as tracking errors / wind-shake, or effects
originating in the sensors.

.. _questions:

Open questions
==============

A number of questions regarding donut-inferred wavefront analysis remain, which we list here:

- What is the impact of inferring the pupil obscuration function from the pinhole images?

- What is the impact of modeling donuts and PSFs monochromatically?

- How does distortion affect donuts or inferred in-focus PSFs?

- Would truncating the Zernike series at a larger order improve the fits?  Would this improve the
  whisker plots?

- Can we verify in a ray-tracing package that misalignments really do only introduce changes in
  Zernike coefficients that vary slowly with field angle?

- How does one infer the degrees of freedom in the wavefront due to misalignments.  How does one
  then use this information to infer the optical PSF of a generic in-focus image?

.. Add content here.

.. .. rubric:: References

.. Make in-text citations with: :cite:`bibkey`.

.. .. bibliography:: local.bib lsstbib/books.bib lsstbib/lsst.bib lsstbib/lsst-dm.bib lsstbib/refs.bib lsstbib/refs_ads.bib
..    :encoding: latex+latin
..    :style: lsst_aa
