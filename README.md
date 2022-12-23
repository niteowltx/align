This code is a C implementation of the avisynth/avxsynth depan plugin. It
takes as input a list of jpg files that are the frame-by-frame capture
of a 16mm, 8mm or super8 film.  These frames are roughly aligned, but
may not be 'exactly' aligned due to mechanical variation during capture.
The depan algorithm is used to look at a particular spot on the image that
is usually the edge of the image that contains the sprocket holes.  A config
file is supplied that tells the program where to find a reference frame and
where on that image the sprocket holes are 'supposed' to be.  All input
frames are then compared to the reference frame and an estimate of how much
to shift the frame is made using a 2D FFT algorithm.  Each input frame is then
rotated to match the reference frame.  Once rotated, each output frame is cropped according
to parameters in the config file and an aligned frame is written.  Aligned frames
can then be combined into a video using the make_mp4 script.

If the config file specifies debug=true, the 2D correlation array is written
as a jpg with the gradient shown as various shades of red(worst) to green(best).

The align_dir script works in a directory containing all of the input jpg's
and a config file to produce a final mp4 output file.
