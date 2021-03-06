# AUG
Access and process data from the AUG tokamak.

The software developed here is ready to be used in the `toki##` machines, providing you have an MPCDF account and `ipfnpytools` installed.
___
## Setting up

### On a toki## machine

Clone the `ipfnpytools` git repository and add it to your python path.

### On your own PC

You can use `toki.yml` to create a conda environment with the basic packages needed to run these notebooks. You will also need to install ipfnpytools.


### Plotting styles

In the 'Styles' directory you can find several plotting styles.

* "darklab": A style suitable to use with jupyter-lab's dark theme. Use transparent=True when using plt.savefig().
* "helvet2": Uses latex fonts and helvetica.
* "helvet2dark": Same as previous, but for use with jupyter-lab's dark theme.

Most notebooks are set to use the dark theme from jupyter-lab and they will use the "darklab" theme as the matplotlib style.

