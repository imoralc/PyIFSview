<center> <h1>PyIFSview</h1> </center>

PyIFSview is a tool developed to interactively visualize integral field spectroscopy (IFS) data,  such as CALIFA, MaNGA, SAMI or MUSE data.

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Method</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#results">Results</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#cite">Cite</a></li>
    <!-- <li><a href="#acknowledgments">Acknowledgments</a></li> -->
  </ol>
</details>

## Getting Started:

We try to make the usage of PyIFSview as simple as possible. For that purpose, we have create a PyPI package to install PyIFSview.

#### Software Requirements:

- `Numpy`: Python library used for working with arrays
- `Matplotlib`: Python 2D plotting library
- `Astropy`: astronomy library

## Installation:

Using pip you can either install the last relase by

```sh
  pip install PyIFSview
```

or you can install the latest version of the code as 

```sh
  git clone git@github.com:imoralc/PyIFSview.git
  cd PyIFSview
  pip install -e .
```

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->

## Examples:

In here we show how to use PyIFSview for a IFS datacube (three dimensional)
First of all we have to import our library previously install and some dependecies

```python
    from PyIFSview import PyIFSview
```

Then we read the [CALIFA datacube]([https://github.com/PabloMSanAla/fabada/blob/master/examples/bubble.png](https://github.com/imoralc/PyIFSview/blob/main/NGC2253.fits.gz)) borrowed from the [CALIFA Legacy Survey](https://califa.caha.es/).

```python
    PyIFSview(cube_fits_file = 'NGC2253.fits.gz')
```

And its done :wink:
It is as easy as one line of code.

The results obtained running this example would be:

![PyIFSview example](Example.png "Example image using CALIFA data")

Screenshot of the code using data from the CALIFA Survey. The spectrum and the slide can be changed just by clicking on any spaxel in the slide plot and moving the horizontal line in the spectrum plot to select any wavelength, respectively. Additionally, the colormap can be also changed interactively just by clicking on any of the list. It can work with any colormap included in Matplotlib or designed by the user.



![PyIFSview example video](Example_PyIFSview.mp4 "Example video using CALIFA data")


<p align="right">(<a href="#top">back to top</a>)</p>
