(installation)=

# Installing SLOPpy

## Setting up an environment

Before proceeding with the installation, I suggest to create an environment dedicated to `SLOPpy` using python\<=3.9 .  At the moment of writing I hve received a few complaints (unrelated to SLOPpy) about Python 3.10, so you may use it at your own risk.
With conda/anaconda:

```{code} bash
conda create --name sloppy python=3.9
```

To list the available environments do:

```{code} bash
conda env list
```

The active environment will be marked with a \*

To activate the `sloppy` environment:

```{code} bash
WINDOWS: activate sloppy
LINUX, macOS: conda activate sloppy
```

## Install using pip

You can then install `SLOPpy` by using `pip` inside the code repository:

```{code} bash
 pip install sloppy-package
```

Note that the name is  `sloppy-package` and not  `sloppy`, as the former was already taken by another package in PyPI. The name for package importation will still be `SLOPpy`:

```{code} bash
 python -c "import SLOPpy"
```

## Install from the repository

Download the latest version from the GitHub repository:

```{code} bash
 git clone https://github.com/LucaMalavolta/SLOPpy.git
```

You can then install `SLOPpy` by using `pip` inside the code repository:

```{code} bash 
 cd SLOPpy
 pip install .
```

Alternatively, you can install `SLOPpy` using the `setup.py` file:

```{code} bash
 cd SLOPpy
 python setup.py install
```

## Requirements

```{admonition} Give people credit for their work

If you are using any of those packages listed above, *please be sure to cite the proper references*, as stated in the relative web page. 
```

These packages are installed automatically when using pip.

- `numpy`, `scipy`, `matplotlib`: pretty standard
- `numba`: open source JIT compiler, actually required as undeclared dependency by some packages ([numba home page])
- `argparse`: Parser for command-line options, arguments and sub-commands, required to pass terminal keywords ([argpares home page])
- `oyaml`: a full-featured YAML framework for the Python programming language.  YAML is the language used for the configuration file ([oyaml home page], [yaml home page])
- `pygtc`: Make a publication-ready giant-triangle-confusogram (GTC) ([pygtc home page])
- `pyDE`: global optimization package ([PyDE home page])
- `emcee`: ensemble sampling toolkit for affine-invariant MCMC ([emcee home page]). 

`emcee` is already included in the requirements, `pyDE` needs to be installed separately as the GitHub version supports multiprocessing: 

```{code} bash
 pip install git+https://github.com/hpparvi/PyDE.git
```

THe alternative ensemble slice sampler `zeus` is supported as well ([zeus home page]).

````{tip}
 `pyDE` and all the additional requirements can be installed by locating the `extra_requirements.txt` file in the [SLOPpy repository] or by downloading it directly [from here](https://raw.githubusercontent.com/LucaMalavolta/SLOPpy/main/extra_requirements.txt) and then running from a terminal:

 ```{code} bash
 pip install -r extra_requirements.txt
 ```

````

[extra_requirements.txt]: https://github.com/LucaMalavolta/SLOPpy/blob/main/extra_requirements.txt

[SLOPpy repository]: https://github.com/LucaMalavolta/SLOPpy

[numba home page]: https://numba.pydata.org/
[pygtc home page]: https://pygtc.readthedocs.io/
[argpares home page]: https://docs.python.org/3/library/argparse.html
[oyaml home page]: https://pypi.org/project/oyaml/
[yaml home page]: https://yaml.org/
[emcee home page]: https://emcee.readthedocs.io/
[h5py home page]: http://docs.h5py.org/
[PyAstronomy home page]: https://pyastronomy.readthedocs.io/
