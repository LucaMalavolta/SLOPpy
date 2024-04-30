### Spectral Lines Of Planets with python
***
### `SLOPpy` version 1.1 - March 2023

**Main contributors**
Daniela Sicilia,
Luca Malavolta,
Lorenzo Pino,
Gaetano Scandariato,
Valerio Nascimbeni,
Giampaolo Piotto,
Isabella Pagano

For **Scientific questions**, please contact Daniela Sicilia
Fro **Technical problems and software maintenance** please contact Luca Malavolta


**News**
Version 1.3 now satisfies the new FITS requirements of `molecfit` version 4.3.1 (as for 30th of April 2024) 
Starting from version 1.1, it is much easier to add new instruments

Basic [documentation is available here](https://sloppy.readthedocs.io/en/latest/).

Most of the information can be found in Sicilia et al. (2022) [A&A link](https://doi.org/10.1051/0004-6361/202244055) [ADS link](https://ui.adsabs.harvard.edu/abs/2022arXiv220813045S/abstract) [arXiv link](https://arxiv.org/abs/2208.13045)

**Note on the use of `molecfitv**

`SLOPpy` supports both the old version `1.5.x` (available as a stand-alone program) and the latest version `>4` through ESO `esorex`.

To use version 1.5.x, you have to use the `molecfit` modules marked with `v1`, and specify the installation path under the `molecfit` section:

```bash
pipeline:
  - ...
  - telluric_molecfit_v1_preparation
  - telluric_molecfit_v1_coadd
  - ...
molecfit:
  installation_path: /Applications/molecfit/bin/
  ...
```

For the newest version, drop the `v1` in the module names and substitute the `installation_path` keyword with the new `esorex_exec' keyword, as in this example:

```bash
pipeline:
  - ...
  - telluric_molecfit_preparation
  - telluric_molecfit_coadd
  - ...
molecfit:
  esorex_exec: esorex
  ...
```

This new keyword specifies the location of the `esorex` executable: if the command is available system-wide (i.e., you can launch esorex by simply writing `esorex` in your terminal), then it is not necessary to specify the full path of the executable.


**Changelog**

- Version 1.3.0
    Default AER version is now "3.8.1.2" are required by Molecfit v. 4.3.1 (April 2024))
    It is now possible to perform Molecfit telluric modeling on individual spectra rather than on coadded spectra
    Temporary FITS files for telluric correction are now compliant with the new requirements of Molecfit 4.3

- Version 1.1.1 to 1.2.2
    Minor bugfixes

- Version 1.1.0
    Easier interface for new instruments
    New module to write transmission spectra in the stellar reference system with (optional) CLV and RM correction: `write_output_transmission_stellarRF`. Note: CLV and RM correction only available through `CLV_RM_models` module

- Version 1.0.2:

    Minor bugfix

- Version 1.0.1:

    Added support to AER v3.8, bundled with Molecfit as in November 2022 \
    To use version 3.6, add this keyword in the configuration file under the molecfit section:\
    NOTE: it must be a string and not a number\

    ```bash
        aer_version: "3.6"
    ```

- Version 1.0: first release after acceptance of the paper.
