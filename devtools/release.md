# How to release


## Setup

Install required packages:
```
$ conda env create -f environment.yml
```

Activate env:
```
$ conda activate seq-to-first-iso
```


## Tests

Before any release, double-check all tests had run successfully:
```
$ make test-coverage
```


## Update version number

We use the tool [bumpversion](https://github.com/peritus/bumpversion) to update and synchronize the version number
across different files:
```
$ pipenv run bumpversion --verbose --config-file devtools/bumpversion.cfg patch
$ git push origin
$ git push origin --tags
```


## Publish in PyPI

Create the file `$HOME/.pypirc`:
```
[distutils]
index-servers = pypi

[pypi]
repository=https://upload.pypi.org/legacy/
username=<your-login-on-PyPI>
```

Build the package:
```
$ make compile
```

Then upload the package to PyPI:
```
$ make upload-to-pypi
```

Enter your password upon request.

The new package should be available here: https://pypi.org/project/seq-to-first-iso/


### Add new release on GitHub

On [GitHub release page](https://github.com/pierrepo/autoclasswrapper/releases) :

- Click the *Draft a release* button.
- Select the latest version as *tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy and paste the content of the `CHANGELOG.md` in the *Describe this release* field.
- Hit the *Publish Release* button.


### Zenodo integration

For Zenodo integration, see [Making Your Code Citable](https://guides.github.com/activities/citable-code/).

After the creation  of the new release in GitHub, check the archive has been creating on [Zenodo](https://zenodo.org/deposit).
