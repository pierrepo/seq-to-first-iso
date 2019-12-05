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

Export fully reproductible conda environment:
```
$ conda env export --no-builds | grep -v "^prefix:" > environment.lock.yml
```


## Tests

Before any release, double-check all tests had run successfully:
```
$ make test-coverage
```


## Linting

Also check code complies with PEP 8 and PEP 257:
```
$ make lint
```


## Update version number

We use the tool [bumpversion](https://github.com/peritus/bumpversion) to update and synchronize the version number
across different files:
```
$ bumpversion --verbose --config-file devtools/bumpversion.cfg patch
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


## Add new release on GitHub

On [GitHub release page](https://github.com/pierrepo/seq-to-first-iso/releases) :

- Click the *Draft a new release* button.
- Enter the latest version as *Tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy and paste the content of the `CHANGELOG.md` in the *Describe this release* field.
- Hit the *Publish Release* button.


## Publish in bioconda

Following this doc: https://bioconda.github.io/contributor/updating.html
Calculate sha256 sum for the latest release archive:
```
$ wget -O- https://github.com/pierrepo/seq-to-first-iso/archive/v1.0.0.tar.gz | shasum -a 256
```

Update the sha256 sum in `meta.yml`, commit and push.

If not already done, create a [personnal access token](https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line).

Export this token in `GITHUB_TOKEN`:
```
$ export GITHUB_TOKEN=token....
```

Send `bioconda-utils` will automatically submit a pull request:
```
$ bioconda-utils update recipes/ meta.yml \
  --packages seq-to-first-iso \
  --create-pr
```

At this point the command file:
```
bioconda-utils: error: invalid choice: 'update'
```

However, it seems bioconda package is automatically updated upon new release on GitHub...


## Zenodo integration

For Zenodo integration, see [Making Your Code Citable](https://guides.github.com/activities/citable-code/).

After the creation  of the new release in GitHub, check the archive has been creating on [Zenodo](https://zenodo.org/deposit).
