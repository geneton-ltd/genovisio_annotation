# Genovisio annotation

[![Python version](https://img.shields.io/badge/python-3.12+-green.svg)](https://www.python.org/downloads/)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

Annotation of genovisio database. Handling data on .json file with data from MongoDB.

## Installation

In python3.12 you can simply use pip to install `annotation`:

```bash
pip install git+https://github.com/geneton-ltd/genovisio_annotation.git
```

Without python3.12, you can install using mamba:

```bash
mamba env create -f conda_example.yaml
```

This gives you 1 entrypoint:

- `annotate` - produce annotation JSON from the mongo database

Also it gives you API operating over the produced JSON annotation file.

## Running

To produce an annotation, running instance of mongo database is required. Mongo URI and database name can be supplied to the entrypoint commands, see `--help`. Default MongoDB URI is `mongodb://localhost:27017/` and the database name 'genovisio'.

To run `annotate` command (if installed using conda, activate it first):

```sh
annotate chr16:34289161-34490212/loss --output annotation.json
```

To use API:

```py
import annotation

annot = annotation.Annotation.load_from_json('annotation.json')
# use annot object ...
```

## Development

There are some changes whether you develop on machines with Python 3.12 and on machines without Python 3.12.
If you cannot install python3.12 and are limited to conda environments, skip to the next guide.

### Development with Python 3.12

Install poetry using pipx:

```sh
pipx install poetry
```

Now in the cloned repository, install the package:

```sh
poetry install
```

Activate the virtual environment where dependencies are installed:

```sh
poetry shell
```

All dependencies are now installed and you can run using:

- entrypoint command `annotate {ARGS}`
- `python {python_script}` like `python annotation/annotation.py {ARGS}`

### Development without python 3.12

Install custom conda environment with python3.12.

Then, install poetry (python packaging management library):

```sh
curl -sSL https://install.python-poetry.org | python3 -
```

Now in the cloned repository, install the package:

```sh
poetry install
```

All dependencies are now installed and you can run using:

- entrypoint command `annotate {ARGS}`
- `python {python_script}` like `python annotation/annotation.py {ARGS}`

### Adding/removing dependencies

When adding or removing dependencies, you need to define them in `pyproject.toml`.

Then, locked versions need to be redefined:

```sh
poetry lock
```

Install again:

```sh
poetry install
```

### Conventional PRs

When committing, you must follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.
This means that there the PR title must be like `feat: XY` or `fix: XY` and so on, and the PR should contain at least one commit named like this.

Further, any push (i.e. after merged PR) to the `main` branch creates in a new PR:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history
