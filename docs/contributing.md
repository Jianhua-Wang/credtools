# Contributing

Contributions are welcome. The best contributions are small, tested, and easy
to review.

## Set Up a Development Environment

```bash
git clone https://github.com/Jianhua-Wang/credtools.git
cd credtools
uv sync
```

If you do not use `uv`:

```bash
python -m pip install -e ".[dev]"
```

## Run Tests

```bash
uv run pytest
```

Run a smaller test file while developing:

```bash
uv run pytest tests/test_qc.py
```

## Format and Lint

Use the project configuration in `pyproject.toml` and `setup.cfg`.

```bash
uv run black credtools tests
uv run isort credtools tests
uv run flake8 credtools tests
```

## Build the Docs

```bash
uv run mkdocs build
```

Docs should be written in clear English. Prefer tutorial language over academic
language. Introduce a concept before asking the reader to use it.

## Pull Request Checklist

- The change is scoped to one clear problem.
- Tests cover the behavior change.
- CLI behavior is documented when a command changes.
- New public functions have docstrings.
- The docs build with `mkdocs build`.

## Reporting Issues

Open an issue on GitHub and include:

- CREDTOOLS version,
- Python version,
- operating system,
- command you ran,
- short input description,
- the relevant log or traceback.

Small reproducible examples are the fastest to debug.
