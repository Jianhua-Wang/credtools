# Repository Guidelines

## Project Structure & Module Organization
- Source: `credtools/` — core library modules and subpackages.
- Tests: `tests/` — mirrors `credtools/` (e.g., `tests/test_utils.py`).
- Docs: user-facing docs in `docs/`; internal notes and manuscript drafts in `@devdocs/`.
- Examples/Tools: `examples/` or `scripts/` if present for runnable snippets.
- Config: `pyproject.toml`, `.ruff.toml`/`.flake8`, and `mypy.ini` (if present).

## Build, Test, and Development Commands
- Install (editable): `python -m pip install -e .[dev]` — local dev with extras.
- Test: `pytest -q` — run unit tests.
- Lint: `ruff check .` — static checks.
- Format: `black .` — apply code formatting.
- Type check: `mypy credtools` — strict typing.
- Build: `python -m build` — create sdist/wheel.

## Coding Style & Naming Conventions
- Python 3; 4‑space indentation; UTF‑8.
- Names: `snake_case` (modules/functions), `PascalCase` (classes), `UPPER_SNAKE_CASE` (constants).
- Public APIs require type hints; prefer explicit return types.
- NumPy-style docstrings for public modules/classes/functions.
- Keep files focused and not too long; factor helpers into small modules.

## Testing Guidelines
- Framework: `pytest` with plain asserts and fixtures in `tests/conftest.py`.
- Naming: files `tests/test_<module>.py`; functions `test_<behavior>()`.
- Add tests for new features and bug fixes; include edge cases.

## Commit & Pull Request Guidelines
- Use Conventional Commits (e.g., `feat:`, `fix:`, `docs:`).
- PRs are small and focused; include description, motivation, and linked issues.
- Update docs, type hints, and tests with code changes; add usage notes if behavior changes.

## Security & Configuration Tips
- Never commit secrets; use environment variables and `.env` (excluded from VCS).
- Validate and sanitize input at module boundaries; prefer safe defaults.
- Pin or constrain dependencies for reproducible builds.

## Agent‑Specific Instructions
- Enforce type linting and resolve warnings quickly.
- Write NumPy-style docstrings for public symbols.
- Avoid long files; refactor when modules grow large.
- Put user docs in `docs/`; internal docs and manuscript/drafts in `@devdocs/`.
