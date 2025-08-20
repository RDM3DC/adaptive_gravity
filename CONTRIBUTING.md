# Contributing

Thanks for helping improve **adaptive_gravity**!

## Dev setup
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```
(Optional) Install pre-commit and enable ruff/black locally.

## Code style
- Run `black .` and `ruff .` before pushing
- Add tests in `tests/` when you add features
- Keep functions documented with short docstrings and units

## PR checklist
- [ ] Includes units in docstrings
- [ ] Adds/updates tests
- [ ] Updates README if flags/CLI changed
- [ ] Mentions SPARC data citation where relevant

By contributing, you agree to license your work under MIT.
