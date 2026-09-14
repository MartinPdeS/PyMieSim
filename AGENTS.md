# Repository instructions

- Format C++ sources with the repository `.clang-format` configuration.
- Use four spaces for indentation; do not use tab characters.
- Keep C++ declarations in headers, implementations in `.cpp` files, and pybind11 bindings in `interface.cpp` files.
- Do not use `from __future__ import annotations`. Quote forward references when runtime evaluation would otherwise fail.
- Prefer keyword arguments in documentation examples when the public API supports them.
- Run relevant pytest tests and `git diff --check` after changes.
