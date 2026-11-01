# Repository instructions

- Format C++ sources with the repository `.clang-format` configuration.
- Use four spaces for indentation; do not use tab characters.
- Keep C++ declarations in headers, implementations in `.cpp` files, and pybind11 bindings in `interface.cpp` files.
- Do not use `from __future__ import annotations`. Quote forward references when runtime evaluation would otherwise fail.
- Prefer keyword arguments in documentation examples when the public API supports them.
- Run relevant pytest tests and `git diff --check` after changes.

## Shared Python style preferences

Apply these conventions to new and edited Python code. Leave unrelated existing code unchanged. Preserve each repository's existing public APIs and project-specific scientific or behavioral requirements.

- Use keyword-only parameters by default for new functions and methods, apart from implicit `self` and `cls`. Preserve existing public calling conventions unless an API change is part of the task.
- Use NumPy-style docstrings for classes, functions, methods, and properties. Document parameters, return values, units, array shapes, and limitations where they apply.
- Keep imports at module scope when practical. Use local imports only when needed to avoid circular dependencies or defer optional, expensive, or side-effectful imports.
- Keep all imports in one contiguous block, without blank lines between standard-library, third-party, and project imports.
- Avoid redundant self-alias imports such as `from package.module import Name as Name`; import the name directly unless the alias clarifies its meaning.
- Keep simple, single-purpose UI component calls on one line when practical, even when a text string is long. Use multiline formatting when several arguments or nested components make the call easier to read.
- Prefer descriptive variable and parameter names over abbreviations, including for local values and callback parameters.
- Avoid extracting a helper used only once unless it improves clarity or separates behavior that belongs in another layer.
- Keep callbacks focused on reading inputs, calling a service or workflow, and returning outputs. Move domain decisions and complicated data conversion into the appropriate service or workflow module.
- Write comments to explain non-obvious reasoning, scientific assumptions, or compatibility constraints; do not use them to narrate code that is already clear.
