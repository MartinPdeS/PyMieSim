# Documentation logo

The documentation uses a cyan particle, an incident wave, and a blue scattering
lobe as its shared logo and favicon. The theme renders the PyMieSim wordmark as
text so it stays sharp in both light and dark modes.

Assets live in `docs/source/_static/`:

- `logo.png`: transparent 512 × 512 header mark.
- `favicon.png`: transparent 256 × 256 icon.
- `favicon.ico`: 16, 32, 48, 64, 128, and 256 pixel browser icons.
- `gallery-default.svg`: editable layout for gallery examples without plots,
  reusing the logo on a slate-blue background.
- `gallery-default.png`: 600 × 600 browser-rendered export used by Sphinx-Gallery.

The mark was generated with the built-in image-generation tool. PNG resizing
and ICO encoding use the same selected artwork without changing its design.

## Generation prompt

```text
Use case: logo-brand.
Asset type: compact PyMieSim favicon and matching documentation logo mark.
Create a single polished, beautifully balanced scientific-software logo symbol, on a genuinely transparent background. This should look like a professional geometric optical brand identity, with a near-circular silhouette, not a literal diagram. A solid circular cyan particle at the center, surrounded by one broad deep-blue curved scattering lobe that wraps around it, with a single short clean teal sinusoidal beam entering from the left. The beam and scattering lobe should feel like parts of one unified emblem, not three separate symbols. Keep the whole symbol compact and visually centered. Broad smooth shapes, clear open negative space, deliberate consistent curves, flat solid cyan/teal and blue fills matching the cyan accents of the PyMieSim documentation. No glow, gradients, white outlines, texture, speckles, fine lines, atom orbits, wifi/RSS stacks, generic app tile, arrows, text, letters, labels, or presentation mockups. It must be highly recognizable at 16 and 32 pixel favicon sizes, with generous stroke weight and only three main shapes. Fill approximately 85 percent of the square canvas with the near-circular mark, keeping equal transparent padding on all sides.
```
