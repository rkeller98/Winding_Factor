# Building the learning manuscript

The existing `main.tex` is now an eight-chapter English learning manuscript.
Its previous version is retained in Git history. The numerical models, GUI,
and atlas are unchanged.

From this directory, run:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The output is `main.pdf`. No bibliography or external literature is needed.

The source reads eight chapter files. Selected Set 2 contour views are retained
in `figures/conventional_*.png`, with field views exported from the existing plot
routine using the cached conventional factor. Every d-flux view puts q on X and
physical or reduced d on Y, including its error heatmap. The dense vector export
was impractically large, so high-resolution PNGs were retained. No factor
identification or optimization was run for this rewrite.

Numerical evidence comes from the existing contour and derivative caches,
the archived whole-group comparison, and the already verified pointwise,
consistent-global, and four-parameter results from 25 October 2026.
These populations and references are deliberately kept separate in the text.
Recent GUI-only results are transcribed from those completed numerical outputs;
they are not described as additional saved result files.

The controller discussion was checked against two available local
implementations. Their active configuration during the measurement campaigns
remains unverified. The requested handwritten `Motor Model.pdf` was not
found in the searched local project and usual document locations.

The chapter files are authored manuscript sources. Older bulk plotting
exports can overwrite `results.tex` and `recommendation.tex`; do not use
those exports merely to build this document. No atlas build is required.
