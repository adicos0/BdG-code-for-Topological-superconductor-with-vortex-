"""Build source-only reading notebooks: original thesis pages and Wolfram code.

No solver, new derivation, altered parameter set or generated numerical plot.
Run extract_original_code.py first. Page PNGs must match the saved manifest.
"""
from pathlib import Path
import json
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
THESIS = 'https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view'
REPO = 'https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-'
RAW = 'https://raw.githubusercontent.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/main'
PIN = 'd88d12cded1a954d53c2d36d515525c8100ded1c'
SOURCES = {stem: json.loads((ROOT/'original_code'/f'{stem}.json').read_text())
           for stem in ['bdg_topological', 'added_vortex']}


def md(text):
    return nbf.v4.new_markdown_cell(text.strip())


def intro(number, title, description):
    return [md(f'''# {number}. {title}

{description}

**Source-only edition.** The equations, derivation and figures below are reproduced as complete pages from [your 2023 MSc thesis]({THESIS}), with the original notation and captions. Code blocks are transcribed directly from your original Mathematica input cells; their parameters and comments are preserved.

**How to use this notebook:** read the thesis pages, then the matching source blocks. This is a reading notebook, not a Python translation. To execute the calculation, open the original [BdG topological.nb]({REPO}/blob/{PIN}/BdG%20topological.nb) or [added vortex.nb]({REPO}/blob/{PIN}/added%20vortex.nb) in Mathematica. Source fragments below are labelled with their original file and input cell. Full readable transcriptions are in [original_code](../original_code/).

**Figure provenance:** the displayed plots are archived thesis figures, not newly computed outputs. The source files contain commented sections and saved output from earlier evaluations; the notebook does not claim that running the currently active cells reproduces every thesis figure.''')]


def page(number, description):
    name = f'thesis_p{number:03d}.png'
    cell = md(f'''### Thesis page {number} — {description}

![Original thesis page {number}: {description}]({RAW}/reference_pages/{name})

[Open the full-size page]({RAW}/reference_pages/{name}) · [Original PDF, viewer page {number+14}]({THESIS})

*Original page, including its equations, wording, figure captions and interpretation. No equation or curve has been redrawn.*''')
    cell.metadata['thesis_source'] = {'printed_page': number, 'pdf_page': number+14,
                                     'asset': f'reference_pages/{name}'}
    return cell


def excerpt(stem, index, title, start=None, end=None):
    original = SOURCES[stem]['cells'][index-1]
    whole = original['code']
    begin = 0 if start is None else whole.index(start)
    stop = len(whole) if end is None else whole.index(end, begin)
    snippet = whole[begin:stop]
    name = SOURCES[stem]['source_file']
    cell = md(f'''### {title}

**Original code:** `{name}`, input {index}, original `.nb` line {original['line']}. {'Complete input cell.' if begin == 0 and stop == len(whole) else 'Exact contiguous excerpt; see the complete source for its surrounding setup.'}

```mathematica
{snippet}
```''')
    cell.metadata['original_code'] = {'file': name, 'source_key': stem, 'input': index,
                                      'start': begin, 'end': stop,
                                      'source_line': original['line']}
    return cell


def save(filename, cells):
    for index, cell in enumerate(cells):
        cell.id = f'source-{index:03d}'
    nb = nbf.v4.new_notebook(cells=cells)
    # Markdown-only reading notebooks deliberately have no executable kernel.
    nb.metadata = {'source_edition': 'Original Wolfram code and thesis page facsimiles',
                   'original_source_commit': PIN}
    nbf.validate(nb)
    nbf.write(nb, ROOT/'notebooks'/filename)


cells = intro('01', 'Lattice and normal-state Hamiltonian',
              'Follow the tight-binding derivation in §2.7 and the topological-superconductor Hamiltonian in §4.4, then read the lattice construction from the original file.')
cells += [md('## 1. Derivation as printed in the thesis')]
for p, title in [(17, 'start of §2.7'), (18, 'tight-binding derivation'),
                 (19, 'square lattice, Figure 2.1 and dispersion'),
                 (20, 'square-lattice neighbour vectors and dispersion'),
                 (45, '§4.4 and Eq. (4.8)'), (46, 'Eqs. (4.9)–(4.13)')]:
    cells.append(page(p, title))
cells += [md('''## 2. Original lattice code

The following blocks are consecutive parts of input 1 in `BdG topological.nb`. Together they preserve the complete input cell, including `Nx = 20`, `Ny = 20`, the hopping dictionaries, periodic-image lookup and `Hop1x`. `added vortex.nb` has its own complete transcription in `original_code/added_vortex.wl`.

**Source distinction:** the hopping dictionaries in the files and the printed Hamiltonian are not identical. Both are shown as supplied. No replacement Hamiltonian or corrected hopping is inserted here.'''),
          excerpt('bdg_topological', 1, 'Parameters and hopping dictionaries', end='Nx = 20;'),
          excerpt('bdg_topological', 1, 'Site numbering and original lattice size', start='Nx = 20;', end='(*Hopping Matrix*)'),
          excerpt('bdg_topological', 1, 'Periodic hopping, spin mixing and Hop1x', start='(*Hopping Matrix*)')]
save('01_lattice_and_normal_hamiltonian.ipynb', cells)

cells = intro('02', 'BdG matrix and diagonalization',
              'Read the BdG construction in the thesis and the same spin/pairing blocks and Eigensystem call used in the original Mathematica file.')
cells += [md('## 1. Thesis equations and derivation'), page(30, '§4.1 Hamiltonian and mean-field construction'),
          page(31, 'Bogoliubov transformation, eigenproblem and self-consistency equations'),
          page(46, 'pairing and spin-resolved BdG matrix, Eqs. (4.10)–(4.13)'),
          md('''## 2. Original code blocks

The original program defines the onsite pairing matrix in input 2 and constructs and diagonalizes `H` inside input 5. The full iteration surrounding the latter excerpt is reproduced in notebook 03. The code below retains the original matrix signs and basis ordering; it adds no canonical-basis transformation.'''),
          excerpt('bdg_topological', 2, 'Onsite pairing matrix SCs'),
          excerpt('bdg_topological', 5, 'Active BdG matrix and Eigensystem', start='  Huu =', end='  fermi ='),
          md('## 3. Original plotting commands and archived thesis spectrum'),
          excerpt('bdg_topological', 6, 'Original sorting input'),
          excerpt('bdg_topological', 7, 'Original ListPlot command — BdG topological.nb'),
          excerpt('added_vortex', 7, 'Original ListPlot command — added vortex.nb'),
          md('''The plot labels, ranges and grid lines above are retained exactly, including their differences between source files. The thesis figure below is an archived result with its own printed caption. The saved plot labels are not asserted to match the parameters in every currently active input cell.'''),
          page(47, 'original Figure 4.20 and its discussion')]
save('02_bdg_matrix_and_diagonalization.ipynb', cells)

cells = intro('03', 'Self-consistency in the original file',
              'Use the thesis mean-field equations and the original occupation, pairing and chemical-potential iteration. The two source files retain their own parameter values.')
cells += [md('## 1. Original self-consistency equations and explanation'),
          page(31, 'self-consistent quantities and eigenproblem'),
          page(32, 'LDOS expression and calculation discussion'),
          page(33, 'boundary conditions and calculation parameters'),
          md('''## 2. Original parameters and initialization

Both complete parameter cells are shown so the different `vs` and `nexp` values remain visible. Earlier commented experiments remain commented. The original 20 × 20 setup is in notebook 01; it has not been replaced by a smaller tutorial lattice.'''),
          excerpt('bdg_topological', 3, 'Fermi function and initialization — BdG topological.nb'),
          excerpt('added_vortex', 3, 'Fermi function and initialization — added vortex.nb'),
          md('''## 3. Complete original iteration

These are the complete saved input cells, including `iter = 100`, `Eigensystem`, the full occupation/gap sums, mixing, chemical-potential update, commented vortex prototype and file exports. No stopping rule, extra diagnostic calculation or numerical output has been added. Running the original exports in Mathematica can overwrite the named `.wdx` files.'''),
          excerpt('bdg_topological', 5, 'Complete iteration — BdG topological.nb'),
          excerpt('added_vortex', 5, 'Complete iteration — added vortex.nb')]
save('03_self_consistent_pairing.ipynb', cells)

cells = intro('04', 'Vortex profile and original figures',
              'Read §4.2 and its original vortex profile and figures, followed by the topological-superconductor vortex results in §4.4.')
cells += [md('## 1. Vortex expression and explanation from the thesis'),
          page(34, '§4.2 and Eq. (4.2)'), page(35, 'Eqs. (4.3)–(4.5) and explanation'),
          md('''## 2. Vortex code exactly as it appears in the files

The following excerpts include their enclosing `(* ... *)` comment delimiters. The vortex prototype is **commented out** in both original files. Its original `ArcTan`, length expression, position and signs are preserved; it has not been activated or replaced by `atan2`.

**Reading note:** preserving the printed/source expression does not resolve its phase-convention issues. The original expression and the archived figures are shown as source material, without asserting that this commented block currently generates those figures.''')]
for stem in ['bdg_topological', 'added_vortex']:
    # Include the original enclosing comment delimiters, not only its contents.
    full = SOURCES[stem]['cells'][4]['code']
    vpos = full.index('vsite =')
    start = full.rfind('(*', 0, vpos)
    # The next source section is outside the vortex comment block.
    end = full.index('  (*spmtcalcnn', vpos)
    cells.append(excerpt(stem, 5, 'Commented vortex prototype — '+SOURCES[stem]['source_file'],
                         start=full[start:vpos+len('vsite =')], end='  (*spmtcalcnn'))
cells += [md('## 3. Original conventional-superconductor vortex plots')]
for p, title in [(36, 'Figure 4.5 and its discussion'), (37, 'Figure 4.6 and explanation'),
                 (38, 'Figures 4.7–4.8'), (39, 'Figure 4.9 and subsequent impurity definitions')]:
    cells.append(page(p, title))
cells += [md('## 4. Original topological-superconductor vortex plots'),
          page(47, 'Figure 4.20 and its interpretation in the thesis'),
          page(48, 'Figures 4.21–4.22')]
save('04_vortex_profile_and_spectrum.ipynb', cells)

cells = intro('05', 'Impurities and LDOS — original thesis material',
              'Use the impurity definitions, LDOS expression and plots already present in the thesis. The original code is shown only where it actually contains the corresponding terms.')
cells += [md('## 1. LDOS and impurity equations as printed'),
          page(32, 'original LDOS expression'), page(39, 'Eqs. (4.6)–(4.7)'),
          page(40, 'nonmagnetic and magnetic impurity discussion; Figure 4.10'),
          page(46, 'impurity Hamiltonian Eq. (4.13)'),
          md('''## 2. What the original code contains

Input 1 sets `Vimp = 0`. In the active loop, `+imp` is inside a comment. These exact source excerpts are reproduced below. The committed input cells do not contain an active impurity-strength sweep or an active routine generating the thesis LDOS curves, so no substitute sweep or LDOS implementation is supplied.'''),
          excerpt('bdg_topological', 1, 'Original impurity parameter', start='Vimp = 0;', end='empty='),
          excerpt('bdg_topological', 5, 'Original commented impurity term and active spin blocks', start='  Huu =', end='    Huv='),
          md('## 3. Original conventional-superconductor impurity figures')]
for p, title in [(41, 'Figures 4.11–4.12'), (42, 'Figures 4.13–4.14'),
                 (43, 'Figures 4.15–4.16'), (44, 'Figures 4.17–4.18'),
                 (45, 'Figure 4.19; beginning of §4.4')]:
    cells.append(page(p, title))
cells += [md('## 4. Original topological-superconductor impurity figures')]
for p, title in [(49, 'Figures 4.23–4.24'), (50, 'Figures 4.25–4.26 and magnetic-impurity discussion'),
                 (51, 'Figures 4.27–4.28'), (52, 'Figure 4.29')]:
    cells.append(page(p, title))
cells += [md('''## 5. Scope of the original results

All plots above come directly from the thesis. Their captions and discussion are preserved, including the author's physical interpretation. They are not new numerical evidence or an independent validation of that interpretation. The original `.nb` files remain the authoritative record of code and saved Mathematica outputs; a complete executable workflow for every impurity and LDOS figure is not present in the committed inputs.''')]
save('05_impurities_and_ldos.ipynb', cells)
print('Built five source-only reading notebooks.')
