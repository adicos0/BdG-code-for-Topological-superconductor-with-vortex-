"""Verify code excerpts and unchanged diagram pixels without running Wolfram."""
from pathlib import Path
import hashlib
import json
import re
import nbformat
from PIL import Image
from extract_original_code import extract

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((ROOT/'reference_pages/manifest.json').read_text())
    for name, expected in manifest['original_files_sha256'].items():
        assert sha256(ROOT/name) == expected, f'Original changed: {name}'
    sources = {}
    for stem in ['bdg_topological', 'added_vortex']:
        saved = json.loads((ROOT/'original_code'/f'{stem}.json').read_text())
        fresh = extract(ROOT/saved['source_file'])
        assert saved['cells'] == fresh, f'Transcription differs: {stem}'
        sources[stem] = fresh
        text = '\n\n'.join(f'(* Original input {c["index"]}; {saved["source_file"]}, line {c["line"]} *)\n'+c['code'] for c in fresh)+'\n'
        assert (ROOT/'original_code'/f'{stem}.wl').read_text() == text
    for page in manifest['pages']:
        assert sha256(ROOT/page['path']) == page['sha256'], page['path']
    diagrams = json.loads((ROOT/'diagrams/manifest.json').read_text())
    for diagram in diagrams:
        assert sha256(ROOT/diagram['path']) == diagram['sha256']
        with Image.open(ROOT/diagram['source']) as source, Image.open(ROOT/diagram['path']) as shown:
            expected = source.crop(diagram['rectangle'])
            assert shown.size == expected.size
            assert shown.mode == expected.mode
            assert shown.tobytes() == expected.tobytes(), diagram['path']
    assert sha256(ROOT/'original_code/vortex_ldos_original.m') == 'a5c27b03f1d6ea45677a2d47d808d641a5c17e633a25a0965be2ac6adf6ae635'
    code_count, diagram_count, adapted_count = 0, 0, 0
    used_diagrams = set()
    for path in sorted((ROOT/'notebooks').glob('*.ipynb')):
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
        for cell in notebook.cells:
            assert cell.cell_type == 'markdown', 'Source-only edition has no execution outputs'
            if 'original_code' in cell.metadata:
                info = cell.metadata['original_code']
                original = sources[info['source_key']][info['input']-1]['code']
                expected = original[info['start']:info['end']]
                actual = cell.source.split('```mathematica\n', 1)[1].rsplit('\n```', 1)[0]
                assert actual == expected, f'Code changed: {path.name}'
                code_count += 1
            if 'adapted_code' in cell.metadata:
                info = cell.metadata['adapted_code']
                expected = (ROOT/info['path']).read_text()[info['start']:info['end']]
                actual = cell.source.split('```mathematica\n', 1)[1].rsplit('\n```', 1)[0]
                assert actual == expected, f'Adapted code differs: {path.name}'
                adapted_count += 1
            assert not re.search(r'thesis|chapter|reference_pages|drive.google', cell.source, re.I), path.name
            if 'diagram' in cell.metadata:
                asset = cell.metadata['diagram']['asset']
                assert (ROOT/asset).is_file(), asset
                used_diagrams.add(asset)
                diagram_count += 1
            for referenced in re.findall(r'/diagrams/([a-z_]+\.png)', cell.source):
                assert (ROOT/'diagrams'/referenced).is_file()
        print('PASS', path.name)
    assert used_diagrams == {d['path'] for d in diagrams}
    assert adapted_count == 11, adapted_count
    print(f'{code_count} original and {adapted_count} adapted code excerpts; {diagram_count} diagrams with pixel-identical source regions.')
    print('No document comparisons or page scans in notebook text. Original .nb/.wdx hashes unchanged.')


if __name__ == '__main__':
    main()
